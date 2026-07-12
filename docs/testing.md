# BCCDC-PHL/gx-threshold: Testing

This page describes the test datasets, how they are generated, the overall testing
strategy, and how to run the tests. It intentionally stays high level — the specific
processes and their assertions will evolve, so look at the test files themselves
(`*.nf.test`) for current, tool-specific detail.

## Philosophy

Tests should be **fast, deterministic, and runnable offline on a laptop**. To that
end the default test suite uses tiny, purpose-built fixtures rather than real
sequencing runs, so a change either clearly passes or clearly fails without waiting
on large downloads or worrying about drifting inputs.

## Test datasets

The default fixtures live in-repo under [`tests/data/`](../tests/data) and total only
a few hundred kilobytes:

```
tests/data/
  generate.py                 # deterministic generator (see below)
  samplesheet.csv             # test-profile input: one Illumina PE + one Nanopore SE sample
  references/
    host_mini.fasta           # small "host" reference
    microbe_mini.fasta        # small, unrelated "microbe" reference
  reads/
    illumina/  SAMPLE_PE_R1/R2.fastq.gz
    nanopore/  SAMPLE_ONT.fastq.gz
  real/                        # test_full profile input - see below
    README.md                  # provenance: source SRA accessions, subsampling, spike-in
    spike_human.py              # deterministic generator: reads_src + human_mt.fasta -> reads/
    samplesheet.csv
    references/
      human_mt.fasta            # real human mitochondrial genome (rCRS, NC_012920.1)
    reads/                       # reads_src + spiked-in human reads - actual pipeline input, committed (~1.7MB)
      illumina/  KPNEUMONIAE_WGS_R1/R2, ECOLI_WGS_R1/R2, SARS2_AMPLICON_ILLUMINA_R1/R2
      nanopore/  SARS2_AMPLICON_ONT.fastq.gz
    reads_src/                   # pristine downloaded+downsampled SRA reads - NOT committed, see README
  species_db/                  # tiny mash/sourmash/sylph DBs - see tests/data/species_db/README.md
    manifest.csv                 # 10 reference genomes: the 3 known species above + 7 decoys
    build_dbs.py                 # regenerates everything below from manifest.csv
    mash/species_mini.msh
    sourmash/species_mini.sig.zip, sourmash/taxonomy.csv
    sylph/species_mini.syldb
    genomes_src/                 # downloaded genome FASTAs - NOT committed, see README
```

The two references are independent pseudo-random sequences, so reads derived from one
do **not** map to the other. Each read set is a deliberate **mix of host- and
microbe-derived reads**, which lets host-removal tests assert that the host reads are
actually removed and the microbe reads are retained — not just that the step ran.

### How the fixtures are generated

All fixtures are produced by a single seeded, dependency-free script,
[`tests/data/generate.py`](../tests/data/generate.py). Regenerate them at any time
with:

```bash
python3 tests/data/generate.py
```

The script is fully deterministic (fixed seed, timestamp-free gzip), so regenerating
produces byte-identical files. On run it prints a short **manifest** of how many
host- vs microbe-derived reads it wrote for each platform; that manifest is the
ground truth the module tests assert against.

### Synthetic vs. real data

The default fixtures are **synthetic** because they give exact, reproducible control
over what each test contains (e.g. a known host fraction). These are complemented by
a small set of **real, downsampled** reads under
[`tests/data/real/`](../tests/data/real) (see that directory's `README.md` for
provenance): real bacterial WGS (*Klebsiella pneumoniae*, *E. coli*) and SARS-CoV-2
tiling-amplicon runs, pulled from SRA/ENA and downsampled to a few thousand reads each
with [`rasusa`](https://github.com/mbhall88/rasusa) - with a small number of reads
simulated from the real human mitochondrial genome (rCRS, NC_012920.1) spiked into
each sample, so `test_full` exercises real host-read removal too (a few percent host
fraction expected), not just a 0%-host check. As the real-data collection grows
further, test inputs may move to versioned URLs (hosted in our own test-datasets
repository) referenced from the test configs, mirroring the nf-core convention,
rather than being committed in-repo.

### Species-ID databases

[`tests/data/species_db/`](../tests/data/species_db) holds tiny mash/sourmash/sylph
databases built from 10 real reference genomes - the true species of every
`test_full` sample, plus 7 decoys (including two deliberately "hard" same-genus/family
cases) - so a species-ID call in `test_full` is checked against a known right answer
rather than a database with only one possible entry. All three tools - mash,
sourmash, sylph - are wired into the pipeline now (`--mash_db`/`--sourmash_db`/
`--sylph_db`, all on by default) so results can be compared side by side; each
tool's call gets normalised (`modules/local/species_id_summary/`,
`bin/parse_species_id.py`) into one row per sample/tool in
`${outdir}/species_id/species_id_summary.tsv`. See that directory's `README.md`
for the genome panel, how each tool performed, and the summary format.

### Library type (amplicon vs. shotgun)

`LIBRARY_TYPE` (`modules/local/library_type/`, `bin/classify_library_type.py`)
classifies each Illumina paired-end sample as amplicon or shotgun from fastp's own
QC JSON - no new tool or reference database needed. The signal: amplicon libraries
show a sharp, narrow peak in the insert-size histogram (reads cluster around the
fixed amplicon length) and elevated PCR duplication, vs. shotgun's broad insert-size
spread and near-zero duplication. Validated against the three real Illumina samples
in `tests/data/real/`: `SARS2_AMPLICON_ILLUMINA` (56-60% of read pairs within 10bp of
the insert-size peak, 30% duplication) vs. `ECOLI_WGS`/`KPNEUMONIAE_WGS` (5-8%
concentration, 0% duplication) - a wide margin either side of the 20% classification
threshold.

**Nanopore/long-read is not classified** - `fastplong`'s JSON reports neither a
duplication rate nor an insert-size histogram (confirmed by inspecting real output,
not assumed), and the raw read-length distribution on the one real ONT sample we have
didn't show a comparably clean signal either (broad, smoothly decreasing, not a tight
single mode - and there's no real ONT shotgun sample to contrast against). Real
long-read amplicon detection likely needs an alignment-based approach (read start/end
position clustering against a reference); out of scope for now. Long-read samples
report `verdict=not_classified` with a reason, rather than a guess. Toggle the whole
stage with `--skip_library_type`.

### Reference genome selection + fetch/cache

Foundational piece for a future alignment-based approach to the Nanopore
library-type gap above (mapping reads to a reference to check coverage
evenness/read-position clustering, rather than relying on QC-JSON fields that
don't exist for long reads). Two new local modules:

- `SELECT_REFERENCE_ACCESSION` (`bin/select_reference_accession.py`) picks one
  accession per sample from the mash/sourmash/sylph calls: majority vote on
  `species_taxid` (2-of-2 or 2-of-3 agreement; a lone hit is trivially its own
  consensus), falling back to a fixed priority - sylph, then mash, then
  sourmash - if there's no majority, based on which gave the cleanest results
  in our own evaluation (see `tests/data/species_db/README.md`). Collected into
  `${outdir}/reference_genome/reference_selection_summary.tsv`.
- `FETCH_REFERENCE_GENOME` fetches the whole-assembly FASTA (all
  contigs/plasmids) for the selected accession via NCBI's E-utilities
  (`esearch | elink | efetch`, resolving assembly → nucleotide sequences),
  cached by accession using Nextflow's `storeDir` in
  `--reference_genome_cache_dir` - a genome is never re-downloaded across
  runs, and never re-downloaded twice in the same run even if multiple samples
  resolve to the same accession (deduplicated before fetching).

**Tool choice, checked rather than assumed**: `ncbi-genome-download` (what the
nf-core `ncbigenomedownload` module wraps) currently fails against NCBI's live
`assembly_summary.txt` - a real, maintainer-closed bug
([kblin/ncbi-genome-download#237](https://github.com/kblin/ncbi-genome-download/issues/237),
closed as "not a bug in ncbi-genome-download... I don't know how I'd work
around it" - NCBI added a field that itself contains unescaped tabs, breaking
the tab-delimited parser). NCBI's own `datasets` CLI would be the obvious
alternative, but every container build we could find (bioconda/biocontainers,
`staphb/ncbi-datasets`) is too old to work against the current API
("No assemblies found that match selection" even for accessions that exist).
`entrez-direct` is what actually works, using the same E-utilities protocol
NCBI has kept backwards-compatible for decades - the same tool nf-core's own
`entrezdirect` modules already use.

**This stage is off by default** (`--skip_reference_genome_fetch` defaults to
`true`) - unlike the species-ID databases, it does real network I/O rather
than using a bundled test database, so it isn't exercised by the standard test
profiles. Verified manually against `-profile test_full,dev,docker` (with
`--skip_reference_genome_fetch=false --reference_genome_cache_dir <dir>`, via
`-params-file` since Nextflow's CLI doesn't reliably coerce
`--flag=false`/`--flag false` to boolean for schema-validated params): all 4
samples correctly resolve with full 3/3 consensus, exactly 3 (not 4) genomes
get fetched since the two SARS-CoV-2 samples share one accession, and a second
run against the same cache directory skips `FETCH_REFERENCE_GENOME` entirely
(`[skipping] Stored process`) rather than re-fetching.

### Library type, platform-unified (alignment-based)

`LIBRARY_TYPE_ALIGNED` (`modules/local/library_type_aligned/`,
`bin/classify_library_type_aligned.awk`) closes the Nanopore gap left by
`LIBRARY_TYPE` above: it classifies **both** platforms by aligning reads to
the reference genome fetched by the stage above and tracking a streaming
**index of dispersion (variance/mean) of per-base alignment depth**.
Amplicon libraries repeatedly re-sequence the same small set of tiling
primer targets, producing wildly uneven depth; shotgun libraries sample
depth close to uniformly (a Poisson process, index of dispersion ~1).

The statistic updates as reads stream in (an awk script consuming SAM
records piped directly from `minimap2`, so nothing is written to disk or
re-read), and a confident amplicon verdict stops the alignment early rather
than waiting for the whole input - `minimap2` gets `SIGPIPE` when the awk
script exits, which is why the module scopes `set +o pipefail` around just
that one pipe (the pipeline-wide `pipefail` default elsewhere is untouched).
Stopping is one-sided: shotgun requires exhausting the read stream (or a
`max_reads` cap) to confirm the *absence* of the amplicon signal, matching
the observed asymmetry between the two - amplicon gives a fast, distinctive
signal; shotgun is "no signal", which needs more data to be sure of.

**This design point wasn't the first one tried.** An earlier version
reformulated the same statistic as a per-read Bernoulli sequential
probability ratio test ("did this read hit an already-covered reference bin,
or a new one?"), specifically so it could run in `awk` without needing
gamma/special functions. That discretization didn't survive contact with
real data, once tested against deeper, more representative samples than the
tiny committed fixtures (see "Dev-only sample datasets" below): fixing a
false amplicon call on real long-read bacterial shotgun data (chance bin
collisions early in the stream, a birthday-paradox effect) required
coarsening the bins, which then broke detection on real SARS-CoV-2 amplicon
tiling data (whose ~90 distinct amplicons span most of a coarse bin grid, so
"new bin" events - genuine tiling behaviour - looked like evidence *against*
amplicon). Tracking actual depth directly avoids that tension.

**Validated against all 4 real Illumina/Nanopore samples in
`tests/data/real/`**, run through the full pipeline
(`-profile test_full,dev,docker`, `--skip_reference_genome_fetch=false`):

| Sample | Platform | Verdict | Index of dispersion | Method |
|---|---|---|---|---|
| `ECOLI_WGS` | Illumina | shotgun | 1.01 | `eof` |
| `KPNEUMONIAE_WGS` | Illumina | shotgun | 1.01 | `eof` |
| `SARS2_AMPLICON_ILLUMINA` | Illumina | amplicon | 2.51 | `eof`/`threshold_stop` |
| `SARS2_AMPLICON_ONT` | **Nanopore** | **amplicon** | 2.5 | `threshold_stop` |

`SARS2_AMPLICON_ONT` is the win this whole stage exists for - the first
confident, non-`not_classified` verdict on a real Nanopore sample anywhere in
this pipeline.

**Known limitation, found (not assumed) via the deeper dev datasets**: real
bacterial WGS shotgun data at realistic low/modest depth can show index of
dispersion elevated into the amplicon range, purely from genome repeat
structure (rRNA operons, IS elements - common in most bacterial genomes)
concentrating coverage onto one representative copy while the rest of the
genome stays sparse - not an artifact of this implementation, since it
reproduces identically with the plain batch `samtools depth -a`
variance/mean on the same alignment. Confirmed on our real ONT bacterial
shotgun dev sample (`KPNEUMONIAE_WGS_ONT_DEV`, ~1.5x mean depth): index of
dispersion 2.5-4.8 depending on how much of the stream is used, solidly
inside the amplicon range, despite being genuine shotgun sequencing. Not
fixable by threshold tuning alone (genuine amplicon and this artifact
overlap). Accepted as a documented caveat for now, same honest-caveats spirit
as the rest of this project - worth revisiting (e.g. masking known repeat
regions before computing the statistic) if it proves to be a problem in
practice.

**Off whenever reference-genome fetch is off** (`--skip_reference_genome_fetch`,
default `true`) - this stage consumes that one's output directly rather than
introducing its own flag. Collected into
`${outdir}/library_type/library_type_aligned_summary.tsv`.

#### Dev-only sample datasets

Four deeper, more representative real samples than the tiny fixtures
above - not committed (`data/` is gitignored) - were used to develop and
validate this stage, one per platform × library-type combination, including
the first real Nanopore bacterial WGS shotgun sample used anywhere in this
project (`KPNEUMONIAE_WGS_ONT_DEV`, the case that surfaced the repeat-region
caveat above). See `data/dev_samples/README.md` for the SRA accessions and
regeneration recipe (same `download_fastq` + `rasusa reads -n <N> -s 42`
approach as `tests/data/real/`).

## Testing strategy

Tests are layered:

1. **Module tests** (`modules/**/tests/*.nf.test`) — exercise a single process in
   isolation with the tiny fixtures. These favour **explicit assertions** (exact read
   counts, invariants) for correctness, plus snapshots to catch unintended drift.
2. **End-to-end / smoke test** — run the whole pipeline on the `test` profile to
   confirm everything wires together and produces the expected outputs, fully offline.
3. **Realistic test** (`test_full`) — real, downsampled SRA reads (see
   [`tests/data/real/README.md`](../tests/data/real/README.md)) for a more
   representative check than the synthetic fixtures, still small and fast enough to
   run routinely rather than reserved for release validation only.

Continuous integration runs the suite across container engines in parallel — **Docker,
Conda, and Singularity/Apptainer** — via `.github/workflows/nf-test.yml`, so a change
is validated on every supported backend.

### A note on snapshots

Some tests use nf-test **snapshots**. A snapshot records a fingerprint of a test's
outputs (file checksums, channel contents, reported tool versions) the first time it
runs; later runs compare against it. A snapshot therefore guards against **unintended
change**, not correctness in the abstract — when a change is intentional (new output,
deliberate version bump), re-bless the baseline with `--update-snapshot` and commit the
updated `.snap` file.

## Running the tests

### End-to-end (pipeline) run

```bash
# Offline, uses the in-repo synthetic fixtures
nextflow run . -profile test,docker --outdir results
```

Swap `docker` for `conda`, `singularity`, or `apptainer` to check another engine.

For the realistic dataset (real, downsampled SRA reads - still fully in-repo, no
download needed at run time):

```bash
nextflow run . -profile test_full,docker --outdir results
```

### On a memory-constrained machine (e.g. a laptop)

Compose the `dev` profile, which lowers the resource ceilings so every process fits:

```bash
nextflow run . -profile test,dev,docker --outdir results
```

### Component tests with nf-test

[Install nf-test](https://www.nf-test.com/installation/), then from the repo root:

```bash
# a single module
nf-test test modules/local/dehost/tests/main.nf.test --profile=+docker

# everything
nf-test test --profile=+docker
```

`--profile=+docker` adds `docker` to the `test` profile already configured in
`nf-test.config`. Add `,dev` (e.g. `--profile=+docker,dev`) on a constrained machine.
Use `--update-snapshot` when you intend to update snapshot baselines.
