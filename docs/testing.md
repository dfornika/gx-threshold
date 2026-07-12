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
