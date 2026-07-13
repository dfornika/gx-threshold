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

### Read QC + dehosting

`READ_QC_AND_DEHOSTING` (`subworkflows/local/read_qc_and_dehosting.nf`) bundles
everything that runs on every sample before any species-ID/library-type
analysis touches the reads: platform tagging, FastQC, fastp/fastplong
trimming, and dehosting (`DEHOST`).

FastQC and fastp/fastplong each run **twice** - once on the raw input reads,
and again on the final (trimmed + dehosted) reads that flow to every
downstream analysis, so the QC report visible in MultiQC reflects what
downstream analyses actually see, not just the raw input. Running the same
process twice with different roles in one workflow uses **module aliasing**
(`include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_FINAL } from ...`, and
likewise for fastp/fastplong) - this is the standard nf-core idiom for this
exact situation (e.g. nf-core/rnaseq's `FASTQC_RAW`/`FASTQC_TRIM`), not a
workaround; there isn't a cleaner alternative.

**The final pass is measurement-only by construction, not by convention.**
FastQC never modifies reads, so re-running it on the final reads is
inherently safe - nothing to guard against. fastp/fastplong do modify reads
by design (that's the entire point of the first pass), so the final pass
uses their built-in `discard_trimmed_pass: true` option (a real, tested,
documented upstream nf-core module feature - "use fastp for the output
report only"): with it set, fastp/fastplong never write a trimmed-reads
file at all for that invocation. There is nothing from this pass that could
accidentally be published or picked up downstream, regardless of what
fastp/fastplong would otherwise have changed - no extra `--disable_*` flags
or output-suppression logic needed.

Both passes publish to the same per-tool directory (`${outdir}/fastqc/`,
`${outdir}/fastp/`, `${outdir}/fastplong/`) with an `_raw`/`_final` filename
suffix (`ext.prefix`, `conf/modules.config`) distinguishing them, rather
than separate subdirectories - keeps the existing one-directory-per-tool
layout unchanged for every other stage.

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

### Species composition (pure culture vs. metagenomic)

`SPECIES_COMPOSITION` (`modules/local/species_composition/`,
`bin/classify_species_composition.py`) answers a different question than
species-ID's own top-hit call: not "what is the most likely organism" but
"is there *one* organism here, or many". Sourmash's own `gather` output
already answers this directly - it decomposes a sample into the set of
reference genomes that best explain it, so a pure culture shows one
genome explaining nearly everything, while a metagenomic sample takes many
genomes to explain a similar fraction, none of them dominant. No new
alignment or tool is needed - `SOURMASH_GATHER`'s existing output
(previously truncated to its top row by `parse_species_id.py`) is reused in
full.

The chosen statistic is `f_unique_weighted` (gather's incremental,
non-overlapping contribution of each match - the same number gather's own
human-readable summary shows as `p_query`), summarized as the *effective
number of genomes* (inverse Simpson index over each hit's share of the
total explained fraction) - low (~1) means one genome dominates; high means
many contribute comparably.

**The Escherichia coli / Shigella problem**: a real, reproducible issue found
while validating this (not assumed) - a pure *E. coli* culture's gather
output includes several *Shigella* species as separate "hits" purely because
they are ~98%+ ANI to *E. coli* (a pre-genomic-era clinical naming split, not
a real genomic distinction - the same root cause behind the mixed-call
behaviour widely seen with Kraken2 on this exact genus pair). Naively counted,
this looks like spurious extra breadth for what is actually one organism.

**Fix, applied principled rather than as a hardcoded exception list**:
`REFERENCE_GENOME_DISTANCES` + `CLUSTER_REFERENCE_GENOMES`
(`bin/cluster_reference_genomes.py`) run once per database (not per sample):
all-vs-all `mash dist` on the database's own genomes, converted to an ANI
estimate, union-find clustered at `--species_composition_ani_threshold`
(default 95%, the conventional operational species boundary - the same one
GTDB uses, which is why GTDB's own taxonomy already merges *Shigella* into
*Escherichia coli*). `classify_species_composition.py` then reports **both**
the naive (raw per-accession) and ANI-adjusted (cluster-collapsed) breadth
side by side, rather than only the adjusted figure - specifically so the
adjustment's effect can be spot-checked over time (a `note` column names
which accessions actually got merged for a given sample) instead of trusted
blindly. The verdict is based on the adjusted numbers.

This isn't specific to *E. coli*/*Shigella* - clustering the 100-genome dev
panel (see below) at 95% ANI also merged *Neisseria gonorrhoeae* with
*N. meningitidis* without being told to, another well-documented close-relatedness
pair. A hardcoded synonym list would need to know about cases like this in
advance; the ANI-based approach doesn't.

**Known limitation, found (not assumed) via the deeper dev metagenomic
samples**: breadth is only meaningful if the reference database actually
covers the sample's real organisms. A real animal-gut shotgun sample whose
actual community wasn't well represented in our dev panel produced only one
weak hit (0.02% of k-mers explained) and would otherwise have been called a
confident "pure culture" - wrong, for the right reason: it isn't pure, the
database just has nothing relevant. `--min-explained-frac` (default 5%) gates
this - below it, the verdict is `inconclusive` rather than a confident but
unreliable guess. The two 16S dev samples hit the same gate for a different
reason (see below).

**Validated against the tiny `tests/data/species_db/` panel (all correctly
`pure_culture`, as expected for `test_full`'s four samples) and the larger
100-genome dev panel** (`data/dev_species_db/` - see below): both pure
cultures correctly `pure_culture` (including `KPNEUMONIAE_WGS_ONT_DEV`, whose
ONT error rate dilutes its own top-hit fraction to just ~30%, but the
effective-genome statistic still resolves it correctly since the remaining
mass is spread across many individually-tiny noise hits, not one real
competitor), the one dev metagenomic sample with good database coverage
(`GUT_SHOTGUN_ONT_DEV`, real human gut) correctly `metagenomic`, and the
three poorly-covered samples (an animal-gut shotgun sample plus both 16S
amplicon samples) correctly `inconclusive` rather than a wrong confident
answer.

**16S rRNA amplicon data needs a different reference entirely, not just a
lower confidence bar**: sourmash recovered under 3% of k-mers for both 16S
dev samples against the whole-genome database, because 16S reads sample
~1.5kb of one gene, not whole genomes - there's almost nothing for a
whole-genome sketch to match. This confirms an earlier discussion: assessing
purity of 16S/marker-gene libraries needs a dedicated 16S reference database
(NCBI's curated 16S rRNA RefSeq targeted-loci project, BioProject PRJNA33175,
is the current plan - free, actively maintained, and stays on NCBI taxonomy
like the rest of this pipeline, unlike SILVA which needs a commercial license
beyond academic use), not an extension of the whole-genome species-ID DBs.
Not yet implemented.

**Needs both mash and sourmash enabled** (`--skip_mash`/`--skip_sourmash`
both `false`) - it clusters mash's reference genomes, then applies that to
sourmash's gather output; toggle the whole stage with
`--skip_species_composition`. Unlike the alignment-based library-type
stages, this needs no fetched reference genome, so it runs by default in
every profile including plain `test`. Collected into
`${outdir}/species_id/species_composition_summary.tsv`.

#### Dev-only species database

A moderately-sized, gitignored mash/sourmash/sylph database
(`data/dev_species_db/`, not `tests/data/species_db/`'s tiny 10-genome
panel) - 100 real bacterial genomes spanning gut commensals, clinical
pathogens, close *E. coli*/*Klebsiella* relatives (deliberately included for
this exact contamination/complex-resolution testing), and oral/environmental
diversity. Composition-breadth questions only mean something against a
database with real taxonomic diversity, which the tiny validation panel
was never built to provide. See `data/dev_species_db/README.md` for the
full panel and regeneration recipe.

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

### Species-ID consensus + reference genome fetch/cache

- `SELECT_REFERENCE_ACCESSION` (`bin/select_reference_accession.py`) picks one
  accession per sample from the mash/sourmash/sylph calls: majority vote on
  `species_taxid` (2-of-2 or 2-of-3 agreement; a lone hit is trivially its own
  consensus), falling back to a fixed priority - sylph, then mash, then
  sourmash - if there's no majority, based on which gave the cleanest results
  in our own evaluation (see `tests/data/species_db/README.md`). This is "the"
  species-ID consensus call - it runs unconditionally as part of `SPECIES_ID`
  itself (as long as at least one of mash/sourmash/sylph is enabled), not
  gated behind reference-genome fetch, so a consensus organism call exists
  even when `--skip_reference_genome_fetch` is set. Collected into
  `${outdir}/species_id/species_id_consensus_summary.tsv`.
- `REFERENCE_GENOME` (`subworkflows/local/reference_genome.nf`) consumes that
  consensus accession and fetches the whole-assembly FASTA (all
  contigs/plasmids) via NCBI's E-utilities (`esearch | elink | efetch`,
  resolving assembly → nucleotide sequences), cached by accession using
  Nextflow's `storeDir` in `--reference_genome_cache_dir` - a genome is never
  re-downloaded across runs, and never re-downloaded twice in the same run
  even if multiple samples resolve to the same accession (deduplicated before
  fetching). This is the piece that unblocks the alignment-based library-type
  approaches below for Nanopore, where the fastp-JSON heuristic doesn't apply.

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

**`FETCH_REFERENCE_GENOME` is on by default** (`--skip_reference_genome_fetch`
defaults to `false`) - the planned deployment environments have network
access, so unlike the earlier off-by-default design this is now assumed
available rather than opt-in. A per-accession fetch failure (unreachable
network, a bad/withdrawn accession, transient NCBI rate-limiting) doesn't fail
the whole run: `conf/modules.config` gives this process
`errorStrategy = { task.attempt <= 2 ? 'retry' : 'ignore' }` - it retries
twice, then gives up gracefully. A sample whose fetch is ultimately ignored
just ends up with no reference genome, identical to today's "stage off" path
for that sample (`ALIGNMENT_BASED_LIBRARY_TYPE` simply doesn't run for it,
and the library-type consensus below has fewer votes for that one sample).
The offline `test`/`test_full` profiles explicitly set
`skip_reference_genome_fetch = true` to stay fast/deterministic/network-free,
since - unlike the species-ID databases - this stage does real I/O rather
than using a bundled test database.

Verified manually against `-profile test_full,dev,docker` (with
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

**On by default, off whenever reference-genome fetch is off**
(`--skip_reference_genome_fetch`, default `false`) - this stage consumes that
one's output directly rather than introducing its own flag. Collected into
`${outdir}/library_type/library_type_aligned_summary.tsv`.

### Library type, reference-free (read clustering)

`LIBRARY_TYPE_CLUSTER` (`modules/local/read_overlap/`,
`modules/local/library_type_cluster/`, `bin/classify_library_type_cluster.py`)
is a second, independent take on the same question, evaluated in parallel
rather than as a replacement for the two approaches above (see "Testing
strategy" below for why several approaches are being kept side by side for
now). Instead of aligning to a fetched reference genome, it clusters reads
against **each other**: `READ_OVERLAP` runs an all-vs-all minimap2 self
alignment (`ava-ont` for Nanopore, `sr` for Illumina - only one file per
fragment, R1 for paired-end), and `classify_library_type_cluster.py` unions
reads into clusters wherever an overlap passes an identity/coverage
threshold (0.95/0.8 for Illumina, 0.85/0.7 for Nanopore - minimap2's own
preset tuning handles the different error profiles, so these don't need to
be as finely tuned as they might look). The summary statistic is the
effective number of clusters (inverse Simpson index) as a fraction of total
reads - low means a few big clusters dominate (amplicon); high means reads
mostly don't cluster with anything (shotgun).

Being reference-free is a real advantage: it needs no species-ID/
reference-fetch dependency, so it runs by default in every profile,
including the plain `test` profile, with no network dependency
(`--skip_library_type_cluster` to turn it off). It also sidesteps the
specific reference-genome-repeat-structure confound found while validating
`LIBRARY_TYPE_ALIGNED` above, since there's no reference to have repeats in.

**Validated against all 4 real Illumina/Nanopore samples in
`tests/data/real/`, plus the deeper dev datasets (see below) - all 8/8
correct** on first pass, including `KPNEUMONIAE_WGS_ONT_DEV` (the case that
broke the alignment-based approach).

**Known limitation, found (not assumed) while testing on realistic,
QC'd data rather than raw reads**: running this on `KPNEUMONIAE_WGS_ONT_DEV`
*after* `FASTPLONG` trimming (i.e. as the pipeline actually feeds it, not the
raw download) flips the verdict to `amplicon`. Traced to `FASTPLONG`'s own
quality/adapter trimming, not `DEHOST` (this sample has no host spike-in):
untrimmed, noisy read ends drag pairwise identity below threshold for reads
that come from the same real genome repeat (rRNA operons etc. - the same
repeat structure implicated in the `LIBRARY_TYPE_ALIGNED` caveat); trimmed to
their higher-confidence core, those reads become similar enough to cluster.
This is the *same underlying confound* (bacterial genome repeat content)
independently defeating a second, unrelated detection approach once
realistic post-QC data is used. Testing on raw, untrimmed reads to sidestep
this isn't a real fix - nothing stops a user from feeding the pipeline
already-trimmed reads, so the pipeline's actual input (post-QC) is the
correct thing to validate against regardless. Accepted as a documented
caveat, same as the analogous one above; a more thorough fix (e.g. masking
known repeat regions) is a bigger undertaking than this evaluation warrants
right now.

Collected into `${outdir}/library_type/library_type_cluster.tsv`.

### Library type, aligned read-position pileup

`LIBRARY_TYPE_PILEUP` (`modules/local/align_reads/`,
`modules/local/library_type_pileup/`, `bin/classify_library_type_pileup.py`)
is a third approach, going back to reasoning from first principles about what
actually distinguishes the two library types structurally: amplicon
libraries repeatedly re-sequence the same PCR product, so most reads/
fragments pile up at a small number of essentially fixed *aligned* start/end
coordinates (the primer sites); shotgun libraries fragment DNA close to
randomly, so two independent fragments sharing both a start **and** end
coordinate is rare. This is exactly the same positional signature real
duplicate-marking tools (Picard/`samtools markdup`) use to flag PCR/optical
duplicates - the difference here is that for amplicon libraries, this
"duplication rate" isn't a side-effect to police, it's most of the data.

`ALIGN_READS` aligns to the fetched reference (same dependency as
`LIBRARY_TYPE_ALIGNED`) and writes plain-text SAM (not BAM - kept parseable
by stdlib-only Python, no samtools/pysam needed in that container).
`classify_library_type_pileup.py` groups fragments by aligned coordinates via
**single-linkage chaining** (sort positions, merge consecutive ones within a
tolerance) rather than rounding to a fixed grid - a fixed grid can split one
true pileup in two when its spread straddles a grid boundary, which is
exactly what happened during development with real data (a few bases of
Nanopore primer-trim slop split what should have been one large pileup
across grid cells). For paired-end fragments, both the start and end
coordinate must match (within tolerance) to share a signature; for
single-end long reads the end coordinate uses a much looser tolerance (100bp
vs. 15bp for the start) since it's noisier (variable soft-clipping, no mate
to cross-check against) - but it still matters: requiring only a shared
start let the same bacterial genome repeat structure above masquerade as an
amplicon pileup during development, since reads from different genomic
copies of a real repeat can genuinely start near each other by chance. The
summary statistic is the same effective-number-of-signatures fraction used
by `LIBRARY_TYPE_CLUSTER`.

**Validated against all 4 real Illumina/Nanopore samples plus the deeper dev
datasets - 8/8 correct**, including `KPNEUMONIAE_WGS_ONT_DEV` **after**
`FASTPLONG` trimming (unlike `LIBRARY_TYPE_CLUSTER`, this approach doesn't
regress on the post-QC reads - if anything the verdict comes out more
clearly separated on trimmed reads than raw). This is currently the only one
of the three alignment/clustering-based approaches with no known failure
case against the test data gathered so far.

On by default, off whenever reference-genome fetch is off, same as
`LIBRARY_TYPE_ALIGNED`. Collected into
`${outdir}/library_type/library_type_pileup_summary.tsv`.

### Library-type consensus

`LIBRARY_TYPE_CONSENSUS_ANALYSIS` (`subworkflows/local/library_type_consensus.nf`,
`modules/local/library_type_consensus/`, `bin/library_type_consensus.py`) fuses
the (up to 4) verdicts above into one consensus call per sample, and is what
now sets `meta.library_type` for downstream stages to gate on - previously
only `LIBRARY_TYPE_CLUSTER`'s own verdict was tagged, leaving the other three
methods as comparison-only outputs nobody consumed.

**Rule**: majority vote among whichever methods produced a real verdict
(`amplicon`/`shotgun`) for a given sample - a method that reported
`not_classified`/`inconclusive` doesn't count as a vote either way, the same
filter idea `select_reference_accession.py` already uses for "no hit" tools.
Ties (and the trivial "exactly one method ran" case) are broken by a fixed,
**provisional** fallback priority: `pileup` > `aligned` > `cluster` >
`fastp_json` - based on the empirical picture above: `LIBRARY_TYPE_PILEUP` is
currently the only one of the four with no known failure case;
`LIBRARY_TYPE_ALIGNED` and `LIBRARY_TYPE_CLUSTER` share the same
bacterial-genome-repeat-structure confound; `LIBRARY_TYPE` (fastp JSON) is
Illumina-only and validated on the fewest samples. This ordering is expected
to be revisited once more real validation data
(`data/dev_samples/`, `data/dev_metagenomic_samples/`) shows whether any one
or two methods are consistently correct on their own, rather than needing a
full vote.

Since reference-genome fetch is now on by default (see above),
`LIBRARY_TYPE_ALIGNED`/`LIBRARY_TYPE_PILEUP` typically run too, so this vote
is usually a full 4-way vote rather than the 1-2 votes (`LIBRARY_TYPE_CLUSTER`
alone for Nanopore; plus `LIBRARY_TYPE` for Illumina) available with fetch
off. A sample with zero methods able to classify it (e.g. every method
skipped) gets verdict `no_data` rather than disappearing from the pipeline.
Collected into `${outdir}/library_type/library_type_consensus_summary.tsv`;
the `method` column records `majority(n/total)` (a single vote is always its
own trivial majority), `tie_break:<method>`, or `no_data`.

### Why four approaches are being kept side by side

`LIBRARY_TYPE`, `LIBRARY_TYPE_ALIGNED`, `LIBRARY_TYPE_CLUSTER`, and
`LIBRARY_TYPE_PILEUP` are deliberately being accumulated in the pipeline
rather than replacing each other as each new one is built - the shared test
data so far (4 real fixtures + 4 dev samples) is small enough that "8/8
correct" doesn't yet distinguish a robust method from one that's simply
gotten lucky, especially given two of the four already turned out to have a
real, reproducible failure mode on the *same* difficult sample. The plan is
to run all of them against substantially more real data before deciding
whether to keep one, keep several as cross-checks, or retire any.
`LIBRARY_TYPE_CONSENSUS_ANALYSIS` above means `meta.library_type` is already
populated pipeline-wide from a real vote in the meantime, without needing to
prune any method first.

#### Dev-only sample datasets

Four deeper, more representative real samples than the tiny fixtures
above - not committed (`data/dev_samples/`, gitignored) - were used to
develop and validate these alignment/clustering-based stages, one per
platform × library-type combination, including the first real Nanopore
bacterial WGS shotgun sample used anywhere in this project
(`KPNEUMONIAE_WGS_ONT_DEV`, the sample that surfaced the repeat-region
caveats above). See `data/dev_samples/README.md` for the SRA accessions and
regeneration recipe (same `download_fastq` + `rasusa reads -n <N> -s 42`
approach as `tests/data/real/`).

A second, similarly gitignored set (`data/dev_metagenomic_samples/`) covers
real shotgun-metagenomic and 16S rRNA amplicon samples (Illumina + Nanopore
each) - used for the species-composition work above. See
`data/dev_metagenomic_samples/README.md`.

### 16S rRNA amplicon detection (gated, alignment-based)

`SIXTEEN_S_DETECTION` (`subworkflows/local/sixteen_s_detection.nf`,
`modules/local/subsample_reads_head/`, `modules/local/classify_16s_amplicon/`,
`bin/classify_16s_amplicon.py`) answers the specific question the species
composition section above left open: for a sample the whole-genome
species-ID databases can't explain (`meta.composition == 'inconclusive'`),
is that because it's 16S/marker-gene data, or something else entirely
(insufficient database coverage, a novel organism, etc)?

**Gating**: this stage only runs for reads whose `meta.composition` was
already tagged `inconclusive` by `SPECIES_COMPOSITION_ANALYSIS` earlier in
the pipeline - a sample confidently called `pure_culture` or `metagenomic`
from whole-genome k-mer content doesn't need a 16S-specific check. An
earlier design also required `meta.library_type == 'amplicon'` before
running this check, but that was dropped: `LIBRARY_TYPE_CLUSTER` got a real
16S Nanopore sample's library-type call wrong during validation, so gating
on it would have skipped the 16S check for a sample that actually needed it.
Gating on `composition` alone is more conservative and doesn't depend on a
different classifier being right first.

**Why alignment, not k-mer/minhash, specifically for 16S**: the same
sourmash/mash machinery used for whole-genome species-ID was tried first
against a dedicated 16S database, at two different k-mer sizes (k=21,
scaled=1000 and k=31, scaled=200) - both recovered very little (11.2%/3.7%)
and the organisms that did match were implausible (marine/extremophile taxa
for a real gut sample). This isn't a tuning problem: 16S carries long,
near-universally-conserved regions (the property that makes it useful as a
universal marker in the first place), so short k-mers spuriously match many
unrelated references regardless of size. Full alignment doesn't have this
problem - identity and coverage are computed over the whole read, not a
k-mer at a time - so this stage aligns a read subsample to a dedicated 16S
database with `ALIGN_READS` (the same module `LIBRARY_TYPE_PILEUP` uses,
reused as-is since it already accepts any reference FASTA) and computes,
per read, identity from the SAM `NM:i:` tag over the CIGAR alignment-block
length, and coverage as reference-consumed length over read length.
`classify_16s_amplicon.py` reports the fraction of reads clearing both bars
(default 90%/80% identity, 80%/70% coverage for Illumina/Nanopore) against
`--min-passed-frac` (default 50%) for the verdict. Reads with zero
alignments at all count in the denominator, not just mapped ones - a read
that aligns to nothing in the 16S database is itself evidence against 16S,
not an ambiguous non-answer.

`SUBSAMPLE_READS_HEAD` takes the first N reads (`--sixteen_s_max_reads`,
default 500) with plain `head`, not a random subsample - a coarse yes/no
check doesn't need one, and this avoids adding a new dependency just for
subsampling.

**Reference database**: NCBI's 16S ribosomal RNA (Bacteria and Archaea type
strains) BLAST DB (BioProject PRJNA33175) - confirmed small (68MB
compressed, 27,648 sequences) and actively maintained before committing to
it. Not bundled with the pipeline (`--sixteen_s_db`, a FASTA); a dev-only
copy extracted via `blastdbcmd` lives at `data/dev_16s_db/` (gitignored,
~270MB uncompressed) with a `README.md` covering provenance and the
k-mer/alignment finding above. Species-level taxonomy was deliberately not
resolved for all 27,648 entries - this check only needs "16S or not", not
which organism, and can be revisited if the database is ever used for
anything more granular.

**Validated against all 4 real dev metagenomic samples**, both standalone
(direct `minimap2`+script, full read set) and through the full pipeline
(subsampled, post-fastp/dehost reads):

| Sample | Expected | Standalone passed_frac | Pipeline passed_frac | Verdict |
|---|---|---|---|---|
| `GUT_16S_ONT_DEV` | 16S | 88.0% | 95.6% | `16S_amplicon` (correct) |
| `GUT_16S_ILLUMINA_DEV` | 16S | 90.6% | 51.6% | `16S_amplicon` (correct, but a much thinner margin - see below) |
| `GUT_SHOTGUN_ILLUMINA_DEV` | not 16S | 0% | 0% | `other` (correct) |
| `GUT_SHOTGUN_ONT_DEV` | not 16S | 0% | n/a | not reached - `SPECIES_COMPOSITION_ANALYSIS` correctly called this sample `metagenomic` (real database coverage), so the composition gate skipped the 16S check entirely |

`GUT_16S_ILLUMINA_DEV`'s pipeline-run margin (51.6%) is notably thinner than
its standalone one (90.6%) - both clear the 50% bar and land on the correct
verdict, but the gap is large enough to flag rather than ignore. The
standalone check ran the sample's full untrimmed read set directly; the
pipeline run feeds `SUBSAMPLE_READS_HEAD` the first 500 read pairs *after*
`FASTP` and `DEHOST`, which may simply differ in composition from the full
set. Worth re-checking once more real 16S samples are available, rather
than concluding anything from a single data point.

Off by default (`--skip_sixteen_s_detection`, default `true`) - like
reference-genome fetch, it needs an external resource (`--sixteen_s_db`) not
bundled in any test profile. Collected into
`${outdir}/library_type/sixteen_s_detection_summary.tsv`. On a positive or
negative call, the sample's `meta.sixteen_s` tag is set to the verdict for
any downstream stage to use.

### Sample summary

`SAMPLE_SUMMARY` (`modules/local/sample_summary/`,
`bin/build_sample_summary.py`) builds one CSV with one row per sample -
`${outdir}/sample_summary.csv` - so a reviewer can see the main verdict/metric
from every stage above without opening eight different TSVs. It doesn't
replace those TSVs (which keep every field for each stage); it's a single
at-a-glance table across stages, deliberately kept to one or two headline
columns per stage (verdict plus one supporting metric) rather than every
field - see the header row in `bin/build_sample_summary.py` for the exact
columns.

The sample manifest emitted by `READ_QC_AND_DEHOSTING` (sample, platform) is
the only required input and anchors which rows exist; every other stage's
summary TSV is optional. If a stage was skipped (by flag, e.g.
`--skip_reference_genome_fetch`) or simply produced nothing for this run
(e.g. no sample needed the 16S check because none had inconclusive
composition), its columns are just `NA` for the affected rows - the row
itself is never dropped.

**A real, non-obvious bug found while building this**: `collectFile()`
emits nothing at all - not even the seed/header row - when its input channel
is completely empty, rather than falling back to a header-only file. This
first surfaced as `SAMPLE_SUMMARY` silently running **zero times** (no
error, just absent from the process list) whenever any upstream stage's
collated TSV happened to have zero rows for a given run - which is exactly
what happens under the plain `test` profile, where sourmash gather
correctly finds no hit for the tiny synthetic fixtures against the real
species DB (expected behaviour, not a bug in that stage). Every conditional
summary channel now ends in `.ifEmpty([])`, falling back to the same "stage
produced nothing" placeholder used when a stage is toggled off - documented
inline at each site (see `subworkflows/local/species_id.nf` for the
first occurrence). This is the third time in this project a channel has
gone silently empty rather than erroring (see the whole-meta-map-join
footgun in the "Species-ID consensus + reference genome fetch/cache" section
above, and the earlier alignment-based library-type work) - worth remembering
as a recurring
Nextflow failure mode: **missing output shows up as "nothing happened," not
an error**.

**Validated** against both `test` (2 synthetic samples, only a few stages
producing real data) and `test_full` (4 real samples, every column
populated and matching known ground truth - e.g. `ECOLI_WGS` correctly
`pure_culture` with all three species-ID tools agreeing on *Escherichia
coli*).

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
