# Species-ID test databases (mash / sourmash / sylph)

Tiny reference databases for evaluating species-ID tools against the known-species
samples in [`tests/data/real/`](../real) (`KPNEUMONIAE_WGS` → *K. pneumoniae*,
`ECOLI_WGS` → *E. coli*, `SARS2_AMPLICON_*` → SARS-CoV-2).

```
species_db/
  manifest.csv       accession, organism, strain, category, taxonomy for all 10 genomes
  build_dbs.py       regenerates everything below from manifest.csv
  mash/species_mini.msh
  sourmash/species_mini.sig.zip, sourmash/taxonomy.csv
  sylph/species_mini.syldb
  genomes_src/       *not committed* (gitignored) - downloaded genome FASTAs
```

## Panel

10 real RefSeq/GenBank complete genome **assemblies**: the 3 organisms actually
present in `tests/data/real/` samples, plus 7 decoys so a correct call is a real
discrimination test rather than "the only entry in the database":

| Accession | Organism | Strain | Role |
|---|---|---|---|
| `GCF_000005845.2` | *Escherichia coli* | K-12 MG1655 | known (`ECOLI_WGS`) |
| `GCF_058435815.1` | *Klebsiella pneumoniae* | Kp0179 | known (`KPNEUMONIAE_WGS`) |
| `GCF_009858895.2` | SARS-CoV-2 | Wuhan-Hu-1 | known (`SARS2_AMPLICON_*`) |
| `GCF_002310435.1` | *Staphylococcus aureus* | Newman | decoy, easy (distant genus) |
| `GCF_051027675.1` | *Pseudomonas aeruginosa* | PAE981 | decoy, easy (distant genus) |
| `GCF_058160315.1` | *Listeria monocytogenes* | 19L270 | decoy, easy (distant genus) |
| `GCF_030295805.1` | *Enterococcus faecalis* | L6D | decoy, easy (distant genus) |
| `GCF_000858645.1` | Human adenovirus type 1 | - | decoy, easy (distant, viral) |
| `GCF_058745625.1` | *Salmonella enterica* | Kentucky D2054 | decoy, **hard** (same family, Enterobacteriaceae, as *E. coli*/*K. pneumoniae*) |
| `GCF_040267715.1` | *Klebsiella oxytoca* | KN | decoy, **hard** (same genus as *K. pneumoniae*) |

### Why assembly accessions, not nucleotide accessions

These are **assembly accessions** (`GCF_...`), matching the identifier convention
used by the real, official/canonical databases for all three tools - not plain
nucleotide accessions (`NC_.../NZ_...`), which is what an earlier version of this
panel used. Checked directly rather than assumed: the official Mash RefSeq sketch
(`test_input/refs/RefSeqSketches_235.msh.gz`) IDs entries like
`Nostoc_azollae_GCF_000196515.1`; sourmash's GTDB-based prepared databases name
signatures like `"GCA_000398885.1 Escherichia coli KTE33..."` ([sourmash-bio/sourmash#3006](https://github.com/sourmash-bio/sourmash/issues/3006));
sylph's pre-built GTDB databases report genome identifiers like `GCA_000011.fasta`
([sylph pre-built databases docs](https://sylph-docs.github.io/pre%E2%80%90built-databases/)).
Building this panel the same way means it actually exercises the identifier format
production will see if/when the real databases are swapped in, rather than a format
specific to this toy panel.

Practically, this also means each entry here is a **whole assembly** (all
chromosomes/plasmids), not a single representative record - several of these
organisms have plasmids that a single-nucleotide-record panel would have missed
entirely (e.g. `GCF_058435815.1` has 3 sequences, `GCF_058745625.1` has 4).

`manifest.csv`'s `source_nucleotide_accession` column records the single
representative RefSeq record originally used to pick each organism (resolved to
its assembly via `elink`, `nuccore` → `assembly`) - kept for provenance, not used
by anything downstream.

### Species-level taxonomy

`manifest.csv` also carries `taxid`, `species_taxid`, and `species_name`, resolved
per-assembly from NCBI's own assembly esummary record - which already reports both
the assembly's own taxid (often strain-level) and its species-level parent, so no
manual taxonomy-tree walking is needed. This is more robust than looking up a taxid
from the organism name string: current ICTV virus species names don't match the
informal names historically used in sequence records, and our own panel demonstrates
this directly - `species_name` for the SARS-CoV-2 genome is `Betacoronavirus
pandemicum` (the current, correct species per ICTV's 2023 binomial renaming), not
"Severe acute respiratory syndrome coronavirus 2" (the `organism` column's informal
name, which is what the sample and sequence records still call it). Same story for
the adenovirus (`species_name` is `Mastadenovirus caesari`, not "Human adenovirus 1").
For most bacteria in this panel, `taxid` and `species_taxid` are identical (no
distinct strain-level node was registered for that particular assembly) - there's no
way to know that in advance per-organism, so `build_dbs.py` always resolves and
records both.

## Building

Requires `ncbi-client` (see `~/code/ncbi-client-py`, or `pip install -e .` from that repo)
and the `mash`, `sourmash`, `sylph` CLIs on `PATH` (a `species-db` conda env with
`conda create -n species-db -c conda-forge -c bioconda python=3.11 mash sourmash sylph`
works - **channel order matters**: `conda-forge` before `bioconda`, or `sourmash`'s
`screed` dependency resolves to an ancient bioconda build and the solve fails).

```bash
python3 tests/data/species_db/build_dbs.py
```

Downloads each assembly's whole-genome FASTA (via `ncbi_client`'s Datasets-API-backed
`download_genome`) into `genomes_src/` (skips ones already present), then (re)builds
all three DB artifacts from scratch. Not byte-deterministic like
`tests/data/generate.py` (depends on live NCBI availability + tool versions), but
idempotent in the sense that re-running with the same manifest and genomes produces
equivalent databases. To add a new accession, resolve its assembly accession first
(`esearch`/`elink`/`esummary`, `nuccore` → `assembly`, as above) rather than adding a
nucleotide accession directly.

## Validated against real data

Spot-checked all three tools/DBs against the real `ECOLI_WGS` and
`SARS2_AMPLICON_ILLUMINA` read sets (see `tests/data/real/`) before wiring anything
into the pipeline:

- **mash** (`mash dist species_mini.msh reads.fastq.gz`): correctly picks the true
  organism as the lowest-distance hit for both samples (e.g. *E. coli* at distance
  0.14 vs. the next-closest decoy, *Salmonella*, at 0.26 - the "hard" decoy is
  measurably closer than the "easy" ones, as expected, but still clearly distinguishable).
- **sylph** (`sylph sketch` + `sylph profile`): correctly identifies both samples
  (98.35% ANI for *E. coli*, 99.66% for SARS-CoV-2) even at the low coverage these
  tiny read sets provide - it's designed for exactly this.
- **sourmash** (`sourmash sketch dna` + `sourmash gather`): correctly identifies
  *E. coli*, but **found nothing for the SARS-CoV-2 sample** with default settings -
  `gather`'s default `--threshold-bp` is 50kb, larger than the entire 30kb SARS-CoV-2
  genome. **Update**: fixed by setting `--threshold-bp 1000` for `SOURMASH_GATHER`
  (see `conf/modules.config`) - confirmed working for all 4 `test_full` samples,
  including both SARS-CoV-2 ones, once wired into the pipeline.

## Wired into the pipeline

All three are now wired in and toggleable independently, so results can be compared
side by side per sample:

- **mash**: `MASH_DIST`, toggle with `--skip_mash`/`--mash_db`.
- **sourmash**: `SOURMASH_SKETCH` (pools paired reads into one merged signature via
  `--merge`, k=21/scaled=1000 to match this DB) → `SOURMASH_GATHER` (`--threshold-bp
  1000`, see above), toggle with `--skip_sourmash`/`--sourmash_db`.
- **sylph**: `SYLPH_PROFILE` alone - `sylph profile` sketches and profiles reads
  against the database in one step, so no separate sketch stage is needed (unlike
  sourmash). Toggle with `--skip_sylph`/`--sylph_db`.

All 4 `test_full` samples correctly identify their true species with all three tools.
sylph gives the cleanest output of the three - exactly one row per sample, since it
only reports genomes clearing its containment/ANI threshold (98-99.66% ANI on these
samples) - whereas mash returns a distance to every reference (needs sorting to find
the top hit) and sourmash's `gather` can report minor secondary matches to
same-family decoys after explaining away the primary hit (e.g. small hits to
*Klebsiella*/*Salmonella* on the `ECOLI_WGS` sample) - not wrong, but noisier to read.

## Comparison summary

`SPECIES_ID_SUMMARY` (`modules/local/species_id_summary/`, using
`bin/parse_species_id.py`) normalises each tool's own output format into one row -
sample, platform, tool, accession, organism, species taxid, species name, metric
name, metric value - collected into `${outdir}/species_id/species_id_summary.tsv`.
Each tool keeps its own native confidence metric (mash: distance, lower is better;
sourmash: `f_match_orig`, higher is better; sylph: adjusted ANI, higher is better)
rather than trying to normalise them onto one scale. `--species_id_manifest`
(optional; both `test` and `test_full` point it at `manifest.csv` in this directory)
resolves accessions to organism names and species-level taxonomy; without it the
summary reports raw accessions and `NA` taxonomy instead.

Note the three tools don't always agree on *whether* to report a result: against the
purely synthetic `test`-profile fixtures (no real biological content), `sourmash
gather`'s output is declared optional and produces no file at all when nothing clears
its threshold, so its summary row is silently absent for that sample/tool - whereas
`sylph profile` always writes a (possibly empty) output file, which the parser
reports as an explicit "no hit" row. Not a bug, just a real difference in how each
CLI tool signals "nothing found."
