# Real (downsampled) test data

Small, real-read fixtures for the `test_full` profile, complementing the fully
synthetic fixtures in [`tests/data/`](../) used by the default `test` profile.
Unlike those, these reads are genuine sequencer output subsampled down to a
few hundred/thousand reads per sample - real base composition and quality
profiles, at synthetic-fixture size - with a small number of real-human-derived
reads spiked in so host-read removal is actually exercised, too.

```
real/
  references/
    human_mt.fasta       real human mitochondrial genome (rCRS, NC_012920.1)
  spike_human.py         regenerates reads/ from reads_src/ + human_mt.fasta
  reads/                 reads_src/ + spiked-in human reads - this is the samplesheet input, committed
    illumina/  ...
    nanopore/  ...
  samplesheet.csv
  reads_src/              *not committed* (gitignored) - pristine pre-spike SRA
                          reads; see "Regenerating" below to recreate it
```

Only `reads/` (the final, spiked fixtures - ~1.7MB total) is committed.
`reads_src/` is a local working layer, not present in a fresh clone.

## Provenance: the underlying microbial/viral reads

Each FASTQ pair/file was downloaded from the ENA Portal API (mirrors SRA) and
then randomly subsampled with [`rasusa`](https://github.com/mbhall88/rasusa)
(`reads -n <N> -s 42`, fixed seed).

| Sample | SRA run | Organism | Platform | Library | Subsampled to |
|---|---|---|---|---|---|
| `KPNEUMONIAE_WGS` | [SRR14584974](https://www.ncbi.nlm.nih.gov/sra/SRR14584974) | *Klebsiella pneumoniae* | Illumina MiSeq, paired-end | WGS | 1200 read pairs |
| `ECOLI_WGS` | [SRR12695070](https://www.ncbi.nlm.nih.gov/sra/SRR12695070) | *Escherichia coli* | Illumina MiSeq, paired-end | WGS | 1200 read pairs |
| `SARS2_AMPLICON_ILLUMINA` | [SRR18111140](https://www.ncbi.nlm.nih.gov/sra/SRR18111140) | SARS-CoV-2 | Illumina MiSeq, paired-end | PCR tiling amplicon | 800 read pairs |
| `SARS2_AMPLICON_ONT` | [SRR22939730](https://www.ncbi.nlm.nih.gov/sra/SRR22939730) | SARS-CoV-2 | Oxford Nanopore MinION, single-end | PCR tiling amplicon | 800 reads |

Downloaded via the ENA FASTQ convenience path (`download_fastq`) in
[`ncbi-client-py`](https://github.com/dfornika/ncbi-client-py); subsampled
with `rasusa reads`.

## The human host spike-in (`reads/`, via `spike_human.py`)

`references/human_mt.fasta` is the real human mitochondrial genome (the rCRS,
[NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1), 16,569 bp),
fetched via `efetch` in `ncbi-client-py`. It stands in for a full masked human
reference: small enough to commit in-repo, but unambiguously real human DNA
(unlike `tests/data/references/host_mini.fasta`, which is synthetic).

[`spike_human.py`](spike_human.py) simulates a small number of Illumina- and
Nanopore-like reads from this reference (same seeded substitution-mutation
approach as `tests/data/generate.py` - stdlib `random` only, fixed seed) and
appends them to each sample from `reads_src/`, writing the result to `reads/`
(what `samplesheet.csv` actually points at):

| Sample | source reads | + human reads | host fraction |
|---|---|---|---|
| `KPNEUMONIAE_WGS` | 1200 pairs | 25 pairs | ~2% |
| `ECOLI_WGS` | 1200 pairs | 25 pairs | ~2% |
| `SARS2_AMPLICON_ILLUMINA` | 800 pairs | 15 pairs | ~1.9% |
| `SARS2_AMPLICON_ONT` | 800 reads | 15 reads | ~1.9% |

The script always regenerates `reads/` from `reads_src/` + the reference from
scratch, so re-running it is idempotent (safe to re-run; it never compounds
spike-in on top of itself) - **provided `reads_src/` is present locally** (see
"Regenerating" below, since it isn't committed).

## What this dataset exercises

Real read-quality profiles (real error/adapter/quality-score distributions)
across both platforms and two library types (WGS vs. amplicon), **and** real
host-read removal against a genuinely human reference - `DEHOST` should report
a small (~2%) but nonzero host fraction for every sample here, with the
microbial/viral reads retained. This is a different check than the synthetic
host+microbe mix in `tests/data/` (which uses exact synthetic sequences purely
for deterministic, exact-count assertions in module tests) - this dataset's
value is realism, not exact reproducible ground truth of the underlying reads
(only the spike-in step is deterministic; the microbial/viral reads themselves
are a frozen real-world snapshot, and not committed in raw form - see below).

## Regenerating / extending

`reads_src/` isn't committed (kept small in git; see `.gitignore`), so it must
be recreated locally before `spike_human.py` can run. To rebuild it from
scratch, or to refresh/add an accession:

1. Pick a small SRA run (check `fastq_ftp` is populated via the ENA Portal
   API filereport endpoint - very recent submissions may not be synced yet).
2. Download via `ncbi_client.sra.download_fastq(client, accession, dest)`.
3. Subsample directly to the target count from the table above, e.g.:
   `rasusa reads -n 1200 -s 42 -o reads_src/illumina/KPNEUMONIAE_WGS_R1.fastq.gz -o reads_src/illumina/KPNEUMONIAE_WGS_R2.fastq.gz in_R1.fastq.gz in_R2.fastq.gz`
4. Repeat for each sample in the table (matching path/filename conventions
   under `reads_src/illumina/` or `reads_src/nanopore/`).
5. Run `python3 tests/data/real/spike_human.py` to produce the final,
   committed `reads/`.

To add a new sample: do the above for the new accession, then add an entry to
`ILLUMINA_SAMPLES`/`NANOPORE_SAMPLES` in `spike_human.py`, update
`samplesheet.csv`, and update the tables in this file.
