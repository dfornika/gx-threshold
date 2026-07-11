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
over what each test contains (e.g. a known host fraction). The longer-term intent is a
**mix**: keep synthetic fixtures for deterministic unit tests, and add a small number
of **downsampled real** samples for a more realistic end-to-end check. As the dataset
grows, test inputs are expected to move to versioned URLs (hosted in our own
test-datasets repository) referenced from the test configs, mirroring the nf-core
convention.

## Testing strategy

Tests are layered:

1. **Module tests** (`modules/**/tests/*.nf.test`) — exercise a single process in
   isolation with the tiny fixtures. These favour **explicit assertions** (exact read
   counts, invariants) for correctness, plus snapshots to catch unintended drift.
2. **End-to-end / smoke test** — run the whole pipeline on the `test` profile to
   confirm everything wires together and produces the expected outputs, fully offline.
3. **Realistic test** (`test_full`, planned) — a larger, more representative dataset
   for release validation.

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
