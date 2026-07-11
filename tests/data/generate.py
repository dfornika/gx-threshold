#!/usr/bin/env python3
"""
Generate tiny, deterministic synthetic test fixtures for the gx-threshold test
profile and nf-test unit tests.

Everything is seeded (stdlib `random` only), so re-running reproduces byte-identical
FASTA/FASTQ. Regenerate with:  python3 tests/data/generate.py

Design:
  - Two unrelated pseudo-random references: a "host" and a "microbe". Because they
    share no meaningful k-mers, reads from one do not map to the other.
  - Illumina PE and Nanopore SE read sets are a MIX of host-derived and
    microbe-derived reads. After dehosting against host_mini.fasta, the host-derived
    reads are removed and the microbe-derived reads are retained - so the fixtures
    exercise the real removal path (not just "keep everything").
  - A small mutation rate keeps reads realistic while still mapping reliably.

The exact host/microbe read counts are printed as a manifest so tests can assert
against them.
"""

import gzip
import os
import random
import textwrap

SEED = 20260711
RNG = random.Random(SEED)

HERE = os.path.dirname(os.path.abspath(__file__))
REF_DIR = os.path.join(HERE, "references")
ILLUMINA_DIR = os.path.join(HERE, "reads", "illumina")
NANOPORE_DIR = os.path.join(HERE, "reads", "nanopore")

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def random_sequence(length):
    return "".join(RNG.choice(BASES) for _ in range(length))


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def mutate(seq, rate=0.01):
    """Introduce substitutions at the given per-base rate."""
    out = list(seq)
    for i, b in enumerate(out):
        if RNG.random() < rate:
            out[i] = RNG.choice([x for x in BASES if x != b])
    return "".join(out)


def write_fasta(path, name, seq, width=70):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        fh.write("\n".join(textwrap.wrap(seq, width)) + "\n")


def write_fastq_gz(path, records):
    """records: list of (read_id, seq, qual). mtime=0 keeps the gzip byte-identical."""
    payload = "".join(f"@{rid}\n{seq}\n+\n{qual}\n" for rid, seq, qual in records)
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
            gz.write(payload.encode())


def make_illumina_pair(ref, read_len=150, insert=350):
    """Return (r1_seq, r2_seq) as a proper FR pair drawn from `ref`."""
    start = RNG.randint(0, len(ref) - insert)
    frag = ref[start:start + insert]
    r1 = mutate(frag[:read_len])
    r2 = mutate(revcomp(frag[-read_len:]))
    return r1, r2


def make_ont_read(ref, min_len=800, max_len=3000):
    length = RNG.randint(min_len, max_len)
    start = RNG.randint(0, len(ref) - length)
    seq = mutate(ref[start:start + length], rate=0.03)  # ONT: higher error
    return seq


def main():
    for d in (REF_DIR, ILLUMINA_DIR, NANOPORE_DIR):
        os.makedirs(d, exist_ok=True)

    # --- references (independent -> no cross-mapping) ---
    host = random_sequence(40000)
    microbe = random_sequence(35000)
    write_fasta(os.path.join(REF_DIR, "host_mini.fasta"), "host_mini_chr", host)
    write_fasta(os.path.join(REF_DIR, "microbe_mini.fasta"), "microbe_mini_chr", microbe)

    # --- Illumina PE: microbe (keep) + host (remove) ---
    n_microbe_pe, n_host_pe = 150, 40
    r1_records, r2_records = [], []
    for i in range(n_microbe_pe):
        r1, r2 = make_illumina_pair(microbe)
        q = "I" * len(r1)
        r1_records.append((f"microbe_pe_{i}", r1, q))
        r2_records.append((f"microbe_pe_{i}", r2, "I" * len(r2)))
    for i in range(n_host_pe):
        r1, r2 = make_illumina_pair(host)
        r1_records.append((f"host_pe_{i}", r1, "I" * len(r1)))
        r2_records.append((f"host_pe_{i}", r2, "I" * len(r2)))
    write_fastq_gz(os.path.join(ILLUMINA_DIR, "SAMPLE_PE_R1.fastq.gz"), r1_records)
    write_fastq_gz(os.path.join(ILLUMINA_DIR, "SAMPLE_PE_R2.fastq.gz"), r2_records)

    # --- Nanopore SE: microbe (keep) + host (remove) ---
    n_microbe_ont, n_host_ont = 60, 15
    ont_records = []
    for i in range(n_microbe_ont):
        s = make_ont_read(microbe)
        ont_records.append((f"microbe_ont_{i}", s, "5" * len(s)))
    for i in range(n_host_ont):
        s = make_ont_read(host)
        ont_records.append((f"host_ont_{i}", s, "5" * len(s)))
    write_fastq_gz(os.path.join(NANOPORE_DIR, "SAMPLE_ONT.fastq.gz"), ont_records)

    # --- manifest (ground truth for assertions) ---
    print("Generated fixtures (seed={}):".format(SEED))
    print(f"  references/host_mini.fasta      {len(host)} bp")
    print(f"  references/microbe_mini.fasta   {len(microbe)} bp")
    print(f"  Illumina PE pairs: {n_microbe_pe + n_host_pe} "
          f"(microbe/keep={n_microbe_pe}, host/remove={n_host_pe})")
    print(f"  Nanopore SE reads: {n_microbe_ont + n_host_ont} "
          f"(microbe/keep={n_microbe_ont}, host/remove={n_host_ont})")


if __name__ == "__main__":
    main()
