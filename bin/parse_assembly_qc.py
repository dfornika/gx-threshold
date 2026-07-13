#!/usr/bin/env python3
"""
Normalise the assembly-QC outputs for one sample into a single summary row, so
QUAST (contiguity) and CheckM2 (completeness/contamination) can be read side by
side per sample - the same one-row-per-sample-per-stage convention every other
stage in this pipeline uses (cf. bin/parse_species_id.py).

QUAST's `report.tsv` has one metric per row (metric name in column 1, value in
column 2). CheckM2's `quality_report.tsv` is a normal header + one data row per
input genome (we give it exactly one assembly, so one row). CheckM2 is
optional - it only runs when a --checkm2_db was provided - so its columns come
back as NA when absent.
"""

import argparse
import csv
import sys

NA = "NA"


def parse_quast(path):
    """QUAST report.tsv -> {metric: value}. Exact-key lookups below pick the
    plain '# contigs'/'Total length' rows, not the '(>= N bp)' threshold variants."""
    metrics = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                metrics[parts[0]] = parts[1]
    return {
        "n_contigs": metrics.get("# contigs", NA),
        "total_length": metrics.get("Total length", NA),
        "largest_contig": metrics.get("Largest contig", NA),
        "n50": metrics.get("N50", NA),
    }


def parse_checkm2(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return {"completeness": NA, "contamination": NA}
    row = rows[0]
    return {
        "completeness": row.get("Completeness", NA),
        "contamination": row.get("Contamination", NA),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("--assembler", required=True, help="shovill (Illumina) or dragonflye (Nanopore)")
    ap.add_argument("--quast", required=True, help="QUAST report.tsv for this sample")
    ap.add_argument("--checkm2", default=None, help="CheckM2 quality_report.tsv (omit if CheckM2 did not run)")
    args = ap.parse_args()

    q = parse_quast(args.quast)
    c = parse_checkm2(args.checkm2) if args.checkm2 else {"completeness": NA, "contamination": NA}

    print(
        f"{args.sample}\t{args.platform}\t{args.assembler}\t"
        f"{q['n_contigs']}\t{q['total_length']}\t{q['largest_contig']}\t{q['n50']}\t"
        f"{c['completeness']}\t{c['contamination']}"
    )


if __name__ == "__main__":
    sys.exit(main())
