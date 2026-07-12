#!/usr/bin/env python3
"""
Classify a sample's library prep as amplicon vs. shotgun from its fastp/fastplong
QC JSON. Illumina paired-end only, validated against tests/data/real/ (two shotgun
WGS samples, one amplicon): amplicon libraries show a sharp, narrow peak in the
insert-size histogram (reads cluster around the fixed amplicon length) and
elevated PCR duplication, vs. shotgun's broad insert-size spread and near-zero
duplication.

Nanopore/long-read is not classified here: fastplong's JSON reports neither a
duplication rate nor an insert-size histogram (fastp does; fastplong doesn't -
confirmed by inspecting real output, not a doc assumption), and raw read-length
distribution on our one real ONT sample didn't show a comparably clean signal
either. Real long-read amplicon detection likely needs an alignment-based
approach (read start/end clustering against a reference) - see
tests/data/species_db/README.md and docs/testing.md.
"""

import argparse
import json
import sys

CONCENTRATION_THRESHOLD = 0.20  # fraction of read pairs within +/- PEAK_WINDOW bp of the insert-size peak
PEAK_WINDOW = 10  # bp


def classify(report):
    insert_size = report.get("insert_size")
    duplication = report.get("duplication")
    if insert_size is None or duplication is None:
        return None, "fastplong reports no duplication rate or insert-size histogram (long-read not classified)"

    histogram = insert_size["histogram"]
    peak = insert_size["peak"]
    known_pairs = sum(histogram)
    if known_pairs == 0:
        return None, "no read pairs with a determined insert size"

    lo, hi = max(0, peak - PEAK_WINDOW), min(len(histogram), peak + PEAK_WINDOW + 1)
    concentration = sum(histogram[lo:hi]) / known_pairs
    dup_rate = duplication["rate"]
    verdict = "amplicon" if concentration > CONCENTRATION_THRESHOLD else "shotgun"
    return (verdict, concentration, dup_rate, peak), None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("qc_json", help="fastp or fastplong JSON report for this sample")
    args = ap.parse_args()

    with open(args.qc_json) as fh:
        report = json.load(fh)

    result, reason = classify(report)
    if result is None:
        print(f"{args.sample}\t{args.platform}\tnot_classified\tNA\tNA\tNA\t{reason}")
        return

    verdict, concentration, dup_rate, peak = result
    print(f"{args.sample}\t{args.platform}\t{verdict}\t{concentration:.4g}\t{dup_rate:.4g}\t{peak}\t-")


if __name__ == "__main__":
    sys.exit(main())
