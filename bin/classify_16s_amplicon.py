#!/usr/bin/env python3
"""
Classify whether a sample's reads look like 16S rRNA amplicon sequencing, by
aligning (a subset of) them to a dedicated 16S reference database and
checking what fraction pass an identity+coverage bar against *any* reference
in it - not which specific organism, just "is this 16S at all".

This exists because sourmash/mash (k-mer/minhash containment) turned out to
be unreliable for this - 16S has long regions conserved across nearly all
bacteria (that's what makes it a useful universal marker gene), so short
k-mers spuriously match many unrelated reference sequences regardless of k
(tried k=21 and k=31 - see data/dev_16s_db/README.md). Full alignment doesn't
have that weakness: identity and coverage are measured over the read's whole
length, so a read genuinely from one organism's 16S gene aligns with high
identity across nearly all of it, while conserved-region-only noise gets
filtered out by the coverage requirement - the same identity+coverage logic
already used in classify_library_type_pileup.py, just applied to "does this
read match the 16S database at all" instead of "does this read share a
position with others".

Only a subset of reads needs checking (this is a coarse yes/no question, not
a composition breakdown), so the caller aligns a small subsample for speed -
see modules/local/subsample_reads_head/.
"""

import argparse
import sys


def hasflag(flag, bit):
    return (flag // bit) % 2 == 1


def reflen(cigar):
    """Sum of CIGAR operations that consume the reference (M/D/N/=/X)."""
    length = 0
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            if ch in "MDN=X":
                length += int(num) if num else 0
            num = ""
    return length


def cigar_matches(cigar):
    """Sum of CIGAR M/=/X operations - an approximate alignment-block length
    (matches + mismatches), for an identity estimate when the SAM has no NM
    tag. Good enough for a coarse yes/no check, not claimed as exact ANI."""
    length = 0
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            if ch in "M=X":
                length += int(num) if num else 0
            num = ""
    return length


def parse_nm_tag(fields):
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            return int(field[5:])
    return None


def evaluate_reads(sam_path, min_identity, min_coverage):
    # `total` counts every primary record (mapped or not) - a read that
    # aligns to nothing in the 16S database is real evidence *against* 16S,
    # not a missing observation, so it belongs in the denominator. Only
    # secondary/supplementary records are excluded (not independent reads).
    total = 0
    passed = 0

    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            rname = fields[2]
            cigar = fields[5]
            seq = fields[9]

            if hasflag(flag, 256) or hasflag(flag, 2048):
                continue  # secondary / supplementary - not an independent read
            total += 1

            if hasflag(flag, 4) or rname == "*" or cigar == "*" or seq == "*":
                continue  # unmapped - counted above, but nothing to score

            ref_len = reflen(cigar)
            aln_block = cigar_matches(cigar)
            if ref_len < 1 or aln_block < 1:
                continue

            nm = parse_nm_tag(fields)
            identity = 1 - (nm / aln_block) if nm is not None else 1.0
            coverage = ref_len / len(seq) if seq else 0.0
            # A read's own aligned block can exceed its raw length once
            # insertions are counted (ref_len only counts M/D/N/=/X, not I) -
            # cap at 1.0 rather than let that read a coverage bar it didn't
            # really clear.
            coverage = min(coverage, 1.0)

            if identity >= min_identity and coverage >= min_coverage:
                passed += 1

    return total, passed


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("--min-identity", type=float, required=True)
    ap.add_argument("--min-coverage", type=float, required=True)
    ap.add_argument("--min-passed-frac", type=float, default=0.5, help="at/above this fraction of reads passing the identity+coverage bar -> 16S_amplicon")
    ap.add_argument("sam", help="minimap2 SAM alignment against the 16S reference database")
    args = ap.parse_args()

    total, passed = evaluate_reads(args.sam, args.min_identity, args.min_coverage)

    if total == 0:
        print(f"{args.sample}\t{args.platform}\tinconclusive\t0\t0\tNA\tno usable alignments")
        return

    passed_frac = passed / total
    verdict = "16S_amplicon" if passed_frac >= args.min_passed_frac else "other"

    print(f"{args.sample}\t{args.platform}\t{verdict}\t{total}\t{passed}\t{passed_frac:.4g}\t-")


if __name__ == "__main__":
    sys.exit(main())
