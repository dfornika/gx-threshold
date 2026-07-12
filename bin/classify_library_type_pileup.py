#!/usr/bin/env python3
"""
Classify a sample's library prep as amplicon vs. shotgun from aligned
read/fragment start-end positions - "does this read start (and end) at the
same place as many other reads?"

Amplicon libraries repeatedly re-sequence the same PCR product, so most
reads pile up at a small number of essentially fixed fragment coordinates
(the primer sites). Shotgun libraries fragment DNA close to randomly, so two
independent fragments sharing the same start AND end coordinate is rare -
this is exactly the same positional signature real duplicate-marking tools
(Picard/samtools markdup) use to flag PCR/optical duplicates. The difference
here is that for amplicon libraries, this "duplication rate" isn't a
side-effect of over-amplification to police - it's most of the data.

Needs reads aligned to a reference genome (see modules/local/align_reads/),
unlike classify_library_type_cluster.py, which clusters reads against each
other with no reference at all.

Nearby coordinates are grouped by single-linkage chaining (sort, then merge
consecutive positions within a per-platform tolerance) rather than requiring
an exact match or rounding to a fixed grid - a fixed grid can arbitrarily
split one true pileup in two when its spread straddles a grid boundary,
which is exactly what happened during development with real Nanopore
amplicon data (primer-trim slop of just a few bases split what should have
been one large pileup across several grid cells - see docs/testing.md).

For paired-end input, one signature is used per *fragment*, not per read:
POS and PNEXT/TLEN give the true fragment span regardless of which mate is
being read, matching how real duplicate-marking tools group pairs. Both the
start and end coordinate are required to match (within tolerance) for two
fragments to share a signature, not just the start - this matters more than
it might seem: requiring only a shared start let bacterial genome repeat
regions (rRNA operons etc. - the same confound found while validating
classify_library_type_aligned.awk) masquerade as amplicon pileups, since
reads from different genomic copies of a real repeat can genuinely start
near each other by chance. The end coordinate is noisier for single-end long
reads (variable soft-clipping, no mate to cross-check against), so its own
tolerance can be set looser than the start coordinate's.
"""

import argparse
import sys
from collections import Counter, defaultdict


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


def parse_fragments(sam_path, single_end):
    """Yield (rname, start, end) for one alignment record per fragment."""
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            cigar = fields[5]
            pnext = int(fields[7])
            tlen = int(fields[8])

            if hasflag(flag, 4) or hasflag(flag, 256) or hasflag(flag, 2048):
                continue  # unmapped / secondary / supplementary
            if rname == "*" or cigar == "*":
                continue

            if single_end:
                length = reflen(cigar)
                if length < 1:
                    continue
                yield rname, pos, pos + length - 1
            else:
                if not (hasflag(flag, 1) and hasflag(flag, 2)):
                    continue  # need a properly-paired read
                if hasflag(flag, 128):
                    continue  # one signature per fragment: first-in-pair only
                if tlen == 0:
                    continue
                start = min(pos, pnext)
                yield rname, start, start + abs(tlen) - 1


def chain_group(items, tolerance):
    """Single-linkage-cluster (value, id) pairs by sorted value: consecutive
    items merge into the same group while the gap to the previous one stays
    within `tolerance`. Returns {id: group_index}."""
    group_of = {}
    group_id = -1
    prev_value = None
    for value, item_id in sorted(items, key=lambda vi: vi[0]):
        if prev_value is None or value - prev_value > tolerance:
            group_id += 1
        group_of[item_id] = group_id
        prev_value = value
    return group_of


def fragment_signatures(sam_path, single_end, start_tolerance, end_tolerance):
    fragments = list(parse_fragments(sam_path, single_end))

    by_rname = defaultdict(list)
    for idx, (rname, _start, _end) in enumerate(fragments):
        by_rname[rname].append(idx)

    signature_of = {}
    for rname, idxs in by_rname.items():
        start_group = chain_group([(fragments[i][1], i) for i in idxs], start_tolerance)

        by_start_group = defaultdict(list)
        for i in idxs:
            by_start_group[start_group[i]].append(i)
        for sg, sub_idxs in by_start_group.items():
            end_group = chain_group([(fragments[i][2], i) for i in sub_idxs], end_tolerance)
            for i in sub_idxs:
                signature_of[i] = (rname, sg, end_group[i])

    return Counter(signature_of.values()), len(fragments)


def summarize(signatures, total_reads):
    sizes = list(signatures.values())
    n_reads_piled = sum(s for s in sizes if s > 1)
    largest = max(sizes) if sizes else 0

    largest_pileup_frac = largest / total_reads
    piled_frac = n_reads_piled / total_reads

    sum_sq_frac = sum((s / total_reads) ** 2 for s in sizes)
    effective_signatures = (1 / sum_sq_frac) if sum_sq_frac > 0 else 0.0
    effective_signatures_frac = effective_signatures / total_reads

    return {
        "largest_pileup_frac": largest_pileup_frac,
        "piled_frac": piled_frac,
        "effective_signatures_frac": effective_signatures_frac,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("--single-end", action="store_true")
    ap.add_argument("--start-tolerance", type=int, required=True, help="bp window for chaining nearby fragment start coordinates together")
    ap.add_argument("--end-tolerance", type=int, required=True, help="bp window for chaining nearby fragment end coordinates together")
    ap.add_argument("--effective-signatures-frac-threshold", type=float, default=0.5, help="Below this, call amplicon")
    ap.add_argument("sam", help="Plain-text SAM from modules/local/align_reads/")
    args = ap.parse_args()

    signatures, total_reads = fragment_signatures(args.sam, args.single_end, args.start_tolerance, args.end_tolerance)

    if total_reads == 0:
        print(f"{args.sample}\t{args.platform}\tinconclusive\t0\tNA\tNA\tNA\tno usable alignments")
        return

    stats = summarize(signatures, total_reads)
    verdict = "amplicon" if stats["effective_signatures_frac"] < args.effective_signatures_frac_threshold else "shotgun"

    print(
        f"{args.sample}\t{args.platform}\t{verdict}\t{total_reads}\t"
        f"{stats['largest_pileup_frac']:.4g}\t{stats['piled_frac']:.4g}\t"
        f"{stats['effective_signatures_frac']:.4g}\t-"
    )


if __name__ == "__main__":
    sys.exit(main())
