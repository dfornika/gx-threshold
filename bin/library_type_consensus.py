#!/usr/bin/env python3
"""
Fuse the (up to 4) amplicon-vs-shotgun library-type calls - LIBRARY_TYPE
(fastp/fastplong QC JSON), LIBRARY_TYPE_CLUSTER (read clustering),
LIBRARY_TYPE_ALIGNED (depth dispersion), LIBRARY_TYPE_PILEUP (read-position
pileup) - into one consensus verdict per sample, so meta.library_type is
always populated instead of only ever reflecting one method's own call.

Rule: majority vote among whichever methods produced a real verdict
(amplicon/shotgun) for this sample - a method that reported
not_classified/inconclusive doesn't count as a vote either way, same idea as
select_reference_accession.py's own "no hit" filter. A tie is broken by a
fixed, PROVISIONAL fallback priority - pileup, then aligned, then cluster,
then fastp_json - based on docs/testing.md's current empirical read: pileup
has no known failure case yet; aligned and cluster share the same
bacterial-genome-repeat-structure confound; fastp_json is Illumina-only and
validated on the fewest samples. This priority order is expected to be
revisited once more real validation data (data/dev_samples/,
data/dev_metagenomic_samples/) shows whether any one or two methods are
consistently correct on their own.
"""

import argparse
import sys
from collections import Counter

FALLBACK_PRIORITY = ["pileup", "aligned", "cluster", "fastp_json"]
VALID_VERDICTS = {"amplicon", "shotgun"}


def read_verdict(path):
    with open(path) as fh:
        first_line = fh.readline().rstrip("\n")
    columns = first_line.split("\t")
    return columns[2] if len(columns) > 2 else "NA"


def vote(calls):
    """calls: list of (method, verdict). Returns (verdict, method_note) or None if no valid votes."""
    valid = [(method, verdict) for method, verdict in calls if verdict in VALID_VERDICTS]
    if not valid:
        return None

    counts = Counter(verdict for _method, verdict in valid)
    top_verdict, top_count = counts.most_common(1)[0]
    total = len(valid)

    if top_count > total - top_count:
        return top_verdict, f"majority({top_count}/{total})"

    # Tie: break by fixed priority order among the methods that voted for
    # either verdict. (total == 1 can never land here - a single vote is
    # always its own majority via the check above.)
    by_method = dict(valid)
    for method in FALLBACK_PRIORITY:
        if method in by_method:
            return by_method[method], f"tie_break:{method}"
    return valid[0][1], f"tie_break:{valid[0][0]}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument(
        "--call",
        action="append",
        default=[],
        metavar="METHOD:PATH",
        help="one per method that ran for this sample, e.g. --call cluster:SAMPLE.library_type_cluster.tsv "
        "(repeatable). METHOD should be one of: fastp_json, cluster, aligned, pileup.",
    )
    args = ap.parse_args()

    calls = []
    for call in args.call:
        method, _, path = call.partition(":")
        calls.append((method, read_verdict(path)))

    n_amplicon = sum(1 for _method, verdict in calls if verdict == "amplicon")
    n_shotgun = sum(1 for _method, verdict in calls if verdict == "shotgun")
    n_total = len(calls)

    result = vote(calls)
    if result is None:
        print(f"{args.sample}\t{args.platform}\tno_data\t{n_amplicon}\t{n_shotgun}\t{n_total}\tno_data\tno method produced a classified verdict")
        return

    verdict, method_note = result
    print(f"{args.sample}\t{args.platform}\t{verdict}\t{n_amplicon}\t{n_shotgun}\t{n_total}\t{method_note}\t-")


if __name__ == "__main__":
    sys.exit(main())
