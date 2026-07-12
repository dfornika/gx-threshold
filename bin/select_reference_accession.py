#!/usr/bin/env python3
"""
Pick one reference genome accession per sample from the (up to 3) species-ID
tool calls produced by parse_species_id.py, so downstream alignment has a
single genome to work with.

Rule: majority vote on species_taxid among whichever tools produced a hit for
this sample (2-of-2 or 2-of-3 agreement counts as consensus; a lone hit is
trivially its own consensus). If there's no majority (all tools disagree, or
none produced a hit), fall back to a fixed priority order - sylph, then mash,
then sourmash - based on which gave the cleanest results in our own evaluation
(see tests/data/species_db/README.md). This does not attempt anything fancier
(e.g. weighting by each tool's own confidence metric), since those metrics
aren't on comparable scales across tools.
"""

import argparse
import csv
import sys
from collections import Counter

FALLBACK_PRIORITY = ["sylph", "mash", "sourmash"]

FIELDS = ["sample", "platform", "tool", "accession", "organism", "species_taxid", "species_name", "metric", "value"]


def read_rows(paths):
    rows = []
    for path in paths:
        with open(path) as fh:
            for line in fh:
                rows.append(dict(zip(FIELDS, line.rstrip("\n").split("\t"))))
    return rows


def select(rows):
    valid = [r for r in rows if r["accession"] != "NA" and r["species_taxid"] != "NA"]
    if not valid:
        return None

    votes = Counter(r["species_taxid"] for r in valid)
    top_taxid, top_count = votes.most_common(1)[0]

    if top_count > 1:
        winners = [r for r in valid if r["species_taxid"] == top_taxid]
        for tool in FALLBACK_PRIORITY:
            for r in winners:
                if r["tool"] == tool:
                    return r, f"majority({top_count}/{len(valid)})"
        return winners[0], f"majority({top_count}/{len(valid)})"

    for tool in FALLBACK_PRIORITY:
        for r in valid:
            if r["tool"] == tool:
                return r, f"fallback:{tool}"
    return valid[0], "fallback:only-candidate"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("rows", nargs="+", help="one or more *.species_id.tsv files for this sample (from parse_species_id.py)")
    args = ap.parse_args()

    rows = read_rows(args.rows)
    result = select(rows)

    if result is None:
        print(f"{args.sample}\t{args.platform}\tNA\tNA\tNA\tno_consensus")
        return

    row, method = result
    print(f"{args.sample}\t{args.platform}\t{row['accession']}\t{row['species_taxid']}\t{row['species_name']}\t{method}")


if __name__ == "__main__":
    sys.exit(main())
