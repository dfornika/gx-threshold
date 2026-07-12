#!/usr/bin/env python3
"""
Classify a sample as pure-culture vs. metagenomic from sourmash gather's own
composition breakdown - gather already decomposes a sample into the set of
reference genomes that best explain it, so "does one genome explain nearly
everything, or does it take many" is read directly off its output, no new
alignment needed.

Reports both a "naive" version (raw per-accession gather rows, as sourmash
reports them) and an "adjusted" version (rows first merged by genomic
similarity cluster - see cluster_reference_genomes.py) side by side, with the
verdict based on the adjusted numbers. The naive numbers are kept, not
discarded, specifically so the adjustment's effect can be spot-checked over
time rather than trusted blindly - e.g. a pure E. coli culture's naive output
includes several Shigella spp. as separate "hits" purely because they are
~98%+ ANI to E. coli (a pre-genomic-era clinical naming artifact, not a real
genomic distinction - the same issue widely seen with Kraken2), which would
otherwise look like spurious extra breadth.

The chosen statistic is f_unique_weighted (gather's incremental, non-
overlapping contribution of each match to the query - the same number
gather's own human-readable summary shows as "p_query"), not f_match_orig
(fraction of the *reference* recovered - a different question, already used
by parse_species_id.py for top-hit confidence reporting).
"""

import argparse
import csv
import gzip
import sys
from collections import defaultdict


def load_clusters(path):
    with open(path) as fh:
        return {row["accession"]: row["cluster_id"] for row in csv.DictReader(fh)}


def load_organisms(path):
    if not path:
        return {}
    with open(path) as fh:
        return {row["accession"]: row.get("organism", row["accession"]) for row in csv.DictReader(fh)}


def load_gather_hits(path):
    """Returns [(accession, organism, f_unique_weighted), ...]."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        rows = list(csv.DictReader(fh))
    hits = []
    for row in rows:
        name_field = row["name"]
        accession, _, organism = name_field.partition(" ")
        hits.append((accession, organism or accession, float(row["f_unique_weighted"])))
    return hits


def summarize(contributions):
    """contributions: {key: fraction}. Returns top_hit_frac, n_hits, effective_n,
    total_explained - effective_n is the inverse Simpson index computed over
    the *explained* fraction only (i.e. normalized so hits sum to 1), so DB
    coverage gaps (real reads with no good reference match at all) don't get
    conflated with "many organisms present"."""
    total_explained = sum(contributions.values())
    n_hits = len(contributions)
    if total_explained <= 0 or n_hits == 0:
        return 0.0, 0, 0.0, total_explained
    top_hit_frac = max(contributions.values())
    sum_sq = sum((c / total_explained) ** 2 for c in contributions.values())
    effective_n = 1 / sum_sq if sum_sq > 0 else 0.0
    return top_hit_frac, n_hits, effective_n, total_explained


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("--clusters", required=True, help="output of cluster_reference_genomes.py")
    ap.add_argument("--manifest", default=None, help="species_db manifest.csv, for organism names (optional)")
    ap.add_argument("--effective-n-threshold", type=float, default=1.5, help="adjusted effective_n at/below this -> pure culture")
    ap.add_argument("--min-explained-frac", type=float, default=0.05, help="below this fraction of k-mers explained by any reference, the verdict is unreliable regardless of shape - report inconclusive instead")
    ap.add_argument("gather_csv", help="sourmash gather CSV(.gz) for this sample")
    args = ap.parse_args()

    clusters = load_clusters(args.clusters)
    organisms = load_organisms(args.manifest)
    hits = load_gather_hits(args.gather_csv)

    if not hits:
        print(f"{args.sample}\t{args.platform}\tinconclusive\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tno gather hits")
        return

    naive_contributions = defaultdict(float)
    adjusted_contributions = defaultdict(float)
    adjusted_members = defaultdict(set)
    for accession, organism, contribution in hits:
        naive_contributions[organism] += contribution
        cluster_id = clusters.get(accession, accession)
        adjusted_contributions[cluster_id] += contribution
        adjusted_members[cluster_id].add(organisms.get(accession, organism))

    naive_top, naive_n, naive_eff, naive_explained = summarize(naive_contributions)
    adj_top, adj_n, adj_eff, adj_explained = summarize(adjusted_contributions)

    # A low explained fraction means the reference DB doesn't cover this
    # sample well - found empirically (see docs/testing.md) on a real gut
    # metagenome whose actual organisms weren't in our panel, which produced
    # only one weak hit and so looked exactly like a confident pure culture.
    # Breadth is only a meaningful signal when there's enough evidence to
    # measure it from - otherwise this is a DB-coverage gap, not a purity
    # verdict.
    if adj_explained < args.min_explained_frac:
        verdict = "inconclusive"
    else:
        verdict = "pure_culture" if adj_eff <= args.effective_n_threshold else "metagenomic"

    collapsed = [
        "+".join(sorted(members))
        for members in adjusted_members.values()
        if len(members) > 1
    ]
    notes = collapsed[:]
    if verdict == "inconclusive":
        notes.insert(0, f"only {adj_explained:.2%} of k-mers explained by any reference (min {args.min_explained_frac:.0%})")
    note = ";".join(notes) if notes else "-"

    print(
        f"{args.sample}\t{args.platform}\t{verdict}\t"
        f"{naive_n}\t{naive_top:.4g}\t{naive_eff:.4g}\t"
        f"{adj_n}\t{adj_top:.4g}\t{adj_eff:.4g}\t{adj_explained:.4g}\t{note}"
    )


if __name__ == "__main__":
    sys.exit(main())
