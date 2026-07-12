#!/usr/bin/env python3
"""
Cluster a mash reference database's own genomes by genomic similarity (ANI),
so downstream composition-breadth analysis (classify_species_composition.py)
can tell "many distinct organisms" apart from "one organism whose nomenclature
happens to span several named species" - e.g. Escherichia coli and Shigella
spp., which are ~98%+ ANI to each other (a pre-genomic-era clinical naming
artifact, not a real genomic distinction) and so reliably show up together as
spurious extra "hits" for what is actually a pure E. coli culture.

Approach: given `mash dist`'s all-vs-all output on the database's own sketch
file (the same file already used for species-ID - see
modules/local/reference_genome_distances/, a separate Nextflow process since
`mash` and this script's own union-find logic need different containers),
convert each pairwise distance to an ANI estimate (ANI ~= (1 - distance) * 100
- the standard mash-distance-to-ANI approximation), and union-find genomes
together whenever their estimated ANI clears `--ani-threshold` (default 95%,
the conventional operational species boundary - the same one GTDB uses, which
is why GTDB's own taxonomy already merges Shigella into Escherichia coli).
This runs once per database, not per sample - the output is a lookup table
reused across every sample's gather results.
"""

import argparse
import csv
import sys
from pathlib import Path


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def accession_from_path(path):
    return Path(path).name.removesuffix(".fasta").removesuffix(".fa")


def cluster_genomes(distances_path, ani_threshold):
    uf = UnionFind()
    accessions = set()
    with open(distances_path) as fh:
        for line in fh:
            ref, query, distance, _pvalue, _shared = line.rstrip("\n").split("\t")
            ref_acc = accession_from_path(ref)
            query_acc = accession_from_path(query)
            accessions.add(ref_acc)
            accessions.add(query_acc)
            if ref_acc == query_acc:
                continue
            ani = (1 - float(distance)) * 100
            if ani >= ani_threshold:
                uf.union(ref_acc, query_acc)
    return uf, accessions


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ani-threshold", type=float, default=95.0)
    ap.add_argument("--manifest", default=None, help="species_db manifest.csv, for organism names in the output (optional)")
    ap.add_argument("distances", help="mash dist all-vs-all output on the reference database's own sketch file")
    ap.add_argument("output", help="output CSV: accession,cluster_id,organism")
    args = ap.parse_args()

    manifest = {}
    if args.manifest:
        with open(args.manifest) as fh:
            manifest = {row["accession"]: row.get("organism", "") for row in csv.DictReader(fh)}

    uf, accessions = cluster_genomes(args.distances, args.ani_threshold)

    # Use the lowest accession in each cluster as a stable, deterministic
    # cluster id (rather than an arbitrary union-find root).
    cluster_members = {}
    for acc in accessions:
        cluster_members.setdefault(uf.find(acc), []).append(acc)
    cluster_id_of = {}
    for members in cluster_members.values():
        cluster_id = min(members)
        for acc in members:
            cluster_id_of[acc] = cluster_id

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["accession", "cluster_id", "organism"])
        for acc in sorted(accessions):
            writer.writerow([acc, cluster_id_of[acc], manifest.get(acc, "")])

    n_clusters = len(cluster_members)
    n_merged = sum(1 for members in cluster_members.values() if len(members) > 1)
    print(
        f"{len(accessions)} genomes -> {n_clusters} clusters at {args.ani_threshold}% ANI "
        f"({n_merged} clusters merge >1 genome)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    sys.exit(main())
