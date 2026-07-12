#!/usr/bin/env python3
"""
Classify a sample's library prep as amplicon vs. shotgun by clustering reads
against each other - no reference genome needed, unlike
classify_library_type_aligned.awk.

Rationale: amplicon libraries repeatedly re-sequence the same small set of
PCR products, so a large fraction of reads are near-identical to many other
reads in the same run. Shotgun libraries fragment DNA close to randomly, so
at realistic depth, two reads being near-identical across most of their
length is rare. This sidesteps needing a reference genome at all (and so
also sidesteps the failure mode found while validating
classify_library_type_aligned.awk, where real bacterial genome repeat
structure - not the library prep - could inflate a reference-alignment-based
statistic).

Input is a PAF file from an all-vs-all minimap2 alignment (reads against
themselves - see modules/local/read_overlap/). Two reads are linked into the
same cluster if their overlap passes an identity and coverage threshold
(read length and error profile differ enough between platforms that these
are tunable per platform - Nanopore needs a more permissive identity
threshold). Clusters are then just connected components under those edges
(union-find), and the summary statistic is what fraction of all reads ended
up clustered with at least one other read, plus how concentrated the
resulting cluster-size distribution is (effective number of clusters, via
the inverse Simpson index - low means a few big clusters dominate, which is
what a small set of PCR targets would produce; high means clusters stay
small/singleton, consistent with random shotgun sampling).

This is a first-pass implementation, deliberately not tuned to a single
"final" threshold - see docs/testing.md for validation status and known
caveats.
"""

import argparse
import sys


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


def cluster_reads(paf_path, identity_threshold, coverage_threshold):
    uf = UnionFind()
    seen_reads = set()

    with open(paf_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, qlen, _qstart, _qend, _strand, tname, tlen = fields[0:7]
            nmatch, alnlen = fields[9:11]
            qlen, tlen, nmatch, alnlen = int(qlen), int(tlen), int(nmatch), int(alnlen)

            seen_reads.add(qname)
            seen_reads.add(tname)
            if qname == tname or alnlen == 0:
                continue

            identity = nmatch / alnlen
            coverage = alnlen / min(qlen, tlen)
            if identity >= identity_threshold and coverage >= coverage_threshold:
                uf.union(qname, tname)

    return uf, seen_reads


def summarize(uf, seen_reads, total_reads):
    component_sizes = {}
    for read in seen_reads:
        root = uf.find(read)
        component_sizes[root] = component_sizes.get(root, 0) + 1

    # Reads never mentioned in the PAF at all (shouldn't normally happen,
    # since every read aligns trivially to itself) still count as their own
    # singleton cluster.
    n_unseen = max(0, total_reads - len(seen_reads))
    sizes = list(component_sizes.values()) + [1] * n_unseen

    n_singletons = sum(1 for s in sizes if s == 1)
    largest_cluster = max(sizes) if sizes else 0
    largest_cluster_frac = largest_cluster / total_reads if total_reads else 0.0
    clustered_frac = 1 - (n_singletons / total_reads) if total_reads else 0.0

    sum_sq_frac = sum((s / total_reads) ** 2 for s in sizes) if total_reads else 0.0
    effective_clusters = (1 / sum_sq_frac) if sum_sq_frac > 0 else 0.0
    effective_clusters_frac = (effective_clusters / total_reads) if total_reads else 0.0

    return {
        "largest_cluster_frac": largest_cluster_frac,
        "clustered_frac": clustered_frac,
        "effective_clusters_frac": effective_clusters_frac,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("--total-reads", type=int, required=True, help="Total read count in the input (from the FASTQ, not just what appears in the PAF)")
    ap.add_argument("--identity-threshold", type=float, required=True)
    ap.add_argument("--coverage-threshold", type=float, required=True)
    ap.add_argument("--effective-clusters-frac-threshold", type=float, default=0.5, help="Below this, call amplicon")
    ap.add_argument("paf", help="All-vs-all minimap2 PAF file")
    args = ap.parse_args()

    if args.total_reads == 0:
        print(f"{args.sample}\t{args.platform}\tinconclusive\t0\tNA\tNA\tNA\tno reads")
        return

    uf, seen_reads = cluster_reads(args.paf, args.identity_threshold, args.coverage_threshold)
    stats = summarize(uf, seen_reads, args.total_reads)

    verdict = "amplicon" if stats["effective_clusters_frac"] < args.effective_clusters_frac_threshold else "shotgun"

    print(
        f"{args.sample}\t{args.platform}\t{verdict}\t{args.total_reads}\t"
        f"{stats['largest_cluster_frac']:.4g}\t{stats['clustered_frac']:.4g}\t"
        f"{stats['effective_clusters_frac']:.4g}\t-"
    )


if __name__ == "__main__":
    sys.exit(main())
