#!/usr/bin/env python3
"""
Normalise one species-ID tool's raw output for one sample into a single summary
row, so mash/sourmash/sylph calls can be compared side by side. Each tool reports
confidence on its own native scale (mash: distance, lower is better; sourmash:
fraction of the matched reference recovered, higher is better; sylph: adjusted
ANI, higher is better) - this deliberately doesn't try to normalise those onto one
scale, just surfaces each tool's own top call and its own metric.
"""

import argparse
import csv
import gzip
import os
import sys


def load_manifest(path):
    if not path or not os.path.isfile(path):
        return {}
    with open(path) as fh:
        return {row["accession"]: row for row in csv.DictReader(fh)}


def accession_from_path(path):
    return os.path.basename(path).removesuffix(".fasta")


def parse_mash(path):
    """Tab-separated, no header: ref-ID, query-ID, distance, p-value, shared-hashes.
    Multiple rows per reference (one per mate file for paired-end) - take the
    single best (lowest-distance) row across the whole file."""
    best = None
    with open(path) as fh:
        for line in fh:
            ref, _query, distance, _pvalue, _shared = line.rstrip("\n").split("\t")
            distance = float(distance)
            if best is None or distance < best[1]:
                best = (ref, distance)
    if best is None:
        return None
    accession = accession_from_path(best[0])
    return accession, "distance", f"{best[1]:.6g}"


def parse_sourmash(path):
    """Gzipped CSV from `sourmash gather`, already rank-ordered - first row is
    the best hit. `name` is "<accession> <organism>" (set at DB-build time)."""
    with gzip.open(path, "rt") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return None
    top = rows[0]
    accession = top["name"].split(" ", 1)[0]
    return accession, "f_match_orig", f"{float(top['f_match_orig']):.6g}"


def parse_sylph(path):
    """TSV from `sylph profile` - usually one row (only genomes clearing sylph's
    own containment/ANI threshold are reported); take the highest-abundance row."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return None
    top = max(rows, key=lambda r: float(r["Taxonomic_abundance"]))
    accession = accession_from_path(top["Genome_file"])
    return accession, "adjusted_ani", top["Adjusted_ANI"]


PARSERS = {"mash": parse_mash, "sourmash": parse_sourmash, "sylph": parse_sylph}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tool", required=True, choices=sorted(PARSERS))
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("--manifest", default=None, help="species_db manifest.csv (accession -> organism/taxonomy); omit to report raw accessions and NA taxonomy")
    ap.add_argument("result", help="the tool's raw output file for this sample")
    args = ap.parse_args()

    manifest = load_manifest(args.manifest)
    parsed = PARSERS[args.tool](args.result)

    if parsed is None:
        print(f"{args.sample}\t{args.platform}\t{args.tool}\tNA\tno hit\tNA\tNA\tNA\tNA")
        return

    accession, metric_name, metric_value = parsed
    entry = manifest.get(accession, {})
    organism = entry.get("organism", accession)
    species_taxid = entry.get("species_taxid", "NA")
    species_name = entry.get("species_name", "NA")
    print(
        f"{args.sample}\t{args.platform}\t{args.tool}\t{accession}\t{organism}\t"
        f"{species_taxid}\t{species_name}\t{metric_name}\t{metric_value}"
    )


if __name__ == "__main__":
    sys.exit(main())
