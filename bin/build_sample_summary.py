#!/usr/bin/env python3
"""
Build one CSV with one row per sample, pulling the main verdict/metric from
each per-stage summary TSV this pipeline already produces. Not a replacement
for those TSVs (which keep every field for each stage) - this is a single
at-a-glance table across stages, so a reviewer doesn't need to open several
different files to see what happened to one sample.

Every stage input is optional: if a stage was skipped (or never reached a
given sample, e.g. 16S detection only runs for composition-inconclusive
samples), its columns are just "NA" for that row rather than the row being
dropped - `--samples` (sample, platform) is the anchor and always determines
which rows exist.
"""

import argparse
import csv
import sys

NA = "NA"

HEADER = [
    "sample",
    "platform",
    "dehost_percent_host",
    "library_type",
    "library_type_consensus_method",
    "library_type_largest_cluster_frac",
    "library_type_fastp_verdict",
    "library_type_fastp_duplication_rate",
    "library_type_aligned_verdict",
    "library_type_aligned_index_of_dispersion",
    "library_type_pileup_verdict",
    "library_type_pileup_piled_frac",
    "species_id_mash_organism",
    "species_id_sourmash_organism",
    "species_id_sylph_organism",
    "species_id_consensus_accession",
    "species_id_consensus_species_name",
    "species_id_consensus_method",
    "composition",
    "composition_adjusted_effective_n",
    "reference_accession",
    "reference_species_name",
    "sixteen_s",
    "sixteen_s_passed_frac",
    "assembler",
    "assembly_n_contigs",
    "assembly_total_length",
    "assembly_n50",
    "assembly_completeness",
    "assembly_contamination",
]


def load_rows(path):
    if not path:
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def index_by_sample(rows):
    return {row["sample"]: row for row in rows}


def index_species_id_by_sample(rows):
    """species_id_summary.tsv has one row per sample *per tool* - pivot to
    {sample: {tool: organism}}."""
    by_sample = {}
    for row in rows:
        by_sample.setdefault(row["sample"], {})[row["tool"]] = row["organism"]
    return by_sample


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples", required=True, help="sample, platform - the anchor row list")
    ap.add_argument("--dehost")
    ap.add_argument("--library-type")
    ap.add_argument("--library-type-cluster")
    ap.add_argument("--library-type-aligned")
    ap.add_argument("--library-type-pileup")
    ap.add_argument("--library-type-consensus")
    ap.add_argument("--species-id")
    ap.add_argument("--species-composition")
    ap.add_argument("--species-id-consensus")
    ap.add_argument("--sixteen-s")
    ap.add_argument("--assembly")
    args = ap.parse_args()

    samples = load_rows(args.samples)
    dehost = index_by_sample(load_rows(args.dehost))
    library_type = index_by_sample(load_rows(args.library_type))
    library_type_cluster = index_by_sample(load_rows(args.library_type_cluster))
    library_type_aligned = index_by_sample(load_rows(args.library_type_aligned))
    library_type_pileup = index_by_sample(load_rows(args.library_type_pileup))
    library_type_consensus = index_by_sample(load_rows(args.library_type_consensus))
    species_id = index_species_id_by_sample(load_rows(args.species_id))
    species_composition = index_by_sample(load_rows(args.species_composition))
    species_id_consensus = index_by_sample(load_rows(args.species_id_consensus))
    sixteen_s = index_by_sample(load_rows(args.sixteen_s))
    assembly = index_by_sample(load_rows(args.assembly))

    writer = csv.DictWriter(sys.stdout, fieldnames=HEADER)
    writer.writeheader()

    for row in samples:
        sample = row["sample"]
        d = dehost.get(sample, {})
        lt = library_type.get(sample, {})
        ltc = library_type_cluster.get(sample, {})
        lta = library_type_aligned.get(sample, {})
        ltp = library_type_pileup.get(sample, {})
        ltcon = library_type_consensus.get(sample, {})
        sid = species_id.get(sample, {})
        sc = species_composition.get(sample, {})
        sidcon = species_id_consensus.get(sample, {})
        s16 = sixteen_s.get(sample, {})
        asm = assembly.get(sample, {})

        writer.writerow(
            {
                "sample": sample,
                "platform": row["platform"],
                "dehost_percent_host": d.get("percent_host", NA),
                "library_type": ltcon.get("verdict", NA),
                "library_type_consensus_method": ltcon.get("method", NA),
                "library_type_largest_cluster_frac": ltc.get("largest_cluster_frac", NA),
                "library_type_fastp_verdict": lt.get("verdict", NA),
                "library_type_fastp_duplication_rate": lt.get("duplication_rate", NA),
                "library_type_aligned_verdict": lta.get("verdict", NA),
                "library_type_aligned_index_of_dispersion": lta.get("index_of_dispersion", NA),
                "library_type_pileup_verdict": ltp.get("verdict", NA),
                "library_type_pileup_piled_frac": ltp.get("piled_frac", NA),
                "species_id_mash_organism": sid.get("mash", NA),
                "species_id_sourmash_organism": sid.get("sourmash", NA),
                "species_id_sylph_organism": sid.get("sylph", NA),
                "species_id_consensus_accession": sidcon.get("accession", NA),
                "species_id_consensus_species_name": sidcon.get("species_name", NA),
                "species_id_consensus_method": sidcon.get("method", NA),
                "composition": sc.get("verdict", NA),
                "composition_adjusted_effective_n": sc.get("adjusted_effective_n", NA),
                "reference_accession": sidcon.get("accession", NA),
                "reference_species_name": sidcon.get("species_name", NA),
                "sixteen_s": s16.get("verdict", NA),
                "sixteen_s_passed_frac": s16.get("passed_frac", NA),
                "assembler": asm.get("assembler", NA),
                "assembly_n_contigs": asm.get("n_contigs", NA),
                "assembly_total_length": asm.get("total_length", NA),
                "assembly_n50": asm.get("n50", NA),
                "assembly_completeness": asm.get("completeness", NA),
                "assembly_contamination": asm.get("contamination", NA),
            }
        )


if __name__ == "__main__":
    main()
