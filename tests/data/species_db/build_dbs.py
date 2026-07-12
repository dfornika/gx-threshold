#!/usr/bin/env python3
"""
Build tiny mash/sourmash/sylph reference databases from a small panel of real
reference genomes, for evaluating species-ID tools against the known-species
samples in tests/data/real/ (KPNEUMONIAE_WGS, ECOLI_WGS, SARS2_AMPLICON_*).

Panel: manifest.csv - 3 "known" genomes matching those samples' true species,
plus ~7 "decoy" genomes (a mix of easy/distant-genus and hard/same-genus cases)
so a correct call is a real discrimination test, not a trivial one-entry lookup.

Identifiers are NCBI *assembly* accessions (GCF_.../GCA_...), matching what the
official/canonical mash RefSeq sketch, sourmash GTDB databases, and sylph GTDB
databases all use - not plain nucleotide accessions (which is what this script
used before; see git history). `manifest.csv`'s `source_nucleotide_accession`
column records the single representative record originally used to pick each
organism, resolved to its assembly via `elink` (nuccore -> assembly).

Requires (not installed by this script):
  - ncbi-client (https://github.com/dfornika/ncbi-client-py) for genome download
  - mash, sourmash, sylph CLIs on PATH (see the `species-db` conda env)

Downloads each assembly's whole-genome FASTA (all contigs/plasmids in one file)
into a gitignored genomes_src/ scratch dir (not committed - see README.md to
regenerate), then builds the three committed DB artifacts: mash/species_mini.msh,
sourmash/species_mini.sig.zip (+ sourmash/taxonomy.csv), sylph/species_mini.syldb.
Each tool sketches a whole file as one entry by default, so this is the only
part of the pipeline that needed to change - the sketch-building functions
below are accession-format-agnostic.

Regenerate with: python3 tests/data/species_db/build_dbs.py
"""

import csv
import subprocess
import zipfile
from pathlib import Path

from ncbi_client import NCBIClient

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "manifest.csv"
GENOMES_SRC = HERE / "genomes_src"
MASH_DIR = HERE / "mash"
SOURMASH_DIR = HERE / "sourmash"
SYLPH_DIR = HERE / "sylph"


def read_manifest():
    with MANIFEST.open() as fh:
        return list(csv.DictReader(fh))


def download_genomes(client, rows):
    GENOMES_SRC.mkdir(exist_ok=True)
    paths = []
    for row in rows:
        acc = row["accession"]
        dest = GENOMES_SRC / f"{acc}.fasta"
        if dest.exists():
            print(f"  {acc}: already downloaded, skipping")
            paths.append(dest)
            continue
        print(f"  downloading {acc} ({row['organism']})")
        zip_path = GENOMES_SRC / f"{acc}.zip"
        client.download_genome(acc, zip_path)
        with zipfile.ZipFile(zip_path) as zf:
            fna_name = next(n for n in zf.namelist() if n.endswith("_genomic.fna"))
            with zf.open(fna_name) as src, dest.open("wb") as out:
                out.write(src.read())
        zip_path.unlink()
        paths.append(dest)
    return paths


def build_mash(genome_paths):
    MASH_DIR.mkdir(exist_ok=True)
    out_prefix = MASH_DIR / "species_mini"
    subprocess.run(
        ["mash", "sketch", "-o", str(out_prefix), *[str(p) for p in genome_paths]],
        check=True,
    )
    print(f"  wrote {out_prefix}.msh")


def build_sourmash(rows, genome_paths_by_acc):
    SOURMASH_DIR.mkdir(exist_ok=True)
    sigs = []
    for row in rows:
        acc = row["accession"]
        name = f"{acc} {row['organism']}"
        sig_path = SOURMASH_DIR / f"{acc}.sig"
        subprocess.run(
            [
                "sourmash", "sketch", "dna",
                "-p", "k=21,scaled=1000",
                "--name", name,
                "-o", str(sig_path),
                str(genome_paths_by_acc[acc]),
            ],
            check=True,
        )
        sigs.append(sig_path)

    combined = SOURMASH_DIR / "species_mini.sig.zip"
    combined.unlink(missing_ok=True)
    subprocess.run(
        ["sourmash", "sig", "cat", *[str(s) for s in sigs], "-o", str(combined)],
        check=True,
    )
    for s in sigs:
        s.unlink()
    print(f"  wrote {combined}")

    taxonomy_path = SOURMASH_DIR / "taxonomy.csv"
    with taxonomy_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ident", "organism", "strain", "category"])
        for row in rows:
            writer.writerow([row["accession"], row["organism"], row["strain"], row["category"]])
    print(f"  wrote {taxonomy_path}")


def build_sylph(genome_paths):
    SYLPH_DIR.mkdir(exist_ok=True)
    subprocess.run(
        [
            "sylph", "sketch",
            "-g", *[str(p) for p in genome_paths],
            "-o", "species_mini",
        ],
        check=True,
        cwd=SYLPH_DIR,
    )
    print(f"  wrote {SYLPH_DIR}/species_mini.syldb")


def main():
    rows = read_manifest()

    print(f"Downloading {len(rows)} reference genomes:")
    client = NCBIClient()
    genome_paths = download_genomes(client, rows)
    genome_paths_by_acc = {row["accession"]: p for row, p in zip(rows, genome_paths)}

    print("\nBuilding mash database:")
    build_mash(genome_paths)

    print("\nBuilding sourmash database:")
    build_sourmash(rows, genome_paths_by_acc)

    print("\nBuilding sylph database:")
    build_sylph(genome_paths)


if __name__ == "__main__":
    main()
