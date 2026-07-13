<p align="center">
  <img src="https://raw.githubusercontent.com/dfornika/gx-threshold/refs/heads/main/docs/images/threshold.png" width="300" alt="Accessibility Description" />
</p>

## Introduction

**gx-threshold** is a QC and routine-analysis pipeline for microbial genomics,
intended as a standard "entrypoint" that runs on every sequenced library as it
arrives in a genomics analysis platform. It works on both **Illumina** and
**Oxford Nanopore** reads throughout, takes a simple samplesheet of FASTQ
files, and produces a per-sample summary of what each library is and what's in
it. Results from every stage are collated by [`MultiQC`](http://multiqc.info/)
and into a single `sample_summary.csv`.

The pipeline is an active work in progress. Some stages deliberately run
**several independent methods for the same question in parallel** and fuse them
into a consensus verdict, so the individual methods can keep being compared
against real data before any are simplified away. See
[`docs/testing.md`](docs/testing.md) for the detailed, method-by-method
rationale, validation status, and known caveats.

### Analyses

1. **Read QC & trimming** — [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
   plus [`fastp`](https://github.com/OpenGene/fastp) (Illumina) /
   [`fastplong`](https://github.com/OpenGene/fastp) (Nanopore), run on both the
   raw reads and again on the final trimmed/dehosted reads that flow downstream.
2. **Host read removal (dehosting)** — [`minimap2`](https://github.com/lh3/minimap2)
   against a host reference (`--dehost_reference`); toggle with `--skip_dehosting`.
3. **Library type: amplicon vs. shotgun** — four independent detectors
   (fastp QC-JSON insert-size/duplication, reference-free read clustering,
   reference-aligned depth dispersion, and aligned read-position pileup) fused
   into one majority-vote **consensus** verdict (`meta.library_type`).
4. **Species identification** — [`mash`](https://github.com/marbl/Mash),
   [`sourmash`](https://github.com/sourmash-bio/sourmash), and
   [`sylph`](https://github.com/bluenote-1577/sylph) run in parallel against a
   reference-genome database, reconciled into a single consensus call via a
   majority vote on species taxid (each tool toggleable via
   `--skip_mash`/`--skip_sourmash`/`--skip_sylph`).
5. **Pure culture vs. metagenomic** — composition-breadth analysis of the
   sourmash results (how many reference genomes it takes to explain the sample),
   with ANI-based collapsing of nomenclature artifacts; low database coverage
   is reported as `inconclusive` rather than a false confident call.
6. **Reference genome fetch & cache** — fetches the consensus species' assembly
   from NCBI (via [`entrez-direct`](https://www.ncbi.nlm.nih.gov/books/NBK179288/)),
   cached by accession across runs. On by default; a fetch failure degrades
   gracefully rather than failing the run. Toggle with
   `--skip_reference_genome_fetch` / `--reference_genome_cache_dir`.
7. **16S rRNA amplicon detection** — for samples whose composition breadth came
   back inconclusive, an alignment-based check against a dedicated 16S database
   to distinguish 16S/marker-gene libraries from other cases. Off by default
   (needs `--sixteen_s_db`).

Planned but not yet implemented: draft-genome assembly + annotation, and
assembly-level QC / contamination detection for pure-culture shotgun libraries.

## Usage

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

```bash
nextflow run BCCDC-PHL/gx-threshold \
   -profile <conda/docker/singularity/apptainer> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use BCCDC-PHL/gx-threshold for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
