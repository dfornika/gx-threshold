//
// Read QC + dehosting: everything that runs on every sample before any
// species-ID/library-type analysis touches the reads.
//
// FASTQC/FASTP/FASTPLONG each run twice - once on the raw input reads, and
// again on the final (trimmed + dehosted) reads that flow to every
// downstream analysis - via module aliasing (FASTQC_RAW/FASTQC_FINAL etc.),
// the standard nf-core idiom for invoking the same process twice with
// different roles in one workflow (see e.g. nf-core/rnaseq's
// FASTQC_RAW/FASTQC_TRIM). There's no cleaner alternative to aliasing for
// this - it's the recommended pattern, not a workaround.
//
// The final pass is measurement-only by construction, not by convention:
// FASTQC never modifies reads, so re-running it on the final reads is
// inherently safe. FASTP/FASTPLONG do modify reads by design (that's the
// point of the first pass), so the final pass uses their built-in
// `discard_trimmed_pass: true` (a real, tested, documented nf-core module
// option - "use fastp for the output report only") - no trimmed-reads
// output is ever written, so there is nothing from this pass that could
// accidentally be published or picked up downstream, regardless of what
// fastp/fastplong would otherwise have changed.
//

include { FASTQC as FASTQC_RAW     } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_FINAL   } from '../../modules/nf-core/fastqc/main'
include { FASTP as FASTP_TRIM      } from '../../modules/nf-core/fastp/main'
include { FASTP as FASTP_FINAL_QC  } from '../../modules/nf-core/fastp/main'
include { FASTPLONG as FASTPLONG_TRIM     } from '../../modules/nf-core/fastplong/main'
include { FASTPLONG as FASTPLONG_FINAL_QC } from '../../modules/nf-core/fastplong/main'
include { DEHOST                  } from '../../modules/local/dehost/main'

workflow READ_QC_AND_DEHOSTING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    outdir

    main:

    def ch_multiqc_files = channel.empty()

    //
    // Tag each sample's meta with its platform. `single_end` samples (no
    // fastq_2) are assumed to be Nanopore; paired-end samples are assumed to
    // be Illumina. This is the only place that assumption is encoded -
    // everything downstream reads `meta.platform` instead.
    //
    def ch_reads = ch_samplesheet.map { meta, reads ->
        tuple(meta + [platform: meta.single_end ? 'nanopore' : 'illumina'], reads)
    }

    //
    // Anchor row list for SAMPLE_SUMMARY - every sample gets a row regardless
    // of which optional stages ran later.
    //
    def ch_sample_manifest = ch_reads
        .map { meta, _reads -> "${meta.id}\t${meta.platform}" }
        .collectFile(
            name: 'sample_manifest.tsv',
            storeDir: "${outdir}/pipeline_info",
            sort: true,
            seed: "sample\tplatform\n",
            newLine: true
        )

    //
    // MODULE: Run FastQC on the raw input reads.
    //
    FASTQC_RAW(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.map{ _meta, file -> file })

    //
    // Split reads by platform for platform-specific trimming/QC tools.
    //
    def ch_reads_by_platform = ch_reads.branch { meta, _reads ->
        long_reads:  meta.platform == 'nanopore'
        short_reads: meta.platform == 'illumina'
    }

    //
    // MODULE: Run fastp on Illumina paired-end reads
    //
    FASTP_TRIM(ch_reads_by_platform.short_reads.map { meta, reads -> tuple(meta, reads, []) }, false, false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP_TRIM.out.json.map{ _meta, file -> file })

    //
    // MODULE: Run fastplong on Nanopore single-end reads
    //
    FASTPLONG_TRIM(ch_reads_by_platform.long_reads, [], false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG_TRIM.out.json.map{ _meta, file -> file })

    //
    // Recombine platform-specific trimmed reads into single channels.
    //
    def ch_trimmed   = FASTP_TRIM.out.reads.mix(FASTPLONG_TRIM.out.reads)
    def ch_trim_json = FASTP_TRIM.out.json.mix(FASTPLONG_TRIM.out.json)

    //
    // MODULE: Dehost - remove host (human) reads by aligning to a host reference
    // and keeping the unmapped reads. This runs early so that every downstream
    // step (and any shared reads) is host-depleted. Toggle with --skip_dehosting.
    //
    def ch_clean_reads    = ch_trimmed
    def ch_dehost_summary = channel.value([])
    if (!params.skip_dehosting) {
        if (!params.dehost_reference) {
            error("Dehosting is enabled but --dehost_reference was not set. Provide a host reference (FASTA or minimap2 .mmi) or run with --skip_dehosting.")
        }
        def ch_host_reference = channel.value(file(params.dehost_reference, checkIfExists: true))
        DEHOST(ch_trimmed, ch_host_reference, params.dehost_scrub_headers)
        ch_clean_reads = DEHOST.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(DEHOST.out.stats.map { _meta, file -> file })

        // .ifEmpty([]): see subworkflows/local/species_id.nf for why -
        // collectFile emits nothing at all if its input is completely empty.
        ch_dehost_summary = DEHOST.out.stats
            .map { _meta, file -> file }
            .collectFile(
                name: 'dehost_summary.tsv',
                storeDir: "${outdir}/dehost",
                sort: true,
                seed: "sample\tplatform\tinput_reads\thost_reads\tdehosted_reads\tpercent_host\n"
            )
            .ifEmpty([])
    }

    //
    // Final QC pass on the reads that flow to every downstream analysis -
    // see the measurement-only note above.
    //
    FASTQC_FINAL(ch_clean_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_FINAL.out.zip.map{ _meta, file -> file })

    def ch_clean_reads_by_platform = ch_clean_reads.branch { meta, _reads ->
        long_reads:  meta.platform == 'nanopore'
        short_reads: meta.platform == 'illumina'
    }

    FASTP_FINAL_QC(ch_clean_reads_by_platform.short_reads.map { meta, reads -> tuple(meta, reads, []) }, true, false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP_FINAL_QC.out.json.map{ _meta, file -> file })

    FASTPLONG_FINAL_QC(ch_clean_reads_by_platform.long_reads, [], true, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG_FINAL_QC.out.json.map{ _meta, file -> file })

    emit:
    reads          = ch_clean_reads     // tuple(meta, reads) - trimmed + dehosted, ready for downstream analysis
    trim_json      = ch_trim_json       // tuple(meta, json)  - pre-dehost fastp/fastplong QC json
    multiqc_files  = ch_multiqc_files
    sample_manifest = ch_sample_manifest // path - sample, platform - anchor row list for SAMPLE_SUMMARY
    dehost_summary  = ch_dehost_summary  // path (or []) - for SAMPLE_SUMMARY
}
