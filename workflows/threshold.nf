/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { FASTPLONG              } from '../modules/nf-core/fastplong/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { MASH_DIST              } from '../modules/nf-core/mash/dist/main'
include { SOURMASH_SKETCH        } from '../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_GATHER        } from '../modules/nf-core/sourmash/gather/main'
include { SYLPH_PROFILE          } from '../modules/nf-core/sylph/profile/main'
include { DEHOST                 } from '../modules/local/dehost/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_threshold_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow THRESHOLD {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
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
    // MODULE: Run FastQC
    //
    FASTQC(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

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
    FASTP(ch_reads_by_platform.short_reads.map { meta, reads -> tuple(meta, reads, []) }, false, false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map{ _meta, file -> file })

    //
    // MODULE: Run fastplong on Nanopore single-end reads
    //
    FASTPLONG(ch_reads_by_platform.long_reads, [], false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG.out.json.map{ _meta, file -> file })

    //
    // Recombine platform-specific trimmed reads into a single channel.
    //
    def ch_trimmed = FASTP.out.reads.mix(FASTPLONG.out.reads)

    //
    // MODULE: Dehost - remove host (human) reads by aligning to a host reference
    // and keeping the unmapped reads. This runs early so that every downstream
    // step (and any shared reads) is host-depleted. Toggle with --skip_dehosting.
    //
    def ch_clean_reads = ch_trimmed
    if (!params.skip_dehosting) {
        if (!params.dehost_reference) {
            error("Dehosting is enabled but --dehost_reference was not set. Provide a host reference (FASTA or minimap2 .mmi) or run with --skip_dehosting.")
        }
        def ch_host_reference = channel.value(file(params.dehost_reference, checkIfExists: true))
        DEHOST(ch_trimmed, ch_host_reference, params.dehost_scrub_headers)
        ch_clean_reads = DEHOST.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(DEHOST.out.stats.map { _meta, file -> file })
    }

    //
    // MODULE: Mash - species identification via k-mer distance to a panel of
    // reference genome sketches. One of three species-ID tools currently under
    // evaluation (mash/sourmash/sylph) - toggle with --skip_mash.
    //
    if (!params.skip_mash) {
        if (!params.mash_db) {
            error("Mash species ID is enabled but --mash_db was not set. Provide a Mash sketch (.msh) or run with --skip_mash.")
        }
        def ch_mash_db = channel.value(file(params.mash_db, checkIfExists: true))
        MASH_DIST(ch_clean_reads, ch_mash_db)
        ch_multiqc_files = ch_multiqc_files.mix(MASH_DIST.out.dist.map { _meta, file -> file })
    }

    //
    // MODULE: Sourmash - species identification via FracMinHash containment
    // (sketch reads, then gather against the reference panel). Second of
    // three species-ID tools under evaluation - toggle with --skip_sourmash.
    //
    if (!params.skip_sourmash) {
        if (!params.sourmash_db) {
            error("Sourmash species ID is enabled but --sourmash_db was not set. Provide a sourmash signature collection (.sig/.sig.zip) or run with --skip_sourmash.")
        }
        def ch_sourmash_db = channel.value(file(params.sourmash_db, checkIfExists: true))
        SOURMASH_SKETCH(ch_clean_reads)
        SOURMASH_GATHER(SOURMASH_SKETCH.out.signatures, ch_sourmash_db, false, false, false, false)
        ch_multiqc_files = ch_multiqc_files.mix(SOURMASH_GATHER.out.result.map { _meta, file -> file })
    }

    //
    // MODULE: Sylph - species identification via containment ANI. Sketches
    // and profiles reads against the reference panel in one step (no
    // separate sketch stage needed). Third of three species-ID tools under
    // evaluation - toggle with --skip_sylph.
    //
    if (!params.skip_sylph) {
        if (!params.sylph_db) {
            error("Sylph species ID is enabled but --sylph_db was not set. Provide a Sylph genome database (.syldb) or run with --skip_sylph.")
        }
        def ch_sylph_db = channel.value(file(params.sylph_db, checkIfExists: true))
        SYLPH_PROFILE(ch_clean_reads, ch_sylph_db)
        ch_multiqc_files = ch_multiqc_files.mix(SYLPH_PROFILE.out.profile_out.map { _meta, file -> file })
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'gx-threshold_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'gx-threshold'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
