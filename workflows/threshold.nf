/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { FASTPLONG              } from '../modules/nf-core/fastplong/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { DEHOST                 } from '../modules/local/dehost/main'
include { LIBRARY_TYPE_REFERENCE_FREE   } from '../subworkflows/local/library_type_reference_free'
include { SPECIES_ID                    } from '../subworkflows/local/species_id'
include { SPECIES_COMPOSITION_ANALYSIS  } from '../subworkflows/local/species_composition_analysis'
include { REFERENCE_GENOME              } from '../subworkflows/local/reference_genome'
include { ALIGNMENT_BASED_LIBRARY_TYPE  } from '../subworkflows/local/alignment_based_library_type'
include { SIXTEEN_S_DETECTION           } from '../subworkflows/local/sixteen_s_detection'
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
    // Recombine platform-specific trimmed reads into single channels.
    //
    def ch_trimmed   = FASTP.out.reads.mix(FASTPLONG.out.reads)
    def ch_trim_json = FASTP.out.json.mix(FASTPLONG.out.json)

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
    // SUBWORKFLOW: Library type, reference-free - the two amplicon-vs-shotgun
    // approaches that need no fetched reference genome (LIBRARY_TYPE,
    // LIBRARY_TYPE_CLUSTER). Tags meta.library_type from the read-clustering
    // call - see docs/testing.md and subworkflows/local/library_type_reference_free.nf.
    //
    LIBRARY_TYPE_REFERENCE_FREE(ch_trim_json, ch_clean_reads, outdir)
    ch_clean_reads   = LIBRARY_TYPE_REFERENCE_FREE.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(LIBRARY_TYPE_REFERENCE_FREE.out.multiqc_files)

    //
    // SUBWORKFLOW: Species identification - mash/sourmash/sylph, evaluated
    // in parallel; see tests/data/species_db/README.md.
    //
    def ch_species_id_manifest = params.species_id_manifest ? file(params.species_id_manifest, checkIfExists: true) : []
    SPECIES_ID(ch_clean_reads, ch_species_id_manifest, outdir)
    ch_multiqc_files = ch_multiqc_files.mix(SPECIES_ID.out.multiqc_files)

    //
    // SUBWORKFLOW: Pure-culture vs. metagenomic detection from species-ID
    // composition breadth. Tags meta.composition - see docs/testing.md and
    // subworkflows/local/species_composition_analysis.nf.
    //
    SPECIES_COMPOSITION_ANALYSIS(ch_clean_reads, SPECIES_ID.out.sourmash_gather_result, ch_species_id_manifest, outdir)
    ch_clean_reads   = SPECIES_COMPOSITION_ANALYSIS.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(SPECIES_COMPOSITION_ANALYSIS.out.multiqc_files)

    //
    // SUBWORKFLOW: Reference genome selection + fetch/cache - see
    // docs/testing.md and subworkflows/local/reference_genome.nf.
    //
    REFERENCE_GENOME(SPECIES_ID.out.species_id_rows, outdir)
    ch_multiqc_files = ch_multiqc_files.mix(REFERENCE_GENOME.out.multiqc_files)

    //
    // SUBWORKFLOW: Library type, alignment-based - the two approaches that
    // need the fetched reference genome (LIBRARY_TYPE_ALIGNED,
    // LIBRARY_TYPE_PILEUP) - see docs/testing.md and
    // subworkflows/local/alignment_based_library_type.nf.
    //
    ALIGNMENT_BASED_LIBRARY_TYPE(ch_clean_reads, REFERENCE_GENOME.out.sample_reference, outdir)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT_BASED_LIBRARY_TYPE.out.multiqc_files)

    //
    // SUBWORKFLOW: 16S rRNA amplicon detection - only for samples whose
    // composition breadth came back inconclusive (amplicon or not); see
    // docs/testing.md and subworkflows/local/sixteen_s_detection.nf.
    //
    SIXTEEN_S_DETECTION(ch_clean_reads, outdir)
    ch_clean_reads   = SIXTEEN_S_DETECTION.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(SIXTEEN_S_DETECTION.out.multiqc_files)

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
