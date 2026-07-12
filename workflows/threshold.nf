/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                       } from '../modules/nf-core/multiqc/main'
include { READ_QC_AND_DEHOSTING         } from '../subworkflows/local/read_qc_and_dehosting'
include { LIBRARY_TYPE_REFERENCE_FREE   } from '../subworkflows/local/library_type_reference_free'
include { SPECIES_ID                    } from '../subworkflows/local/species_id'
include { SPECIES_COMPOSITION_ANALYSIS  } from '../subworkflows/local/species_composition_analysis'
include { REFERENCE_GENOME              } from '../subworkflows/local/reference_genome'
include { ALIGNMENT_BASED_LIBRARY_TYPE  } from '../subworkflows/local/alignment_based_library_type'
include { SIXTEEN_S_DETECTION           } from '../subworkflows/local/sixteen_s_detection'
include { SAMPLE_SUMMARY                } from '../modules/local/sample_summary/main'
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
    // SUBWORKFLOW: Read QC + dehosting - FastQC/fastp/fastplong on the raw
    // reads, dehosting, then a measurement-only FastQC/fastp/fastplong pass
    // on the final reads that flow to every downstream analysis - see
    // docs/testing.md and subworkflows/local/read_qc_and_dehosting.nf.
    //
    READ_QC_AND_DEHOSTING(ch_samplesheet, outdir)
    def ch_clean_reads = READ_QC_AND_DEHOSTING.out.reads
    def ch_trim_json   = READ_QC_AND_DEHOSTING.out.trim_json
    ch_multiqc_files   = ch_multiqc_files.mix(READ_QC_AND_DEHOSTING.out.multiqc_files)

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
    // MODULE: One-row-per-sample summary CSV, pulling the main verdict/metric
    // from each per-stage summary TSV above - see docs/testing.md and
    // modules/local/sample_summary/main.nf. Every stage input is optional
    // (path or [] if that stage was skipped); the sample manifest from
    // READ_QC_AND_DEHOSTING is the only required one, and anchors which rows
    // exist regardless of which optional stages ran.
    //
    SAMPLE_SUMMARY(
        READ_QC_AND_DEHOSTING.out.sample_manifest,
        READ_QC_AND_DEHOSTING.out.dehost_summary,
        LIBRARY_TYPE_REFERENCE_FREE.out.library_type_summary,
        LIBRARY_TYPE_REFERENCE_FREE.out.library_type_cluster_summary,
        ALIGNMENT_BASED_LIBRARY_TYPE.out.library_type_aligned_summary,
        ALIGNMENT_BASED_LIBRARY_TYPE.out.library_type_pileup_summary,
        SPECIES_ID.out.species_id_summary,
        SPECIES_COMPOSITION_ANALYSIS.out.species_composition_summary,
        REFERENCE_GENOME.out.reference_selection_summary,
        SIXTEEN_S_DETECTION.out.sixteen_s_summary
    )

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
