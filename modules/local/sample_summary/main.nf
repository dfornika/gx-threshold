process SAMPLE_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    path samples                     // sample manifest (sample, platform) - always present, anchors the row list
    path dehost, stageAs: 'dehost_summary.tsv?'
    path library_type, stageAs: 'library_type.tsv?'
    path library_type_cluster, stageAs: 'library_type_cluster.tsv?'
    path library_type_aligned, stageAs: 'library_type_aligned_summary.tsv?'
    path library_type_pileup, stageAs: 'library_type_pileup_summary.tsv?'
    path library_type_consensus, stageAs: 'library_type_consensus_summary.tsv?'
    path species_id, stageAs: 'species_id_summary.tsv?'
    path species_composition, stageAs: 'species_composition_summary.tsv?'
    path species_id_consensus, stageAs: 'species_id_consensus_summary.tsv?'
    path sixteen_s, stageAs: 'sixteen_s_detection_summary.tsv?'
    path assembly, stageAs: 'assembly_summary.tsv?'

    output:
    path "sample_summary.csv", emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def dehost_arg                    = dehost                    ? "--dehost ${dehost}"                                           : ''
    def library_type_arg              = library_type              ? "--library-type ${library_type}"                             : ''
    def library_type_cluster_arg      = library_type_cluster      ? "--library-type-cluster ${library_type_cluster}"              : ''
    def library_type_aligned_arg      = library_type_aligned      ? "--library-type-aligned ${library_type_aligned}"              : ''
    def library_type_pileup_arg       = library_type_pileup       ? "--library-type-pileup ${library_type_pileup}"                : ''
    def library_type_consensus_arg    = library_type_consensus    ? "--library-type-consensus ${library_type_consensus}"          : ''
    def species_id_arg                = species_id                ? "--species-id ${species_id}"                                  : ''
    def species_composition_arg       = species_composition       ? "--species-composition ${species_composition}"                : ''
    def species_id_consensus_arg      = species_id_consensus      ? "--species-id-consensus ${species_id_consensus}"              : ''
    def sixteen_s_arg                 = sixteen_s                 ? "--sixteen-s ${sixteen_s}"                                    : ''
    def assembly_arg                  = assembly                  ? "--assembly ${assembly}"                                      : ''
    """
    build_sample_summary.py \\
        --samples ${samples} \\
        ${dehost_arg} \\
        ${library_type_arg} \\
        ${library_type_cluster_arg} \\
        ${library_type_aligned_arg} \\
        ${library_type_pileup_arg} \\
        ${library_type_consensus_arg} \\
        ${species_id_arg} \\
        ${species_composition_arg} \\
        ${species_id_consensus_arg} \\
        ${sixteen_s_arg} \\
        ${assembly_arg} \\
        > sample_summary.csv
    """

    stub:
    """
    touch sample_summary.csv
    """
}
