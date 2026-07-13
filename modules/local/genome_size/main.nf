process GENOME_SIZE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(reference), path(consensus)   // reference/consensus optional (staged empty when absent)
    path  table
    val   default_genome_size

    output:
    tuple val(meta), path("*.genome_size.txt"), emit: size

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_arg = reference ? "--reference ${reference}" : ''
    def consensus_arg = consensus ? "--consensus ${consensus}" : ''
    def table_arg     = table     ? "--table ${table}"         : ''
    """
    resolve_genome_size.py \\
        ${reference_arg} \\
        ${consensus_arg} \\
        ${table_arg} \\
        --default ${default_genome_size} > ${prefix}.genome_size.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\tstub\\n' "${default_genome_size}" > ${prefix}.genome_size.txt
    """
}
