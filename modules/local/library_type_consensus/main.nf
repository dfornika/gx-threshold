process LIBRARY_TYPE_CONSENSUS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(calls), val(methods)   // calls/methods: parallel lists, one entry per library-type method that ran for this sample

    output:
    tuple val(meta), path("*.library_type_consensus.tsv"), emit: result

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def call_args = [methods, calls].transpose().collect { method, call -> "--call ${method}:${call}" }.join(' \\\n        ')
    """
    library_type_consensus.py \\
        --sample ${meta.id} \\
        --platform ${meta.platform} \\
        ${call_args} > ${prefix}.library_type_consensus.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\t%s\\tstub\\t0\\t0\\t0\\tstub\\t-\\n' "${meta.id}" "${meta.platform}" > ${prefix}.library_type_consensus.tsv
    """
}
