process LIBRARY_TYPE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(qc_json)   // fastp (illumina) or fastplong (nanopore) JSON report

    output:
    tuple val(meta), path("*.library_type.tsv"), emit: result

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    classify_library_type.py \\
        --sample ${meta.id} \\
        --platform ${meta.platform} \\
        ${qc_json} > ${prefix}.library_type.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\t%s\\tstub\\tNA\\tNA\\tNA\\t-\\n' "${meta.id}" "${meta.platform}" > ${prefix}.library_type.tsv
    """
}
