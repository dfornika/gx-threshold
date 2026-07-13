process ASSEMBLY_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(quast), path(checkm2)   // checkm2 optional (staged as '' when CheckM2 did not run)

    output:
    tuple val(meta), path("*.assembly_summary.tsv"), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def checkm2_arg = checkm2 ? "--checkm2 ${checkm2}" : ''
    """
    parse_assembly_qc.py \\
        --sample ${meta.id} \\
        --platform ${meta.platform} \\
        --assembler ${meta.assembler} \\
        --quast ${quast} \\
        ${checkm2_arg} > ${prefix}.assembly_summary.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\t%s\\t%s\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\n' "${meta.id}" "${meta.platform}" "${meta.assembler}" > ${prefix}.assembly_summary.tsv
    """
}
