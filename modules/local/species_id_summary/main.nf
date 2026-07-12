process SPECIES_ID_SUMMARY {
    tag "$meta.id:$tool"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(result), val(tool)
    path  manifest, stageAs: 'manifest.csv?'   // optional: accession,organism,strain,category - omit to report raw accessions

    output:
    tuple val(meta), path("*.species_id.tsv"), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${tool}"
    def manifest_arg = manifest ? "--manifest ${manifest}" : ''
    """
    parse_species_id.py \\
        --tool ${tool} \\
        --sample ${meta.id} \\
        --platform ${meta.platform} \\
        ${manifest_arg} \\
        ${result} > ${prefix}.species_id.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${tool}"
    """
    printf '%s\\t%s\\t%s\\tNA\\tstub\\tNA\\tNA\\n' "${meta.id}" "${meta.platform}" "${tool}" > ${prefix}.species_id.tsv
    """
}
