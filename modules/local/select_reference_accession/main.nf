process SELECT_REFERENCE_ACCESSION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(species_id_rows)   // one or more *.species_id.tsv files for this sample (from SPECIES_ID_SUMMARY)

    output:
    tuple val(meta), path("*.reference_selection.tsv"), emit: selection

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_reference_accession.py \\
        --sample ${meta.id} \\
        --platform ${meta.platform} \\
        ${species_id_rows} > ${prefix}.reference_selection.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\t%s\\tNA\\tNA\\tNA\\tstub\\n' "${meta.id}" "${meta.platform}" > ${prefix}.reference_selection.tsv
    """
}
