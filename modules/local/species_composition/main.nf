process SPECIES_COMPOSITION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(gather_csv)   // from SOURMASH_GATHER
    path clusters                       // from CLUSTER_REFERENCE_GENOMES
    path manifest, stageAs: 'manifest.csv?'   // optional: accession,organism,... - omit to report bare accessions

    output:
    tuple val(meta), path("*.species_composition.tsv"), emit: result

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def manifest_arg = manifest ? "--manifest ${manifest}" : ''
    """
    classify_species_composition.py \\
        --sample ${meta.id} \\
        --platform ${meta.platform} \\
        --clusters ${clusters} \\
        ${manifest_arg} \\
        ${gather_csv} > ${prefix}.species_composition.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\t%s\\tstub\\t0\\tNA\\tNA\\t0\\tNA\\tNA\\tNA\\t-\\n' "${meta.id}" "${meta.platform}" > ${prefix}.species_composition.tsv
    """
}
