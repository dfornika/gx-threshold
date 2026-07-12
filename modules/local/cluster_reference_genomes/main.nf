process CLUSTER_REFERENCE_GENOMES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    path distances
    path manifest, stageAs: 'manifest.csv?'   // optional: accession,organism,... - omit to report bare accessions in collapsed-group notes

    output:
    path("reference_genome_clusters.csv"), emit: clusters

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def manifest_arg = manifest ? "--manifest ${manifest}" : ''
    """
    cluster_reference_genomes.py \\
        ${args} \\
        ${manifest_arg} \\
        ${distances} reference_genome_clusters.csv
    """

    stub:
    """
    printf 'accession,cluster_id,organism\\n' > reference_genome_clusters.csv
    """
}
