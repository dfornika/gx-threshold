process REFERENCE_GENOME_DISTANCES {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

    input:
    path mash_db

    output:
    path("reference_genome_distances.tsv"), emit: distances

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mash dist ${mash_db} ${mash_db} > reference_genome_distances.tsv
    """

    stub:
    """
    touch reference_genome_distances.tsv
    """
}
