process ALIGN_READS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads), path(reference)   // reference: FASTA, optionally gzipped

    output:
    tuple val(meta), path("*.sam.gz"), emit: sam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Same platform -> preset mapping DEHOST/LIBRARY_TYPE_ALIGNED already use.
    def preset = meta.platform == 'nanopore' ? 'map-ont' : 'sr'
    """
    case "${reference}" in
        *.gz) zcat "${reference}" > ref.fasta ;;
        *)    ln -s "${reference}" ref.fasta ;;
    esac

    minimap2 -ax ${preset} -t ${task.cpus} ${args} ref.fasta ${reads} \\
        | gzip > ${prefix}.sam.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.sam.gz
    """
}
