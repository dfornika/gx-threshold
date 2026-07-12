process LIBRARY_TYPE_ALIGNED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads), path(reference)   // reference: FASTA, optionally gzipped

    output:
    tuple val(meta), path("*.library_type_aligned.tsv"), emit: result

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''          // extra minimap2 args
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Same platform -> preset mapping DEHOST already uses.
    def preset = meta.platform == 'nanopore' ? 'map-ont' : 'sr'
    """
    case "${reference}" in
        *.gz) zcat "${reference}" > ref.fasta ;;
        *)    ln -s "${reference}" ref.fasta ;;
    esac
    samtools faidx ref.fasta
    genome_length=\$(awk '{sum+=\$2} END{print sum}' ref.fasta.fai)

    # A confident early verdict makes classify_library_type_aligned.awk exit
    # before consuming all of minimap2's output, which sends minimap2 SIGPIPE
    # (exit 141) - that's the intended early stop, not a failure, so pipefail
    # is scoped off around just this one pipe (the pipeline-wide default set
    # in nextflow.config stays on for everything else in this script).
    set +o pipefail
    minimap2 -ax ${preset} -t ${task.cpus} ${args} ref.fasta ${reads} \\
        | classify_library_type_aligned.awk \\
            -v sample=${meta.id} \\
            -v platform=${meta.platform} \\
            -v genome_length="\$genome_length" \\
            > ${prefix}.library_type_aligned.tsv
    set -o pipefail
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\t%s\\tstub\\t0\\t0\\tstub\\n' "${meta.id}" "${meta.platform}" > ${prefix}.library_type_aligned.tsv
    """
}
