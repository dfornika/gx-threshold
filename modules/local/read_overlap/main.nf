process READ_OVERLAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.overlap.paf.gz"), path("*.n_reads.txt"), emit: overlap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // All-vs-all self overlap - no reference genome needed. Only one file
    // per fragment (R1 for paired-end) goes in: mates of the same fragment
    // don't overlap each other's sequence, so including R2 would only add
    // noise/cost, not signal.
    def r1     = meta.single_end ? reads : reads[0]
    // Long-read (ava-ont) vs short, highly-accurate reads: minimap2's own
    // preset tuning handles the different error/indel profiles, so per-read
    // identity thresholds don't need to be hand-tuned here - only the
    // downstream clustering thresholds in LIBRARY_TYPE_CLUSTER do.
    def preset = meta.platform == 'nanopore' ? 'ava-ont' : 'sr'
    """
    minimap2 -x ${preset} -t ${task.cpus} ${args} ${r1} ${r1} \\
        | gzip > ${prefix}.overlap.paf.gz

    zcat ${r1} | wc -l | awk '{print int(\$1/4)}' > ${prefix}.n_reads.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.overlap.paf.gz
    echo '0' > ${prefix}.n_reads.txt
    """
}
