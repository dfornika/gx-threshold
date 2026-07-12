process SUBSAMPLE_READS_HEAD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads)
    val   max_reads   // take the first N records (not a random subsample - a coarse yes/no check doesn't need one, and this avoids a new dependency)

    output:
    tuple val(meta), path("*.sub.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def n_lines = max_reads.toInteger() * 4
    // `head` exits as soon as it has n_lines, closing the pipe - zcat/gzip
    // upstream of it then get SIGPIPE (exit 141) once there's more input
    // than requested, which the pipeline's global pipefail turns into a
    // task failure. That's the intended early stop, not a real error - see
    // modules/local/library_type_aligned/main.nf for the same fix applied
    // to a similar case.
    if (meta.single_end) {
        """
        set +o pipefail
        zcat -f ${reads} | head -n ${n_lines} | gzip > ${prefix}.sub.fastq.gz
        set -o pipefail
        """
    } else {
        """
        set +o pipefail
        zcat -f ${reads[0]} | head -n ${n_lines} | gzip > ${prefix}_R1.sub.fastq.gz
        zcat -f ${reads[1]} | head -n ${n_lines} | gzip > ${prefix}_R2.sub.fastq.gz
        set -o pipefail
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        echo '' | gzip > ${prefix}.sub.fastq.gz
        """
    } else {
        """
        echo '' | gzip > ${prefix}_R1.sub.fastq.gz
        echo '' | gzip > ${prefix}_R2.sub.fastq.gz
        """
    }
}
