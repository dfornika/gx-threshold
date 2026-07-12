process LIBRARY_TYPE_CLUSTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(paf), path(n_reads_file)   // paf + n_reads.txt: from READ_OVERLAP

    output:
    tuple val(meta), path("*.library_type_cluster.tsv"), emit: result

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Nanopore's higher error/indel rate needs a more permissive identity
    // threshold than short, highly-accurate Illumina reads; both still
    // require most of the shorter read's length to be covered, so a small
    // shared end (as assembly-style overlap detection would allow) doesn't
    // count as "these reads are from the same PCR product".
    def identity_threshold = meta.platform == 'nanopore' ? 0.85 : 0.95
    def coverage_threshold = meta.platform == 'nanopore' ? 0.7  : 0.8
    """
    zcat ${paf} > overlap.paf
    n_reads=\$(cat ${n_reads_file})
    classify_library_type_cluster.py \\
        --sample ${meta.id} \\
        --platform ${meta.platform} \\
        --total-reads \$n_reads \\
        --identity-threshold ${identity_threshold} \\
        --coverage-threshold ${coverage_threshold} \\
        overlap.paf > ${prefix}.library_type_cluster.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\t%s\\tstub\\t0\\tNA\\tNA\\tNA\\t-\\n' "${meta.id}" "${meta.platform}" > ${prefix}.library_type_cluster.tsv
    """
}
