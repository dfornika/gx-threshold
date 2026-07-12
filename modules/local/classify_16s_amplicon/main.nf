process CLASSIFY_16S_AMPLICON {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(sam)   // from ALIGN_READS, aligned against the 16S reference database

    output:
    tuple val(meta), path("*.16s_amplicon.tsv"), emit: result

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Nanopore's higher error rate needs a more permissive identity
    // threshold than short, highly-accurate Illumina reads - same reasoning
    // (and similar values) as LIBRARY_TYPE_PILEUP's platform split.
    def min_identity = meta.platform == 'nanopore' ? 0.80 : 0.90
    def min_coverage = meta.platform == 'nanopore' ? 0.70 : 0.80
    """
    zcat ${sam} > aln.sam
    classify_16s_amplicon.py \\
        --sample ${meta.id} \\
        --platform ${meta.platform} \\
        --min-identity ${min_identity} \\
        --min-coverage ${min_coverage} \\
        aln.sam > ${prefix}.16s_amplicon.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\t%s\\tstub\\t0\\t0\\tNA\\t-\\n' "${meta.id}" "${meta.platform}" > ${prefix}.16s_amplicon.tsv
    """
}
