process LIBRARY_TYPE_PILEUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(sam)   // sam: from ALIGN_READS

    output:
    tuple val(meta), path("*.library_type_pileup.tsv"), emit: result

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single_end_flag = meta.single_end ? '--single-end' : ''
    // Nanopore primer-trim/basecalling slop needs a looser start tolerance
    // than precisely-trimmed Illumina reads, and its 3' end is noisier still
    // (no mate to cross-check against), hence the much looser end tolerance.
    def start_tolerance = meta.platform == 'nanopore' ? 15  : 5
    def end_tolerance   = meta.platform == 'nanopore' ? 100 : 5
    """
    zcat ${sam} > aln.sam
    classify_library_type_pileup.py \\
        --sample ${meta.id} \\
        --platform ${meta.platform} \\
        ${single_end_flag} \\
        --start-tolerance ${start_tolerance} \\
        --end-tolerance ${end_tolerance} \\
        aln.sam > ${prefix}.library_type_pileup.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\t%s\\tstub\\t0\\tNA\\tNA\\tNA\\t-\\n' "${meta.id}" "${meta.platform}" > ${prefix}.library_type_pileup.tsv
    """
}
