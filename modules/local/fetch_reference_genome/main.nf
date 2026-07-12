process FETCH_REFERENCE_GENOME {
    tag "$accession"
    label 'process_low'
    // storeDir (not publishDir - deliberately no entry in conf/modules.config
    // for this process) caches into --reference_genome_cache_dir directly, a
    // persistent, shared location outside any single run's outdir.
    storeDir "${params.reference_genome_cache_dir}/${accession}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1' :
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    val accession

    output:
    path("${accession}.fasta.gz"), emit: fasta
    // No versions-topic emit here: storeDir only supports plain `val`/`path`
    // outputs, not `tuple` (which the versions-topic convention used
    // elsewhere in this pipeline requires) - confirmed by Nextflow refusing
    // to run this process otherwise ("storeDir can only be used with `val`
    // and `path` outputs"). entrez-direct's version is pinned in
    // environment.yml/the container tag instead.

    when:
    task.ext.when == null || task.ext.when

    script:
    // Whole-assembly FASTA (all contigs/plasmids) via NCBI's E-utilities:
    // resolve the assembly accession to its nucleotide sequences
    // (assembly -> nuccore, INSDC set) and fetch them all in one go.
    // NCBI's own `datasets` CLI covers this in one call, but no current,
    // working container build exists for it (tested and confirmed against
    // the live API); `ncbi-genome-download` has a live, maintainer-closed
    // parsing bug against NCBI's current assembly_summary.txt (upstream
    // added an unescaped-tab field). entrez-direct is the option that
    // actually works, using the same stable E-utilities protocol NCBI has
    // kept backwards-compatible for decades.
    """
    esearch -db assembly -query "${accession}[Assembly Accession]" \\
        | elink -target nuccore -name assembly_nuccore_insdc \\
        | efetch -format fasta \\
        | gzip > ${accession}.fasta.gz

    if [ ! -s ${accession}.fasta.gz ] || [ \$(zcat ${accession}.fasta.gz | head -c1 | wc -c) -eq 0 ]; then
        echo "No sequences retrieved for assembly accession ${accession}" >&2
        exit 1
    fi
    """

    stub:
    """
    echo '' | gzip > ${accession}.fasta.gz
    """
}
