process DEHOST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads)
    path  reference       // host reference: FASTA (optionally .gz) or a prebuilt minimap2 .mmi index
    val   scrub_headers   // Boolean: replace read headers with anonymised indices (off by default)

    output:
    tuple val(meta), path("*.dehosted*.fastq.gz"), emit: reads
    tuple val(meta), path("*.dehost.stats.tsv")  , emit: stats
    tuple val("${task.process}"), val("minimap2"), eval("minimap2 --version"),                        emit: versions_minimap2, topic: versions
    tuple val("${task.process}"), val("samtools"), eval("samtools --version | head -n1 | sed 's/samtools //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''          // extra minimap2 args
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Platform-specific minimap2 preset. Everything downstream keys off meta.platform,
    // which is set once in the main workflow (nanopore vs illumina).
    def preset = meta.platform == 'nanopore' ? 'map-ont' : 'sr'

    if (meta.single_end) {
        // Long-read / single-end: keep primary unmapped reads (-f 4), drop secondary/supplementary (-F 0x900).
        def scrub = scrub_headers
            ? "| awk 'NR%4==1{print \"@${prefix}_\"(++c)} NR%4==2{print} NR%4==3{print \"+\"} NR%4==0{print}'"
            : ""
        """
        minimap2 -ax ${preset} -t ${task.cpus} ${args} ${reference} ${reads} \\
            | samtools view -@ ${task.cpus} -b -o aln.bam -

        primary=\$(samtools view -c -F 0x900 aln.bam)
        retained=\$(samtools view -c -f 4 -F 0x900 aln.bam)
        removed=\$((primary - retained))

        samtools fastq -@ ${task.cpus} -f 4 -F 0x900 -n aln.bam ${scrub} | gzip > ${prefix}.dehosted.fastq.gz

        write_stats() {
            printf 'sample\\tplatform\\tinput_reads\\thost_reads\\tdehosted_reads\\tpercent_host\\n' > ${prefix}.dehost.stats.tsv
            pct=\$(awk -v r=\$removed -v p=\$primary 'BEGIN{ if (p>0) printf "%.4f", (r/p)*100; else print "0.0000" }')
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "${meta.platform}" "\$primary" "\$removed" "\$retained" "\$pct" >> ${prefix}.dehost.stats.tsv
        }
        write_stats
        """
    } else {
        // Short-read / paired-end: keep a pair only if BOTH mates are unmapped (-f 12),
        // i.e. drop the whole pair if either mate hits the host. Collate so mates are adjacent for samtools fastq.
        def scrub_r1 = scrub_headers ? "&& zcat ${prefix}.dehosted_R1.raw.fastq.gz | awk 'NR%4==1{print \"@${prefix}_\"(++c)} NR%4==2{print} NR%4==3{print \"+\"} NR%4==0{print}' | gzip > ${prefix}.dehosted_R1.fastq.gz" : "&& mv ${prefix}.dehosted_R1.raw.fastq.gz ${prefix}.dehosted_R1.fastq.gz"
        def scrub_r2 = scrub_headers ? "&& zcat ${prefix}.dehosted_R2.raw.fastq.gz | awk 'NR%4==1{print \"@${prefix}_\"(++c)} NR%4==2{print} NR%4==3{print \"+\"} NR%4==0{print}' | gzip > ${prefix}.dehosted_R2.fastq.gz" : "&& mv ${prefix}.dehosted_R2.raw.fastq.gz ${prefix}.dehosted_R2.fastq.gz"
        """
        minimap2 -ax ${preset} -t ${task.cpus} ${args} ${reference} ${reads} \\
            | samtools view -@ ${task.cpus} -b -o aln.bam -

        primary=\$(samtools view -c -F 0x900 aln.bam)
        retained=\$(samtools view -c -f 12 -F 0x900 aln.bam)
        removed=\$((primary - retained))

        samtools collate -@ ${task.cpus} -u -O aln.bam \\
            | samtools fastq -@ ${task.cpus} -f 12 -F 0x900 -n \\
                -1 ${prefix}.dehosted_R1.raw.fastq.gz \\
                -2 ${prefix}.dehosted_R2.raw.fastq.gz \\
                -0 /dev/null -s /dev/null -
        true ${scrub_r1}
        true ${scrub_r2}
        rm -f ${prefix}.dehosted_R1.raw.fastq.gz ${prefix}.dehosted_R2.raw.fastq.gz

        printf 'sample\\tplatform\\tinput_reads\\thost_reads\\tdehosted_reads\\tpercent_host\\n' > ${prefix}.dehost.stats.tsv
        pct=\$(awk -v r=\$removed -v p=\$primary 'BEGIN{ if (p>0) printf "%.4f", (r/p)*100; else print "0.0000" }')
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "${meta.platform}" "\$primary" "\$removed" "\$retained" "\$pct" >> ${prefix}.dehost.stats.tsv
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        echo '' | gzip > ${prefix}.dehosted.fastq.gz
        printf 'sample\\tplatform\\tinput_reads\\thost_reads\\tdehosted_reads\\tpercent_host\\n' > ${prefix}.dehost.stats.tsv
        printf '%s\\t%s\\t0\\t0\\t0\\t0.0000\\n' "${meta.id}" "${meta.platform}" >> ${prefix}.dehost.stats.tsv
        """
    } else {
        """
        echo '' | gzip > ${prefix}.dehosted_R1.fastq.gz
        echo '' | gzip > ${prefix}.dehosted_R2.fastq.gz
        printf 'sample\\tplatform\\tinput_reads\\thost_reads\\tdehosted_reads\\tpercent_host\\n' > ${prefix}.dehost.stats.tsv
        printf '%s\\t%s\\t0\\t0\\t0\\t0.0000\\n' "${meta.id}" "${meta.platform}" >> ${prefix}.dehost.stats.tsv
        """
    }
}
