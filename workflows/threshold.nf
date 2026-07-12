/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { FASTPLONG              } from '../modules/nf-core/fastplong/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { MASH_DIST              } from '../modules/nf-core/mash/dist/main'
include { SOURMASH_SKETCH        } from '../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_GATHER        } from '../modules/nf-core/sourmash/gather/main'
include { SYLPH_PROFILE          } from '../modules/nf-core/sylph/profile/main'
include { DEHOST                 } from '../modules/local/dehost/main'
include { LIBRARY_TYPE           } from '../modules/local/library_type/main'
include { SPECIES_ID_SUMMARY as SPECIES_ID_SUMMARY_MASH     } from '../modules/local/species_id_summary/main'
include { SPECIES_ID_SUMMARY as SPECIES_ID_SUMMARY_SOURMASH } from '../modules/local/species_id_summary/main'
include { SPECIES_ID_SUMMARY as SPECIES_ID_SUMMARY_SYLPH    } from '../modules/local/species_id_summary/main'
include { SELECT_REFERENCE_ACCESSION } from '../modules/local/select_reference_accession/main'
include { FETCH_REFERENCE_GENOME     } from '../modules/local/fetch_reference_genome/main'
include { LIBRARY_TYPE_ALIGNED       } from '../modules/local/library_type_aligned/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_threshold_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow THRESHOLD {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Tag each sample's meta with its platform. `single_end` samples (no
    // fastq_2) are assumed to be Nanopore; paired-end samples are assumed to
    // be Illumina. This is the only place that assumption is encoded -
    // everything downstream reads `meta.platform` instead.
    //
    def ch_reads = ch_samplesheet.map { meta, reads ->
        tuple(meta + [platform: meta.single_end ? 'nanopore' : 'illumina'], reads)
    }

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // Split reads by platform for platform-specific trimming/QC tools.
    //
    def ch_reads_by_platform = ch_reads.branch { meta, _reads ->
        long_reads:  meta.platform == 'nanopore'
        short_reads: meta.platform == 'illumina'
    }

    //
    // MODULE: Run fastp on Illumina paired-end reads
    //
    FASTP(ch_reads_by_platform.short_reads.map { meta, reads -> tuple(meta, reads, []) }, false, false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map{ _meta, file -> file })

    //
    // MODULE: Run fastplong on Nanopore single-end reads
    //
    FASTPLONG(ch_reads_by_platform.long_reads, [], false, false)
    ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG.out.json.map{ _meta, file -> file })

    //
    // Recombine platform-specific trimmed reads into a single channel.
    //
    def ch_trimmed = FASTP.out.reads.mix(FASTPLONG.out.reads)

    //
    // MODULE: Library type - classify amplicon vs. shotgun from the fastp/
    // fastplong QC JSON (duplication rate + insert-size histogram shape).
    // Illumina paired-end only for now - fastplong (Nanopore) reports neither
    // field, so long-read samples come back "not_classified"; see
    // docs/testing.md for why. Toggle with --skip_library_type.
    //
    if (!params.skip_library_type) {
        def ch_trim_json = FASTP.out.json.mix(FASTPLONG.out.json)
        LIBRARY_TYPE(ch_trim_json)
        ch_multiqc_files = ch_multiqc_files.mix(LIBRARY_TYPE.out.result.map { _meta, file -> file })
        LIBRARY_TYPE.out.result
            .map { _meta, file -> file }
            .collectFile(
                name: 'library_type.tsv',
                storeDir: "${outdir}/library_type",
                sort: true,
                seed: "sample\tplatform\tverdict\tinsert_concentration\tduplication_rate\tinsert_peak\tnote\n"
            )
    }

    //
    // MODULE: Dehost - remove host (human) reads by aligning to a host reference
    // and keeping the unmapped reads. This runs early so that every downstream
    // step (and any shared reads) is host-depleted. Toggle with --skip_dehosting.
    //
    def ch_clean_reads = ch_trimmed
    if (!params.skip_dehosting) {
        if (!params.dehost_reference) {
            error("Dehosting is enabled but --dehost_reference was not set. Provide a host reference (FASTA or minimap2 .mmi) or run with --skip_dehosting.")
        }
        def ch_host_reference = channel.value(file(params.dehost_reference, checkIfExists: true))
        DEHOST(ch_trimmed, ch_host_reference, params.dehost_scrub_headers)
        ch_clean_reads = DEHOST.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(DEHOST.out.stats.map { _meta, file -> file })
    }

    //
    // Species identification: mash/sourmash/sylph run in parallel (each
    // independently toggleable) so their calls can be compared side by side.
    // SPECIES_ID_SUMMARY normalises each tool's own output format into one
    // comparable row per sample; see tests/data/species_db/README.md.
    //
    def ch_species_id_rows = channel.empty()
    def ch_species_id_manifest = params.species_id_manifest ? file(params.species_id_manifest, checkIfExists: true) : []

    //
    // MODULE: Mash - species identification via k-mer distance to a panel of
    // reference genome sketches. One of three species-ID tools currently under
    // evaluation (mash/sourmash/sylph) - toggle with --skip_mash.
    //
    if (!params.skip_mash) {
        if (!params.mash_db) {
            error("Mash species ID is enabled but --mash_db was not set. Provide a Mash sketch (.msh) or run with --skip_mash.")
        }
        def ch_mash_db = channel.value(file(params.mash_db, checkIfExists: true))
        MASH_DIST(ch_clean_reads, ch_mash_db)
        ch_multiqc_files = ch_multiqc_files.mix(MASH_DIST.out.dist.map { _meta, file -> file })
        SPECIES_ID_SUMMARY_MASH(MASH_DIST.out.dist.map { meta, file -> tuple(meta, file, 'mash') }, ch_species_id_manifest)
        ch_species_id_rows = ch_species_id_rows.mix(SPECIES_ID_SUMMARY_MASH.out.summary)
    }

    //
    // MODULE: Sourmash - species identification via FracMinHash containment
    // (sketch reads, then gather against the reference panel). Second of
    // three species-ID tools under evaluation - toggle with --skip_sourmash.
    //
    if (!params.skip_sourmash) {
        if (!params.sourmash_db) {
            error("Sourmash species ID is enabled but --sourmash_db was not set. Provide a sourmash signature collection (.sig/.sig.zip) or run with --skip_sourmash.")
        }
        def ch_sourmash_db = channel.value(file(params.sourmash_db, checkIfExists: true))
        SOURMASH_SKETCH(ch_clean_reads)
        SOURMASH_GATHER(SOURMASH_SKETCH.out.signatures, ch_sourmash_db, false, false, false, false)
        ch_multiqc_files = ch_multiqc_files.mix(SOURMASH_GATHER.out.result.map { _meta, file -> file })
        SPECIES_ID_SUMMARY_SOURMASH(SOURMASH_GATHER.out.result.map { meta, file -> tuple(meta, file, 'sourmash') }, ch_species_id_manifest)
        ch_species_id_rows = ch_species_id_rows.mix(SPECIES_ID_SUMMARY_SOURMASH.out.summary)
    }

    //
    // MODULE: Sylph - species identification via containment ANI. Sketches
    // and profiles reads against the reference panel in one step (no
    // separate sketch stage needed). Third of three species-ID tools under
    // evaluation - toggle with --skip_sylph.
    //
    if (!params.skip_sylph) {
        if (!params.sylph_db) {
            error("Sylph species ID is enabled but --sylph_db was not set. Provide a Sylph genome database (.syldb) or run with --skip_sylph.")
        }
        def ch_sylph_db = channel.value(file(params.sylph_db, checkIfExists: true))
        SYLPH_PROFILE(ch_clean_reads, ch_sylph_db)
        ch_multiqc_files = ch_multiqc_files.mix(SYLPH_PROFILE.out.profile_out.map { _meta, file -> file })
        SPECIES_ID_SUMMARY_SYLPH(SYLPH_PROFILE.out.profile_out.map { meta, file -> tuple(meta, file, 'sylph') }, ch_species_id_manifest)
        ch_species_id_rows = ch_species_id_rows.mix(SPECIES_ID_SUMMARY_SYLPH.out.summary)
    }

    //
    // Collate the three tools' per-sample calls into one comparison TSV.
    //
    ch_species_id_rows
        .map { _meta, file -> file }
        .collectFile(
            name: 'species_id_summary.tsv',
            storeDir: "${outdir}/species_id",
            sort: true,
            seed: "sample\tplatform\ttool\taccession\torganism\tspecies_taxid\tspecies_name\tmetric\tvalue\n"
        )

    //
    // Reference genome selection + fetch/cache, for future alignment-based
    // work (e.g. amplicon-vs-shotgun detection for Nanopore, where the
    // fastp-JSON-based heuristic doesn't apply - see docs/testing.md).
    // Picks one accession per sample from the species-ID tools' calls
    // (majority vote on species_taxid, falling back to sylph > mash >
    // sourmash - see bin/select_reference_accession.py), then fetches and
    // caches the genome, keyed by accession so it's never re-downloaded
    // across runs (Nextflow storeDir). Off by default - unlike the
    // species-ID stages, this does real network I/O rather than using a
    // bundled test database, so it isn't exercised by the standard test
    // profiles. Toggle with --skip_reference_genome_fetch /
    // --reference_genome_cache_dir.
    //
    if (!params.skip_reference_genome_fetch) {
        if (!params.reference_genome_cache_dir) {
            error("Reference genome fetch is enabled but --reference_genome_cache_dir was not set. Provide a cache directory or run with --skip_reference_genome_fetch.")
        }
        def ch_species_id_by_sample = ch_species_id_rows.groupTuple()
        SELECT_REFERENCE_ACCESSION(ch_species_id_by_sample)
        ch_multiqc_files = ch_multiqc_files.mix(SELECT_REFERENCE_ACCESSION.out.selection.map { _meta, file -> file })

        def ch_selected_accession = SELECT_REFERENCE_ACCESSION.out.selection
            .map { meta, file -> tuple(meta, file.text.trim().split('\t')[2]) }
            .branch { _meta, accession ->
                found: accession != 'NA'
                no_hit: true
            }

        def ch_accessions_to_fetch = ch_selected_accession.found
            .map { _meta, accession -> accession }
            .unique()
        FETCH_REFERENCE_GENOME(ch_accessions_to_fetch)

        def ch_reference_by_accession = FETCH_REFERENCE_GENOME.out.fasta
            .map { fasta -> tuple(fasta.name.replace('.fasta.gz', ''), fasta) }
        def ch_sample_reference = ch_selected_accession.found
            .map { meta, accession -> tuple(accession, meta) }
            .combine(ch_reference_by_accession, by: 0)
            .map { _accession, meta, fasta -> tuple(meta, fasta) }
        // ch_sample_reference: tuple(meta, cached reference FASTA) per sample
        // that had a species-ID consensus.

        //
        // Alignment-based, platform-unified amplicon-vs-shotgun detection
        // (streaming index of dispersion of per-base depth against the
        // fetched reference) - see docs/testing.md. Unlike LIBRARY_TYPE
        // (Illumina-only, fastp-JSON-based), this works on Nanopore too,
        // since it's just consuming this stage's own dependency
        // (ch_sample_reference) rather than a new flag.
        //
        def ch_reads_with_reference = ch_clean_reads
            .combine(ch_sample_reference, by: 0)
        LIBRARY_TYPE_ALIGNED(ch_reads_with_reference)
        ch_multiqc_files = ch_multiqc_files.mix(LIBRARY_TYPE_ALIGNED.out.result.map { _meta, file -> file })

        LIBRARY_TYPE_ALIGNED.out.result
            .map { _meta, file -> file }
            .collectFile(
                name: 'library_type_aligned_summary.tsv',
                storeDir: "${outdir}/library_type",
                sort: true,
                seed: "sample\tplatform\tverdict\tn_reads_used\tindex_of_dispersion\tmethod\n"
            )

        SELECT_REFERENCE_ACCESSION.out.selection
            .map { _meta, file -> file }
            .collectFile(
                name: 'reference_selection_summary.tsv',
                storeDir: "${outdir}/reference_genome",
                sort: true,
                seed: "sample\tplatform\taccession\tspecies_taxid\tspecies_name\tmethod\n"
            )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'gx-threshold_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'gx-threshold'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
