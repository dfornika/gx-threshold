//
// Reference genome selection + fetch/cache, for the alignment-based
// library-type approaches (e.g. amplicon-vs-shotgun detection for Nanopore,
// where the fastp-JSON-based heuristic doesn't apply - see docs/testing.md).
// Picks one accession per sample from the species-ID tools' calls (majority
// vote on species_taxid, falling back to sylph > mash > sourmash - see
// bin/select_reference_accession.py), then fetches and caches the genome,
// keyed by accession so it's never re-downloaded across runs (Nextflow
// storeDir). Off by default - unlike the species-ID stages, this does real
// network I/O rather than using a bundled test database, so it isn't
// exercised by the standard test profiles. Toggle with
// --skip_reference_genome_fetch / --reference_genome_cache_dir.
//

include { SELECT_REFERENCE_ACCESSION } from '../../modules/local/select_reference_accession/main'
include { FETCH_REFERENCE_GENOME     } from '../../modules/local/fetch_reference_genome/main'

workflow REFERENCE_GENOME {

    take:
    ch_species_id_rows   // tuple(meta, species_id.tsv) - from SPECIES_ID, one row per tool per sample
    outdir

    main:

    def ch_multiqc_files    = channel.empty()
    def ch_sample_reference = channel.empty()

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
        ch_sample_reference = ch_selected_accession.found
            .map { meta, accession -> tuple(accession, meta) }
            .combine(ch_reference_by_accession, by: 0)
            .map { _accession, meta, fasta -> tuple(meta, fasta) }
        // ch_sample_reference: tuple(meta, cached reference FASTA) per sample
        // that had a species-ID consensus.

        SELECT_REFERENCE_ACCESSION.out.selection
            .map { _meta, file -> file }
            .collectFile(
                name: 'reference_selection_summary.tsv',
                storeDir: "${outdir}/reference_genome",
                sort: true,
                seed: "sample\tplatform\taccession\tspecies_taxid\tspecies_name\tmethod\n"
            )
    }

    emit:
    sample_reference = ch_sample_reference   // tuple(meta, cached reference FASTA) - empty channel if this stage is off
    multiqc_files    = ch_multiqc_files
}
