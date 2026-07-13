//
// Reference genome fetch/cache, for the alignment-based library-type
// approaches (e.g. amplicon-vs-shotgun detection for Nanopore, where the
// fastp-JSON-based heuristic doesn't apply - see docs/testing.md). Consumes
// the species-ID consensus accession picked by SPECIES_ID
// (SELECT_REFERENCE_ACCESSION - majority vote on species_taxid, falling back
// to sylph > mash > sourmash - see bin/select_reference_accession.py), then
// fetches and caches the genome, keyed by accession so it's never
// re-downloaded across runs (Nextflow storeDir). Toggle with
// --skip_reference_genome_fetch / --reference_genome_cache_dir. On by
// default (network access assumed in deployment) - a fetch failure for one
// accession degrades gracefully rather than failing the run, see
// conf/modules.config's FETCH_REFERENCE_GENOME errorStrategy.
//

include { FETCH_REFERENCE_GENOME } from '../../modules/local/fetch_reference_genome/main'

workflow REFERENCE_GENOME {

    take:
    ch_species_id_consensus   // tuple(meta, reference_selection.tsv) - from SPECIES_ID (SELECT_REFERENCE_ACCESSION), one row per sample
    outdir

    main:

    def ch_multiqc_files    = channel.empty()
    def ch_sample_reference = channel.empty()

    if (!params.skip_reference_genome_fetch) {
        if (!params.reference_genome_cache_dir) {
            error("Reference genome fetch is enabled but --reference_genome_cache_dir was not set. Provide a cache directory or run with --skip_reference_genome_fetch.")
        }
        def ch_selected_accession = ch_species_id_consensus
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
        // that had a species-ID consensus and whose fetch succeeded (a
        // fetch that's retried and then ignored after failure just leaves
        // that sample out of this channel, same as "no consensus").
    }

    emit:
    sample_reference = ch_sample_reference   // tuple(meta, cached reference FASTA) - empty channel if this stage is off
    multiqc_files     = ch_multiqc_files
}
