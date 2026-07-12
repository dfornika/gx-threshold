//
// Pure-culture vs. metagenomic detection from species-ID composition
// breadth - sourmash gather already decomposes a sample into the set of
// reference genomes that best explain it, so "does one genome explain
// nearly everything, or does it take many" is read directly off its
// output. Needs mash too: its reference genomes are clustered by ANI
// first (REFERENCE_GENOME_DISTANCES + CLUSTER_REFERENCE_GENOMES), so
// nomenclature artifacts like Escherichia coli/Shigella spp. (~98%+ ANI
// to each other - a pre-genomic-era clinical naming split, not a real
// genomic distinction, the same issue widely seen with Kraken2) don't
// look like spurious extra "species" in a pure culture. Both the naive
// (pre-collapse) and adjusted (post-collapse) breadth are reported side
// by side - see docs/testing.md. Toggle with --skip_species_composition.
//

include { REFERENCE_GENOME_DISTANCES } from '../../modules/local/reference_genome_distances/main'
include { CLUSTER_REFERENCE_GENOMES  } from '../../modules/local/cluster_reference_genomes/main'
include { SPECIES_COMPOSITION        } from '../../modules/local/species_composition/main'
include { tagMetaFromVerdict         } from './meta_utils'

workflow SPECIES_COMPOSITION_ANALYSIS {

    take:
    ch_clean_reads             // tuple(meta, reads) - to be tagged with meta.composition
    ch_sourmash_gather_result  // tuple(meta, gather csv.gz) - from SPECIES_ID
    ch_species_id_manifest     // path (or []) - shared with SPECIES_ID
    outdir

    main:

    def ch_multiqc_files = channel.empty()
    def ch_tagged_reads  = ch_clean_reads
    def ch_species_composition_summary = channel.value([])

    if (!params.skip_species_composition) {
        if (params.skip_mash || params.skip_sourmash) {
            error("Species composition detection needs both mash and sourmash enabled (it clusters mash's reference genomes by ANI, then applies that to sourmash's gather output). Run with --skip_species_composition to turn this stage off, or enable both mash and sourmash.")
        }
        def ch_mash_db = channel.value(file(params.mash_db, checkIfExists: true))
        REFERENCE_GENOME_DISTANCES(ch_mash_db)
        CLUSTER_REFERENCE_GENOMES(REFERENCE_GENOME_DISTANCES.out.distances, ch_species_id_manifest)
        SPECIES_COMPOSITION(ch_sourmash_gather_result, CLUSTER_REFERENCE_GENOMES.out.clusters, ch_species_id_manifest)
        ch_multiqc_files = ch_multiqc_files.mix(SPECIES_COMPOSITION.out.result.map { _meta, file -> file })

        // .ifEmpty([]): see subworkflows/local/reference_genome.nf for why -
        // collectFile emits nothing at all if e.g. sourmash gather found no
        // hit for any sample this run (a real, legitimate case - the plain
        // `test` profile's synthetic fixtures don't match anything in the
        // real species DB), so SAMPLE_SUMMARY still gets a usable value.
        ch_species_composition_summary = SPECIES_COMPOSITION.out.result
            .map { _meta, file -> file }
            .collectFile(
                name: 'species_composition_summary.tsv',
                storeDir: "${outdir}/species_id",
                sort: true,
                seed: "sample\tplatform\tverdict\tnaive_n_hits\tnaive_top_hit_frac\tnaive_effective_n\tadjusted_n_hits\tadjusted_top_hit_frac\tadjusted_effective_n\tadjusted_fraction_explained\tnote\n"
            )
            .ifEmpty([])

        ch_tagged_reads = tagMetaFromVerdict(ch_tagged_reads, SPECIES_COMPOSITION.out.result, 'composition')
    }

    emit:
    reads                       = ch_tagged_reads   // tuple(meta, reads) - meta.composition set if this stage ran
    multiqc_files               = ch_multiqc_files
    species_composition_summary = ch_species_composition_summary // path (or []) - for SAMPLE_SUMMARY
}
