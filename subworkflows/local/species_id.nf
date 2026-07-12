//
// Species identification: mash/sourmash/sylph run in parallel (each
// independently toggleable) so their calls can be compared side by side.
// SPECIES_ID_SUMMARY normalises each tool's own output format into one
// comparable row per sample; see tests/data/species_db/README.md.
//

include { MASH_DIST                                       } from '../../modules/nf-core/mash/dist/main'
include { SOURMASH_SKETCH                                 } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_GATHER                                 } from '../../modules/nf-core/sourmash/gather/main'
include { SYLPH_PROFILE                                    } from '../../modules/nf-core/sylph/profile/main'
include { SPECIES_ID_SUMMARY as SPECIES_ID_SUMMARY_MASH     } from '../../modules/local/species_id_summary/main'
include { SPECIES_ID_SUMMARY as SPECIES_ID_SUMMARY_SOURMASH } from '../../modules/local/species_id_summary/main'
include { SPECIES_ID_SUMMARY as SPECIES_ID_SUMMARY_SYLPH    } from '../../modules/local/species_id_summary/main'

workflow SPECIES_ID {

    take:
    ch_clean_reads          // tuple(meta, reads)
    ch_species_id_manifest  // path (or []) - accession -> organism/taxonomy CSV, shared with SPECIES_COMPOSITION_ANALYSIS
    outdir

    main:

    def ch_multiqc_files           = channel.empty()
    def ch_species_id_rows         = channel.empty()
    def ch_sourmash_gather_result  = channel.empty()

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
    // Its gather output is also reused directly by SPECIES_COMPOSITION_ANALYSIS.
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
        ch_sourmash_gather_result = SOURMASH_GATHER.out.result
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
    // .ifEmpty([]): see reference_genome.nf - collectFile emits nothing at
    // all if every species-ID tool is off, so SAMPLE_SUMMARY still gets a
    // usable value.
    def ch_species_id_summary = ch_species_id_rows
        .map { _meta, file -> file }
        .collectFile(
            name: 'species_id_summary.tsv',
            storeDir: "${outdir}/species_id",
            sort: true,
            seed: "sample\tplatform\ttool\taccession\torganism\tspecies_taxid\tspecies_name\tmetric\tvalue\n"
        )
        .ifEmpty([])

    emit:
    species_id_rows        = ch_species_id_rows        // tuple(meta, species_id.tsv) - one per tool per sample, for REFERENCE_GENOME
    sourmash_gather_result  = ch_sourmash_gather_result // tuple(meta, gather csv.gz) - empty channel if sourmash is off, for SPECIES_COMPOSITION_ANALYSIS
    multiqc_files           = ch_multiqc_files
    species_id_summary      = ch_species_id_summary     // path - for SAMPLE_SUMMARY
}
