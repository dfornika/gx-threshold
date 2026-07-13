//
// Library type, reference-free: the two amplicon-vs-shotgun approaches that
// need no fetched reference genome, so they can run right after DEHOST with
// no network dependency - LIBRARY_TYPE (fastp/fastplong QC JSON) and
// LIBRARY_TYPE_CLUSTER (read clustering, via READ_OVERLAP). The other two
// approaches under evaluation, LIBRARY_TYPE_ALIGNED and LIBRARY_TYPE_PILEUP,
// need a fetched reference genome (itself dependent on species-ID), so they
// live in ALIGNMENT_BASED_LIBRARY_TYPE instead, run later once one is
// available - see docs/testing.md for how all four approaches compare.
//
// Neither result is folded into `meta.library_type` here: that now happens
// once, after all four methods have had a chance to run, via
// LIBRARY_TYPE_CONSENSUS (subworkflows/local/library_type_consensus.nf) -
// this subworkflow just emits its two raw per-sample results outward for
// that stage to consume.
//

include { LIBRARY_TYPE          } from '../../modules/local/library_type/main'
include { READ_OVERLAP          } from '../../modules/local/read_overlap/main'
include { LIBRARY_TYPE_CLUSTER  } from '../../modules/local/library_type_cluster/main'

workflow LIBRARY_TYPE_REFERENCE_FREE {

    take:
    ch_trim_json    // tuple(meta, fastp/fastplong QC json) - pre-dehost
    ch_clean_reads  // tuple(meta, reads) - post-dehost
    outdir

    main:

    def ch_multiqc_files                = channel.empty()
    def ch_library_type_summary         = channel.value([])
    def ch_library_type_cluster_summary = channel.value([])
    def ch_library_type_result          = channel.empty()
    def ch_library_type_cluster_result  = channel.empty()

    //
    // MODULE: Library type - classify amplicon vs. shotgun from the fastp/
    // fastplong QC JSON (duplication rate + insert-size histogram shape).
    // Illumina paired-end only for now - fastplong (Nanopore) reports neither
    // field, so long-read samples come back "not_classified"; see
    // docs/testing.md for why. Toggle with --skip_library_type.
    //
    if (!params.skip_library_type) {
        LIBRARY_TYPE(ch_trim_json)
        ch_multiqc_files = ch_multiqc_files.mix(LIBRARY_TYPE.out.result.map { _meta, file -> file })
        ch_library_type_result = LIBRARY_TYPE.out.result
        // .ifEmpty([]): see subworkflows/local/species_id.nf for why -
        // collectFile emits nothing at all if its input is completely empty.
        ch_library_type_summary = LIBRARY_TYPE.out.result
            .map { _meta, file -> file }
            .collectFile(
                name: 'library_type.tsv',
                storeDir: "${outdir}/library_type",
                sort: true,
                seed: "sample\tplatform\tverdict\tinsert_concentration\tduplication_rate\tinsert_peak\tnote\n"
            )
            .ifEmpty([])
    }

    //
    // MODULE: Library type, reference-free (amplicon vs. shotgun via read
    // clustering) - a second, independent approach to the same question
    // LIBRARY_TYPE/LIBRARY_TYPE_ALIGNED ask, evaluated in parallel rather
    // than as a replacement (see docs/testing.md for how the approaches
    // compare). Clusters reads via all-vs-all self-overlap (READ_OVERLAP,
    // minimap2) rather than aligning to a fetched reference genome, so -
    // unlike LIBRARY_TYPE_ALIGNED/LIBRARY_TYPE_PILEUP - it needs no
    // species-ID/reference-fetch dependency and works offline on the
    // bundled test fixtures. Toggle with --skip_library_type_cluster.
    //
    if (!params.skip_library_type_cluster) {
        READ_OVERLAP(ch_clean_reads)
        LIBRARY_TYPE_CLUSTER(READ_OVERLAP.out.overlap)
        ch_multiqc_files = ch_multiqc_files.mix(LIBRARY_TYPE_CLUSTER.out.result.map { _meta, file -> file })
        ch_library_type_cluster_result = LIBRARY_TYPE_CLUSTER.out.result

        // .ifEmpty([]): see subworkflows/local/species_id.nf for why.
        ch_library_type_cluster_summary = LIBRARY_TYPE_CLUSTER.out.result
            .map { _meta, file -> file }
            .collectFile(
                name: 'library_type_cluster.tsv',
                storeDir: "${outdir}/library_type",
                sort: true,
                seed: "sample\tplatform\tverdict\tn_reads\tlargest_cluster_frac\tclustered_frac\teffective_clusters_frac\tnote\n"
            )
            .ifEmpty([])
    }

    emit:
    reads                        = ch_clean_reads   // tuple(meta, reads) - unchanged pass-through
    multiqc_files                = ch_multiqc_files
    library_type_summary         = ch_library_type_summary         // path (or []) - for SAMPLE_SUMMARY
    library_type_cluster_summary = ch_library_type_cluster_summary // path (or []) - for SAMPLE_SUMMARY
    library_type_result          = ch_library_type_result          // tuple(meta, file) - empty channel if this method is off, for LIBRARY_TYPE_CONSENSUS
    library_type_cluster_result  = ch_library_type_cluster_result  // tuple(meta, file) - empty channel if this method is off, for LIBRARY_TYPE_CONSENSUS
}
