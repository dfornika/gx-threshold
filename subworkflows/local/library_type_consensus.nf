//
// Library-type consensus: fuse the (up to 4) amplicon-vs-shotgun verdicts -
// LIBRARY_TYPE (fastp/fastplong QC JSON), LIBRARY_TYPE_CLUSTER (read
// clustering), LIBRARY_TYPE_ALIGNED (depth dispersion), LIBRARY_TYPE_PILEUP
// (read-position pileup) - into one verdict per sample via majority vote
// (LIBRARY_TYPE_CONSENSUS - bin/library_type_consensus.py), rather than only
// ever tagging `meta.library_type` from one method's own call. Runs after
// both LIBRARY_TYPE_REFERENCE_FREE and ALIGNMENT_BASED_LIBRARY_TYPE so all
// four methods (whichever are enabled) have already had a chance to run -
// see docs/testing.md for the vote rule and provisional fallback order.
//

include { LIBRARY_TYPE_CONSENSUS } from '../../modules/local/library_type_consensus/main'
include { tagMetaFromVerdict     } from './meta_utils'

workflow LIBRARY_TYPE_CONSENSUS_ANALYSIS {

    take:
    ch_clean_reads                   // tuple(meta, reads) - canonical meta/reads, to be tagged with meta.library_type
    ch_library_type_result           // tuple(meta, file) - LIBRARY_TYPE (fastp_json); empty channel if off
    ch_library_type_cluster_result   // tuple(meta, file) - LIBRARY_TYPE_CLUSTER; empty channel if off
    ch_library_type_aligned_result   // tuple(meta, file) - LIBRARY_TYPE_ALIGNED; empty channel if off
    ch_library_type_pileup_result    // tuple(meta, file) - LIBRARY_TYPE_PILEUP; empty channel if off
    outdir

    main:

    def ch_multiqc_files = channel.empty()

    // Each method's meta was captured at a different point in the pipeline
    // (see meta_utils.nf), so - as elsewhere in this codebase - join by
    // meta.id rather than trusting the whole meta map to match across
    // channels. Tag every call with its method name before mixing, so one
    // sample's calls can be grouped back together regardless of how many
    // of the four methods actually ran for it.
    def ch_calls_by_id = channel.empty()
        .mix(ch_library_type_result.map         { meta, file -> tuple(meta.id, file, 'fastp_json') })
        .mix(ch_library_type_cluster_result.map { meta, file -> tuple(meta.id, file, 'cluster') })
        .mix(ch_library_type_aligned_result.map { meta, file -> tuple(meta.id, file, 'aligned') })
        .mix(ch_library_type_pileup_result.map  { meta, file -> tuple(meta.id, file, 'pileup') })
        .groupTuple(by: 0) // tuple(id, [files], [methods])

    // remainder: true + the `?: []` defaults below: a sample with zero
    // library-type calls at all (every method skipped) must still flow
    // through with a "no_data" consensus, not silently disappear from
    // ch_clean_reads - a plain (inner) `.join()` would drop it entirely,
    // since it'd have no matching entry in ch_calls_by_id.
    def ch_calls_with_meta = ch_clean_reads
        .map { meta, _reads -> tuple(meta.id, meta) }
        .join(ch_calls_by_id, remainder: true)
        .map { _id, meta, files, methods -> tuple(meta, files ?: [], methods ?: []) }

    LIBRARY_TYPE_CONSENSUS(ch_calls_with_meta)
    ch_multiqc_files = ch_multiqc_files.mix(LIBRARY_TYPE_CONSENSUS.out.result.map { _meta, file -> file })

    // .ifEmpty([]): see subworkflows/local/species_id.nf for why.
    def ch_library_type_consensus_summary = LIBRARY_TYPE_CONSENSUS.out.result
        .map { _meta, file -> file }
        .collectFile(
            name: 'library_type_consensus_summary.tsv',
            storeDir: "${outdir}/library_type",
            sort: true,
            seed: "sample\tplatform\tverdict\tn_amplicon_votes\tn_shotgun_votes\tn_total_votes\tmethod\tnote\n"
        )
        .ifEmpty([])

    def ch_tagged_reads = tagMetaFromVerdict(ch_clean_reads, LIBRARY_TYPE_CONSENSUS.out.result, 'library_type')

    emit:
    reads                          = ch_tagged_reads   // tuple(meta, reads) - meta.library_type set from the consensus verdict
    multiqc_files                  = ch_multiqc_files
    library_type_consensus_summary = ch_library_type_consensus_summary // path (or []) - for SAMPLE_SUMMARY
}
