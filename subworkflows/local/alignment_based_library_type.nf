//
// Two more amplicon-vs-shotgun approaches under evaluation (see
// docs/testing.md), both needing a fetched reference genome so they run
// after REFERENCE_GENOME rather than alongside LIBRARY_TYPE_REFERENCE_FREE:
//
// - LIBRARY_TYPE_ALIGNED: streaming index of dispersion of per-base depth
//   against the fetched reference. Platform-unified, unlike LIBRARY_TYPE.
// - LIBRARY_TYPE_PILEUP: aligned read/fragment start-end position pileup -
//   does this read start (and end) at the same place as many others, the
//   same positional signature real duplicate-marking tools use to flag PCR
//   duplicates.
//
// Neither currently feeds a meta tag (only LIBRARY_TYPE_CLUSTER does, in
// LIBRARY_TYPE_REFERENCE_FREE - see docs/testing.md for why); these remain
// comparison/validation outputs for now.
//

include { LIBRARY_TYPE_ALIGNED } from '../../modules/local/library_type_aligned/main'
include { ALIGN_READS          } from '../../modules/local/align_reads/main'
include { LIBRARY_TYPE_PILEUP  } from '../../modules/local/library_type_pileup/main'
include { joinByMetaId         } from './meta_utils'

workflow ALIGNMENT_BASED_LIBRARY_TYPE {

    take:
    ch_clean_reads      // tuple(meta, reads)
    ch_sample_reference // tuple(meta, reference fasta) - from REFERENCE_GENOME; empty channel if that stage is off
    outdir

    main:

    def ch_multiqc_files = channel.empty()

    if (!params.skip_reference_genome_fetch) {
        // Not a plain `.combine(by: 0)`: ch_sample_reference's meta was
        // captured before composition tagging happened, so it no longer
        // matches ch_clean_reads' meta exactly - see meta_utils.nf.
        def ch_reads_with_reference = joinByMetaId(ch_clean_reads, ch_sample_reference)

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

        ALIGN_READS(ch_reads_with_reference)
        LIBRARY_TYPE_PILEUP(ALIGN_READS.out.sam)
        ch_multiqc_files = ch_multiqc_files.mix(LIBRARY_TYPE_PILEUP.out.result.map { _meta, file -> file })

        LIBRARY_TYPE_PILEUP.out.result
            .map { _meta, file -> file }
            .collectFile(
                name: 'library_type_pileup_summary.tsv',
                storeDir: "${outdir}/library_type",
                sort: true,
                seed: "sample\tplatform\tverdict\tn_reads\tlargest_pileup_frac\tpiled_frac\teffective_signatures_frac\tnote\n"
            )
    }

    emit:
    multiqc_files = ch_multiqc_files
}
