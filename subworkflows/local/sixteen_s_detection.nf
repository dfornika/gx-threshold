//
// 16S rRNA amplicon detection - the specific "which kind of amplicon is
// this" question left open by LIBRARY_TYPE_REFERENCE_FREE (amplicon vs.
// shotgun) and SPECIES_COMPOSITION_ANALYSIS (pure culture vs. metagenomic):
// an amplicon sample whose composition breadth came back inconclusive
// against the whole-genome species-ID database is either 16S rRNA
// metagenomic profiling, or some other targeted marker-gene amplicon
// (hsp65, an AMR gene panel, ...) - grouped simply as "other" rather than
// further sub-classified, since distinguishing those isn't needed yet.
//
// Only runs for samples where meta.composition == 'inconclusive' - a
// sample with good whole-genome recovery already has its answer (pure
// culture or shotgun metagenomic), and 16S profiling is always PCR-
// amplified, so shotgun samples are never candidates either. See
// docs/testing.md for why the gate is composition-based rather than also
// requiring meta.library_type == 'amplicon': LIBRARY_TYPE_CLUSTER (the tag
// source) turned out unreliable on genuinely multi-organism samples, so
// requiring its agreement would have missed a real 16S sample in testing.
//
// Alignment-based (against a dedicated 16S reference database - see
// data/dev_16s_db/README.md for why sourmash/mash don't work here, and a
// small subsample of reads (SUBSAMPLE_READS_HEAD), since this is a coarse
// yes/no question, not a composition breakdown, and doesn't need many
// reads to answer confidently.
//

include { SUBSAMPLE_READS_HEAD  } from '../../modules/local/subsample_reads_head/main'
include { ALIGN_READS           } from '../../modules/local/align_reads/main'
include { CLASSIFY_16S_AMPLICON } from '../../modules/local/classify_16s_amplicon/main'
include { tagMetaFromVerdict    } from './meta_utils'

workflow SIXTEEN_S_DETECTION {

    take:
    ch_clean_reads   // tuple(meta, reads) - meta.composition already set by SPECIES_COMPOSITION_ANALYSIS
    outdir

    main:

    def ch_multiqc_files = channel.empty()
    def ch_tagged_reads  = ch_clean_reads
    def ch_sixteen_s_summary = channel.value([])

    if (!params.skip_sixteen_s_detection) {
        if (!params.sixteen_s_db) {
            error("16S detection is enabled but --sixteen_s_db was not set. Provide a 16S rRNA reference FASTA or run with --skip_sixteen_s_detection.")
        }

        def ch_branched = ch_clean_reads.branch { meta, _reads ->
            needs_check: meta.composition == 'inconclusive'
            skip: true
        }

        def ch_sixteen_s_db = channel.value(file(params.sixteen_s_db, checkIfExists: true))

        SUBSAMPLE_READS_HEAD(ch_branched.needs_check, params.sixteen_s_max_reads)
        def ch_reads_with_16s_db = SUBSAMPLE_READS_HEAD.out.reads.combine(ch_sixteen_s_db)
        ALIGN_READS(ch_reads_with_16s_db)
        CLASSIFY_16S_AMPLICON(ALIGN_READS.out.sam)
        ch_multiqc_files = ch_multiqc_files.mix(CLASSIFY_16S_AMPLICON.out.result.map { _meta, file -> file })

        // .ifEmpty([]): see reference_genome.nf - collectFile emits nothing
        // if no sample needed the 16S check this run (all composition-
        // conclusive), so SAMPLE_SUMMARY still gets a usable value.
        ch_sixteen_s_summary = CLASSIFY_16S_AMPLICON.out.result
            .map { _meta, file -> file }
            .collectFile(
                name: 'sixteen_s_detection_summary.tsv',
                storeDir: "${outdir}/library_type",
                sort: true,
                seed: "sample\tplatform\tverdict\tn_reads\tn_passed\tpassed_frac\tnote\n"
            )
            .ifEmpty([])

        def ch_checked_reads = tagMetaFromVerdict(ch_branched.needs_check, CLASSIFY_16S_AMPLICON.out.result, 'sixteen_s')
        ch_tagged_reads = ch_branched.skip.mix(ch_checked_reads)
    }

    emit:
    reads             = ch_tagged_reads   // tuple(meta, reads) - meta.sixteen_s set for samples the check ran on
    multiqc_files     = ch_multiqc_files
    sixteen_s_summary = ch_sixteen_s_summary // path (or []) - for SAMPLE_SUMMARY
}
