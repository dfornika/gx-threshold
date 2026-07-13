//
// Small shared helpers for working with `meta` as it accumulates simple
// classification tags (`meta.library_type`, `meta.composition`, ...) while
// flowing through this pipeline's subworkflows - see workflows/threshold.nf.
//
// Important: because `meta` gains fields at different points for different
// channels (e.g. `ch_species_id_rows`' meta was captured before composition
// tagging happened; `ch_clean_reads`' meta has it by the time reference
// fetch runs), two `tuple(meta, ...)` channels almost never have
// byte-for-byte identical `meta` maps for the same sample, even though they
// agree on `meta.id`. Never `.combine(other, by: 0)` or `.join(other)`
// directly on the whole map for that reason - it silently produces zero
// matches instead of an error. Always join by `meta.id` instead (both
// helpers below do this), and decide explicitly which side's `meta` to keep.

// Fold a classification module's verdict into `meta` under `key`, by
// joining its one-line-per-sample result back onto `ch_reads` by sample id.
// Every classify_*.py script in this pipeline writes verdict as the 3rd
// tab-separated column (sample, platform, verdict, ...) - that's the
// convention this relies on.
//
// `remainder: true` + the `?: 'no_data'` default: a classifier can
// legitimately produce zero rows for a sample (e.g. sourmash gather finding
// no hit at all against a database, as happens for the tiny synthetic
// `test`-profile fixtures - see docs/testing.md). A plain (inner) `.join()`
// would silently drop that sample from `ch_reads` entirely rather than just
// leaving `meta[key]` unset - the same "missing output looks like nothing
// happened" footgun documented elsewhere in this codebase (see
// subworkflows/local/species_id.nf) - which is invisible until something
// downstream actually depends on every sample still being present.
def tagMetaFromVerdict(ch_reads, ch_result, key) {
    return ch_reads
        .map { meta, reads -> tuple(meta.id, meta, reads) }
        .join(ch_result.map { meta, file -> tuple(meta.id, file.text.trim().split('\t')[2]) }, remainder: true)
        .map { _id, meta, reads, verdict -> tuple(meta + [(key): (verdict ?: 'no_data')], reads) }
}

// Join two `tuple(meta, value)` channels by `meta.id`, keeping `ch_a`'s
// `meta` (assumed the more currently-enriched side) in the result:
// tuple(meta, a_value, b_value).
def joinByMetaId(ch_a, ch_b) {
    return ch_a
        .map { meta, a_value -> tuple(meta.id, meta, a_value) }
        .join(ch_b.map { meta, b_value -> tuple(meta.id, b_value) })
        .map { _id, meta, a_value, b_value -> tuple(meta, a_value, b_value) }
}
