#!/usr/bin/awk -f
#
# Classify a sample's library prep as amplicon vs. shotgun from a SAM stream
# (piped directly from minimap2), using a streaming index of dispersion
# (variance/mean) of per-base alignment depth.
#
# Rationale: amplicon libraries repeatedly re-sequence the same small set of
# reference regions (tiling primer targets), producing wildly uneven depth
# across the genome, while shotgun libraries sample depth close to uniformly
# (a Poisson process, index of dispersion ~1). This is the same statistic
# validated (in batch, via `samtools depth -a`) during development - see
# docs/testing.md - reimplemented here as a running update so it can be
# computed incrementally as reads stream in, and so a confident amplicon
# verdict can stop the alignment early rather than waiting for the whole
# input.
#
# An earlier version of this script tried to reformulate the statistic as a
# per-read Bernoulli sequential probability ratio test ("did this read hit an
# already-covered bin or a new one?"). That discretization didn't survive
# contact with real data: fixing a false-positive on long-read bacterial
# shotgun data (birthday-paradox chance bin collisions early in the stream)
# required coarsening the bins, which then broke detection on real
# SARS-CoV-2 amplicon tiling data (whose ~90 distinct amplicons span most of
# a coarse bin grid, so "new bin" events - genuine tiling behaviour - looked
# like evidence *against* amplicon). Tracking actual per-base depth directly,
# rather than discretizing into repeat/new-bin events, avoids that tension
# entirely.
#
# Depth is tracked on a fixed-size grid (`grid_size` bp per cell, plain
# position-based averaging, not true per-base resolution) purely so the
# per-read update cost is bounded - this is a performance/memory
# consideration only, not a statistical discretization choice like the old
# bins were: any reasonably fine grid gives essentially the same variance/mean
# ratio. Reference-consumed length is read from the CIGAR string so long
# indel-heavy Nanopore alignments are handled correctly.
#
# Only one alignment record per fragment is scored (unmapped/secondary/
# supplementary/second-in-pair records are skipped): for paired-end input,
# mates of the same fragment land near each other by construction, which
# would otherwise inflate apparent depth unevenness independent of library
# type.
#
# `min_reads` is a burn-in: an early (mid-stream) amplicon call is only made
# once at least this many reads have been scored, since it's a noisy
# estimate from very few observations. Stopping is one-sided - a confident
# amplicon verdict can be reached early (the signature signal), but a
# shotgun verdict requires exhausting the read stream (or `max_reads`) to
# confirm the absence of that signal, matching the observed asymmetry:
# amplicon shows a fast, distinctive signal; shotgun is "no signal", which
# needs more data to be sure of. At end of stream there's no more data
# coming regardless, so the only floor applied there is `MIN_READS_EOF` (a
# small fixed constant, not the mid-stream burn-in) - just enough to rule
# out a call from a literal handful of reads.
#
# Expected -v args: sample, platform, genome_length, grid_size, threshold,
# min_reads, max_reads
#
# Output (one line, tab-separated): sample platform verdict n_reads_used
# index_of_dispersion method

function hasflag(f, bit) {
    return int(f / bit) % 2
}

# Sum of CIGAR operations that consume the reference (M/D/N/=/X).
function reflen(cigar,    len, n, i, c, num) {
    len = 0
    num = ""
    n = length(cigar)
    for (i = 1; i <= n; i++) {
        c = substr(cigar, i, 1)
        if (c ~ /[0-9]/) {
            num = num c
        } else {
            if (c == "M" || c == "D" || c == "N" || c == "=" || c == "X") {
                len += num + 0
            }
            num = ""
        }
    }
    return len
}

function index_of_dispersion(    mean, var) {
    if (n_cells <= 0) return 0
    mean = sum_depth / n_cells
    if (mean <= 0) return 0
    var = sumsq_depth / n_cells - mean * mean
    return var / mean
}

function report() {
    printf "%s\t%s\t%s\t%d\t%.4g\t%s\n", sample, platform, verdict, n_reads, index_of_dispersion(), method
}

BEGIN {
    MIN_READS_EOF = 20

    if (grid_size == "" || grid_size < 1) grid_size = 50
    if (threshold == "")  threshold  = 2.5
    if (min_reads == "")  min_reads  = 200
    if (max_reads == "")  max_reads  = 20000
    if (genome_length == "" || genome_length < 1) genome_length = 1

    n_cells = int(genome_length / grid_size)
    if (n_cells < 1) n_cells = 1

    sum_depth   = 0
    sumsq_depth = 0
    n_reads     = 0
    verdict     = "inconclusive"
    method      = "eof"
}

/^@/ { next }

{
    flag = $2
    if (hasflag(flag, 4))    next  # unmapped
    if (hasflag(flag, 256))  next  # secondary alignment
    if (hasflag(flag, 2048)) next  # supplementary alignment
    if (hasflag(flag, 128))  next  # second-in-pair: score one obs per fragment

    rname = $3
    pos   = $4
    cigar = $6
    if (rname == "*" || cigar == "*") next

    len = reflen(cigar)
    if (len < 1) next

    start_cell = int((pos - 1) / grid_size)
    end_cell   = int((pos - 1 + len - 1) / grid_size)

    for (c = start_cell; c <= end_cell; c++) {
        key = rname "_" c
        d = (key in depth) ? depth[key] : 0
        sumsq_depth += 2 * d + 1
        sum_depth   += 1
        depth[key] = d + 1
    }

    n_reads++

    if (n_reads >= min_reads && index_of_dispersion() >= threshold) {
        verdict = "amplicon"; method = "threshold_stop"; report(); exit
    }
    if (n_reads >= max_reads) {
        verdict = (index_of_dispersion() >= threshold) ? "amplicon" : "shotgun"
        method = "max_reads"; report(); exit
    }
}

END {
    if (verdict == "inconclusive") {
        if (n_reads < MIN_READS_EOF) {
            method = "insufficient_data"
        } else {
            verdict = (index_of_dispersion() >= threshold) ? "amplicon" : "shotgun"
            method = "eof"
        }
        report()
    }
}
