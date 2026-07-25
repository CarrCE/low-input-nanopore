# Filter a FASTQ by mean read quality, matching how ONT tooling defines it.
#
# The mean read qscore is NOT the arithmetic mean of the per-base Phred values.
# It is computed in error-probability space and converted back:
#
#     Q_read = -10 * log10( mean_i( 10^(-Q_i/10) ) )
#
# This is what dorado writes into the `qs:f:` header tag and what fastcat (and
# therefore wf-metagenomics `--min_read_qual`) computes. Averaging Phred scores
# arithmetically -- which is what `seqkit seq --min-qual` does -- gives a
# systematically higher number and would not be comparable to the prior-study
# reanalyses this filter exists to match.
#
# The `qs:f:` tag is used when present (it is, for dorado-basecalled reads) and
# recomputed from the quality string otherwise, so the same threshold can be
# applied to this study's reads and to published datasets from older basecallers.
#
# Usage:
#   awk -v MINQ=10 -f filter_by_qscore.awk in.fastq > out.fastq
#   awk -v MINQ=10 -v STATS=stats.txt -f filter_by_qscore.awk in.fastq > out.fastq

BEGIN {
    if (MINQ == "") MINQ = 0
    # Phred+33 lookup: ASCII code -> 10^(-Q/10)
    for (i = 33; i < 127; i++) {
        c = sprintf("%c", i)
        ord[c] = i - 33
        perr[c] = 10 ^ (-(i - 33) / 10.0)
    }
    kept = 0; total = 0; kept_bases = 0; total_bases = 0
}

NR % 4 == 1 {
    hdr = $0
    qs = ""
    # Header tags are tab-separated in dorado output: e.g. "qs:f:26.7872"
    if (match(hdr, /qs:f:[0-9.]+/)) {
        qs = substr(hdr, RSTART + 5, RLENGTH - 5) + 0
        have_qs = 1
    } else {
        have_qs = 0
    }
    next
}
NR % 4 == 2 { seq = $0; next }
NR % 4 == 3 { plus = $0; next }
NR % 4 == 0 {
    qual = $0
    total++
    total_bases += length(seq)

    if (!have_qs) {
        n = length(qual)
        if (n == 0) {
            qs = 0
        } else {
            s = 0
            for (i = 1; i <= n; i++) s += perr[substr(qual, i, 1)]
            m = s / n
            qs = (m > 0) ? -10 * log(m) / log(10) : 0
        }
    }

    if (qs + 0 >= MINQ + 0) {
        print hdr "\n" seq "\n" plus "\n" qual
        kept++
        kept_bases += length(seq)
    }
}

END {
    # %.0f, never %d. mawk (the default awk in Debian, and so in this image)
    # formats %d through a 32-bit int and saturates at 2147483647, which silently
    # truncated every base count on a real run -- 2.1 Gbases where the true value
    # was tens of Gbases. Counters are held as doubles, so %.0f prints them
    # exactly up to 2^53.
    if (STATS != "") {
        printf "min_qscore\t%s\n", MINQ                  > STATS
        printf "reads_in\t%.0f\n",  total                > STATS
        printf "reads_kept\t%.0f\n", kept                > STATS
        printf "bases_in\t%.0f\n",  total_bases          > STATS
        printf "bases_kept\t%.0f\n", kept_bases          > STATS
        printf "read_pass_fraction\t%.6f\n",  (total ? kept / total : 0)             > STATS
        printf "base_pass_fraction\t%.6f\n",  (total_bases ? kept_bases / total_bases : 0) > STATS
    }
    printf "[filter] MINQ=%s kept %.0f/%.0f reads (%.4f), %.0f/%.0f bases (%.4f)\n",
           MINQ, kept, total, (total ? kept / total : 0),
           kept_bases, total_bases, (total_bases ? kept_bases / total_bases : 0) > "/dev/stderr"
}
