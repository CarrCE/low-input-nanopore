#!/usr/bin/env python3
"""Mask human sequence in a FASTQ, rescuing reads this pipeline positively
attributes to a non-human organism.

Two stages, deliberately separated:

  1. SCREEN.  NCBI's Human Read Removal Tool (HRRT) flags candidate reads by
     k-mer match.  HRRT is sensitive but blunt: used on its own it would have
     destroyed 49.4% of this study's *S. cerevisiae* reads and 92.0% of its
     *C. neoformans* reads, because a conserved region is a conserved region
     whichever genome it sits in.  We therefore use HRRT ONLY for its flag
     list and never for its own masking.

  2. RESCUE.  A flagged read is kept when `assignments.tsv.gz` -- the
     competitive attribution this whole study already rests on -- assigns it
     to an organism.  The reference panel contains no human, so "assigned"
     always means "assigned to something that is not human".

THE RESCUE MUST REQUIRE POSITIVE ATTRIBUTION.  It is written as "keep only if
some organism claims this read", never as "mask unless it failed to look
human".  The inversion is subtle, reads almost identically, and is backwards
for a privacy filter: it keeps every read that aligned to human but fell under
whatever threshold happened to be configured.  A prior implementation shipped
that inversion and it cost a full revision.

Scope of the mask
-----------------
Deciding *whether* to touch a read and deciding *what* to blank are separate
questions.  With H the union of human-aligned query intervals, T the union of
organism-aligned query intervals, and L the read length:

    keep = T \\ H          mask = [0, L) \\ keep

Only bases some organism positively accounts for, and human does not, survive.
That single expression needs no special cases:

    wholly human, no organism hit  ->  T empty, so the whole read is masked
    chimeric                       ->  masks the human part and any gap
    organism read, conserved region->  masks that region only

A read that is `assigned` and NOT chimeric is left completely untouched -- we
do not blank unaligned gaps in reads we have positively attributed, because
the flag was a k-mer coincidence and there is nothing to protect.

Chimeras
--------
A read that is part human and part microbial wins its competitive comparison
on the microbial half, so whole-read best-hit assignment structurally cannot
catch it and the human half would be released intact.  When the human-aligned
sequence that no organism accounts for reaches --chimera-min-bp, the read is
masked by the scope rule above regardless of its assignment.  Rare -- 259
times in 60.4 M reads -- and precisely the case that matters.

Invariants
----------
Quality strings are never touched, trimmed or recomputed.  Output length equals
input length for every record; a fully-masked read is still written, as an
all-N record.  Read counts and base counts are therefore identical before and
after, which is why masking (rather than deletion) leaves this study's
published per-femtogram and enrichment figures undisturbed.

Outputs
-------
    <out-fastq>       masked reads, the only artifact that gets deposited
    <out-manifest>    per-read audit trail: what was masked and why
    <out-stats>       JSON counts, thresholds and provenance
    <out-masked-ids>  IDs of reads carrying any mask
"""

import argparse
import gzip
import io
import json
import os
import sys
from collections import Counter


# --------------------------------------------------------------------------
# interval arithmetic.  Half-open [start, end), query coordinates.
# --------------------------------------------------------------------------

def merge(intervals):
    """Union of a list of half-open intervals, sorted and coalesced."""
    if not intervals:
        return []
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            if end > out[-1][1]:
                out[-1] = (out[-1][0], end)
        else:
            out.append((start, end))
    return out


def subtract(a, b):
    """a \\ b, both assumed already merged."""
    if not b:
        return list(a)
    out = []
    for start, end in a:
        cur = start
        for bs, be in b:
            if be <= cur or bs >= end:
                continue
            if bs > cur:
                out.append((cur, min(bs, end)))
            cur = max(cur, be)
            if cur >= end:
                break
        if cur < end:
            out.append((cur, end))
    return out


def complement(intervals, length):
    """[0, length) \\ intervals, intervals assumed merged."""
    return subtract([(0, length)], intervals)


def total_len(intervals):
    return sum(end - start for start, end in intervals)


# --------------------------------------------------------------------------
# input parsing
# --------------------------------------------------------------------------

def _open(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_flagged_ids(path):
    """Read HRRT's removed-spots output.

    HRRT's `-u` file is a FASTQ of the reads it flagged, NOT an ID list.  Read
    IDs are located BY POSITION -- every 4th line -- and never by a leading
    '@', because a quality string may legitimately begin with '@' (Phred 31).
    Counting sequence lines, '+' lines and quality strings as read IDs inflated
    a prior implementation's human counts roughly fourfold and filled its
    report with nonsense.

    A plain one-ID-per-line list is also accepted, for the test fixtures.
    """
    ids = set()
    with _open(path) as fh:
        first = fh.readline()
        if not first:
            return ids
        fh.seek(0)
        if first.startswith("@"):
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    if not line.startswith("@"):
                        raise RuntimeError(
                            f"{path}: line {i + 1} should be a FASTQ header but "
                            f"does not start with '@'; the file is not 4-line "
                            f"periodic and positional parsing would be wrong"
                        )
                    ids.add(line[1:].split()[0])
        else:
            for line in fh:
                line = line.strip()
                if line:
                    ids.add(line.split()[0])
    return ids


def load_assignments(path, flagged):
    """read_id -> (call, organism, role) for flagged reads only.

    Streams; only flagged reads are retained, so memory is bounded by the flag
    count (tens of thousands) rather than by the read count (tens of millions).
    """
    wanted = {}
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_id = header.index("read_id")
            i_org = header.index("organism")
            i_call = header.index("call")
            i_role = header.index("role")
        except ValueError as exc:
            raise RuntimeError(f"{path}: unexpected assignments header: {exc}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rid = f[i_id]
            if rid in flagged:
                wanted[rid] = (f[i_call], f[i_org], f[i_role])
    return wanted


def load_human_intervals(paf, flagged):
    """read_id -> merged human query intervals, from minimap2 PAF.

    Identity is deliberately NOT computed here and no identity threshold is
    applied: the rescue decision comes from `assignments.tsv.gz`, so any human
    alignment at all contributes to H.

    If anyone later adds an identity filter, do NOT derive it from PAF columns
    10 and 11 -- they are equal when minimap2 runs without `-c`, which makes
    any threshold built on them silently inert (a genuinely 99.92%-identical
    read reported 510/510).  Use the `de:f` tag with `-c`, or `dv:f` without.
    """
    out = {}
    with _open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            rid = f[0]
            if rid not in flagged:
                continue
            out.setdefault(rid, []).append((int(f[2]), int(f[3])))
    return {k: merge(v) for k, v in out.items()}


def load_target_intervals(bam, flagged):
    """read_id -> merged organism-aligned query intervals, from the BAM.

    Query intervals come from the CIGAR, so a read aligning to two organisms
    or in two pieces contributes both spans.  Secondary and supplementary
    alignments are included: for masking we want every base any organism can
    account for, not just the primary call.
    """
    import pysam

    # Accept SAM as well as BAM: the committed test fixture ships a text SAM so
    # its expected query intervals can be read without tooling.
    mode = "r" if bam.endswith((".sam", ".sam.gz")) else "rb"
    out = {}
    with pysam.AlignmentFile(bam, mode, check_sq=False) as fh:
        for aln in fh.fetch(until_eof=True):
            if aln.is_unmapped:
                continue
            rid = aln.query_name
            if rid not in flagged:
                continue
            start = aln.query_alignment_start
            end = aln.query_alignment_end
            if end > start:
                out.setdefault(rid, []).append((start, end))
    return {k: merge(v) for k, v in out.items()}


# --------------------------------------------------------------------------
# the decision
# --------------------------------------------------------------------------

def decide(length, human, target, call, chimera_min_bp):
    """Return (mask_intervals, reason, chimeric).

    `human` and `target` are merged interval lists; `call` is the value from
    assignments.tsv.gz.  An empty mask list means the read is released intact.
    """
    assigned = call == "assigned"

    # A chimera is human sequence that NO organism accounts for, on a read
    # that some organism nonetheless partly explains.  Subtracting target from
    # human first is what makes this specific: a conserved region covered by an
    # organism alignment is not human-exclusive and does not trip the rule.
    #
    # `target` must be non-empty.  Without that clause a wholly-human read
    # trivially satisfies the length test and gets counted as a chimera, which
    # would report every human read as one -- inflating the count from a few
    # hundred to a few hundred thousand and destroying the only statistic that
    # tells us the chimera rule is working.  The masking decision is unaffected
    # either way; this guards the *statistic*.
    human_exclusive = total_len(subtract(human, target))
    chimeric = bool(target) and human_exclusive >= chimera_min_bp

    if assigned and not chimeric:
        return [], "assigned:" + call, False

    keep = subtract(target, human)
    mask = complement(keep, length)

    if not target:
        reason = "no organism alignment" if not human else "human, no organism alignment"
    elif chimeric:
        reason = f"chimeric: {human_exclusive} human-exclusive bp"
    else:
        reason = f"not positively attributed ({call})"
    return mask, reason, chimeric


def apply_mask(seq, mask):
    """Blank the masked intervals with N, preserving length exactly."""
    if not mask:
        return seq
    buf = bytearray(seq, "ascii")
    for start, end in mask:
        buf[start:end] = b"N" * (end - start)
    return buf.decode("ascii")


# --------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fastq", required=True, help="input FASTQ (.gz ok)")
    p.add_argument("--flagged", required=True,
                   help="HRRT removed-spots FASTQ, or a plain read-ID list")
    p.add_argument("--assignments", required=True, help="assignments.tsv.gz")
    p.add_argument("--human-paf", required=True,
                   help="PAF of the flagged reads against the human reference")
    p.add_argument("--bam", required=True,
                   help="organism BAM (may be pre-filtered to the flagged reads)")
    p.add_argument("--sample-id", required=True)
    p.add_argument("--out-fastq", required=True, help="'-' for stdout")
    p.add_argument("--out-manifest", required=True)
    p.add_argument("--out-stats", required=True)
    p.add_argument("--out-masked-ids", required=True)
    p.add_argument("--chimera-min-bp", type=int, default=150)
    p.add_argument("--human-reference", default="", help="recorded in the stats")
    p.add_argument("--scrubber-version", default="", help="recorded in the stats")
    args = p.parse_args()

    flagged = parse_flagged_ids(args.flagged)
    sys.stderr.write(f"[mask] {args.sample_id}: {len(flagged):,} reads flagged by HRRT\n")

    assignments = load_assignments(args.assignments, flagged)
    human = load_human_intervals(args.human_paf, flagged)
    target = load_target_intervals(args.bam, flagged)
    sys.stderr.write(
        f"[mask] {args.sample_id}: {len(assignments):,} flagged reads found in "
        f"assignments, {len(human):,} with human alignments, "
        f"{len(target):,} with organism alignments\n"
    )

    stats = Counter()
    reasons = Counter()
    calls = Counter()
    masked_len_hist = []
    kept_len_hist = []

    out_is_stdout = args.out_fastq == "-"
    if out_is_stdout:
        out_fh = sys.stdout
    elif args.out_fastq.endswith(".gz"):
        out_fh = io.TextIOWrapper(gzip.open(args.out_fastq, "wb", compresslevel=6))
    else:
        out_fh = open(args.out_fastq, "w")

    man_fh = _open(args.out_manifest, "wt") if args.out_manifest.endswith(".gz") \
        else open(args.out_manifest, "w")
    ids_fh = open(args.out_masked_ids, "w")

    man_fh.write("read_id\tread_len\tmasked_bases\tmasked_frac\twhole_read\t"
                 "chimeric\tcall\torganism\tintervals\treason\n")

    seen_flagged = set()
    with _open(args.fastq) as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            seq = fh.readline().rstrip("\n")
            plus = fh.readline()
            qual = fh.readline().rstrip("\n")
            if not qual:
                raise RuntimeError(f"{args.fastq}: truncated record at {head[:60]!r}")
            if not head.startswith("@"):
                raise RuntimeError(
                    f"{args.fastq}: expected a header, got {head[:60]!r}; the "
                    f"file is not 4-line periodic"
                )
            rid = head[1:].split()[0]
            length = len(seq)

            stats["reads"] += 1
            stats["bases"] += length

            if rid in flagged:
                seen_flagged.add(rid)
                call, organism, _role = assignments.get(rid, ("absent", "", ""))
                calls[call] += 1
                mask, reason, chimeric = decide(
                    length, human.get(rid, []), target.get(rid, []),
                    call, args.chimera_min_bp,
                )
                if mask:
                    masked_bases = total_len(mask)
                    whole = 1 if masked_bases == length else 0
                    seq = apply_mask(seq, mask)
                    stats["masked_reads"] += 1
                    stats["masked_bases"] += masked_bases
                    stats["fully_masked_reads" if whole else "partially_masked_reads"] += 1
                    if chimeric:
                        stats["chimeric_masks"] += 1
                        stats["chimeric_human_exclusive_bp"] += masked_bases
                    reasons[reason.split(":")[0]] += 1
                    masked_len_hist.append(length)
                    man_fh.write(
                        f"{rid}\t{length}\t{masked_bases}\t{masked_bases / length:.4f}\t"
                        f"{whole}\t{1 if chimeric else 0}\t{call}\t{organism}\t"
                        f"{','.join(f'{s}-{e}' for s, e in mask)}\t{reason}\n"
                    )
                    ids_fh.write(rid + "\n")
                else:
                    stats["rescued_reads"] += 1
                    reasons["rescued"] += 1
                    kept_len_hist.append(length)

            if len(seq) != length:
                raise RuntimeError(
                    f"{rid}: masking changed the read length "
                    f"({length} -> {len(seq)}); this must never happen"
                )
            out_fh.write(f"{head}{seq}\n{plus}{qual}\n")

    if not out_is_stdout:
        out_fh.close()
    man_fh.close()
    ids_fh.close()

    # ---- guards ----------------------------------------------------------
    # If HRRT and the FASTQ ever key their read IDs differently, every
    # comparison above silently becomes meaningless while all the counts stay
    # entirely plausible.  Fail loudly instead.
    if flagged and not seen_flagged:
        raise RuntimeError(
            "HRRT flagged read IDs do not intersect the FASTQ read IDs at all "
            "-- the two files are keyed differently and every downstream count "
            "would be meaningless"
        )
    missing = len(flagged) - len(seen_flagged)
    if missing:
        sys.stderr.write(
            f"[mask] WARNING {args.sample_id}: {missing:,} flagged IDs were not "
            f"seen in the FASTQ\n"
        )

    def pct(n, d):
        return round(100.0 * n / d, 4) if d else 0.0

    out = {
        "sample_id": args.sample_id,
        "rule": {
            "screen": "NCBI sra-human-scrubber (HRRT), flag list only",
            "rescue": "keep iff assignments.tsv.gz call == 'assigned'",
            "mask_scope": "[0,L) \\ (target \\ human)",
            "chimera_min_bp": args.chimera_min_bp,
            "identity_threshold": None,
        },
        "provenance": {
            "human_reference": args.human_reference,
            "scrubber_version": args.scrubber_version,
            "assignments": os.path.basename(args.assignments),
        },
        "reads": stats["reads"],
        "bases": stats["bases"],
        "hrrt_flagged": len(flagged),
        "hrrt_flagged_pct": pct(len(flagged), stats["reads"]),
        "rescued_reads": stats["rescued_reads"],
        "masked_reads": stats["masked_reads"],
        "masked_reads_pct": pct(stats["masked_reads"], stats["reads"]),
        "masked_bases": stats["masked_bases"],
        "masked_bases_pct": pct(stats["masked_bases"], stats["bases"]),
        "fully_masked_reads": stats["fully_masked_reads"],
        "partially_masked_reads": stats["partially_masked_reads"],
        "chimeric_masks": stats["chimeric_masks"],
        "chimeric_human_exclusive_bp": stats["chimeric_human_exclusive_bp"],
        "flagged_call_breakdown": dict(calls),
        "mask_reason_breakdown": dict(reasons),
        # Adaptive sampling truncates reads, so a human read may be only a few
        # hundred bases.  These distributions are what any future length
        # threshold must be justified against -- do not add one without
        # looking at them first.
        "flagged_read_length": {
            "masked_median": _median(masked_len_hist),
            "masked_min": min(masked_len_hist) if masked_len_hist else 0,
            "masked_max": max(masked_len_hist) if masked_len_hist else 0,
            "rescued_median": _median(kept_len_hist),
            "rescued_min": min(kept_len_hist) if kept_len_hist else 0,
            "rescued_max": max(kept_len_hist) if kept_len_hist else 0,
        },
    }
    with open(args.out_stats, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")

    sys.stderr.write(
        f"[mask] {args.sample_id}: {stats['reads']:,} reads, "
        f"{len(flagged):,} flagged ({out['hrrt_flagged_pct']}%), "
        f"{stats['rescued_reads']:,} rescued, "
        f"{stats['masked_reads']:,} masked ({out['masked_reads_pct']}%), "
        f"{stats['chimeric_masks']:,} chimeric\n"
    )


def _median(xs):
    if not xs:
        return 0
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) // 2


if __name__ == "__main__":
    main()
