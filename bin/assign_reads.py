#!/usr/bin/env python3
"""
Assign each read to one organism in the reference set, competitively.

Every read is mapped once against a single combined index holding the carrier,
the carrier-derived contaminant, and every community member. For each read we
compare the best alignment score (AS) achieved against each *organism* and
assign the read to the winner only when it beats the runner-up by a margin.
Reads whose best two organisms are effectively tied are reported as ambiguous
rather than silently awarded to one of them.

Why this matters for this study
-------------------------------
The lowinput_s1 community contains E. coli B-1109, and the lambda carrier
carries over E. coli K-12 from its production host. The two strains differ
(4,925,141 bp vs 4,641,652 bp) but share a highly similar core genome. A
sequential "map to K-12, delete what sticks" subtraction therefore removes the
community's own E. coli along with the contaminant. Competitive assignment
instead partitions reads into three honest bins -- B-1109-specific,
K-12-specific, and core-genome-ambiguous -- so the contaminant contribution can
be estimated from the strain-specific fraction rather than assumed.

Two assignment rules are available over the same alignments:

  competitive  the rule described above -- best organism wins by margin.
  sequential   the classic subtraction chain, kept so the two can be compared
               directly: any read aligning to the carrier is removed, then any
               read aligning to the contaminant, then the survivors are assigned
               to their best community hit. This is what over-removes the
               community's own E. coli, and running it against the identical
               alignments makes that loss measurable rather than assumed.

Input is a qname-grouped (unsorted, straight-from-minimap2) BAM/SAM so all
alignments for a read arrive together and nothing has to be held in memory.

Outputs
    <prefix>.assignments.tsv.gz  per-read: read_id, organism, call, AS, margin, length
    <prefix>.counts.tsv          per-organism: reads, read_bases, aligned_bases
    <prefix>.readlengths.tsv.gz  per-read length by assignment class, for the
                                 adaptive-sampling ejection-signature analysis
"""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict

import pysam


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bam", help="qname-grouped BAM/SAM from minimap2 (unsorted)")
    p.add_argument("--contig-map", required=True,
                   help="TSV: contig<TAB>organism<TAB>role")
    p.add_argument("--prefix", required=True, help="output prefix")
    p.add_argument("--min-margin-abs", type=int, default=10,
                   help="minimum absolute AS margin over the runner-up organism")
    p.add_argument("--min-margin-frac", type=float, default=0.01,
                   help="minimum AS margin as a fraction of the winning AS")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="discard alignments below this MAPQ before scoring")
    p.add_argument("--mode", choices=["competitive", "sequential"],
                   default="competitive",
                   help="assignment rule; see module docstring")
    p.add_argument("--subtract-order", default="carrier,contaminant",
                   help="sequential mode only: roles to subtract, in order")
    p.add_argument("--min-subtract-as", type=int, default=0,
                   help="sequential mode only: minimum AS for a read to count as "
                        "belonging to a subtracted role")
    return p.parse_args()


def load_contig_map(path):
    contig2org, contig2role = {}, {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig, organism, role = parts[0], parts[1], parts[2]
            contig2org[contig] = organism
            contig2role[organism] = role
    if not contig2org:
        sys.exit(f"error: no contigs read from {path}")
    return contig2org, contig2role


def alignment_score(aln):
    """minimap2 AS tag; fall back to aligned length when absent."""
    try:
        return int(aln.get_tag("AS"))
    except KeyError:
        return aln.query_alignment_length or 0


def read_length(aln):
    """
    Full read length, independent of clipping.

    Secondary alignments carry no SEQ, and supplementary alignments are hard
    clipped, so infer_read_length() (which counts hard clips) is the only value
    that is consistent across every record for a read.
    """
    return aln.infer_read_length() or aln.query_length or 0


def main():
    args = parse_args()
    contig2org, org2role = load_contig_map(args.contig_map)

    counts = defaultdict(lambda: {"reads": 0, "read_bases": 0, "aligned_bases": 0})

    assign_fh = gzip.open(f"{args.prefix}.assignments.tsv.gz", "wt")
    len_fh = gzip.open(f"{args.prefix}.readlengths.tsv.gz", "wt")
    assign_fh.write("read_id\torganism\tcall\trole\tas_best\tas_runnerup\tmargin\tread_length\taligned_bases\n")
    len_fh.write("read_id\tcall\torganism\trole\tread_length\n")

    def flush(qname, per_org, rlen):
        """Resolve one read's alignments into a single call."""
        if not per_org:
            counts["unassigned"]["reads"] += 1
            counts["unassigned"]["read_bases"] += rlen
            assign_fh.write(f"{qname}\tunassigned\tunmapped\tnone\t0\t0\t0\t{rlen}\t0\n")
            len_fh.write(f"{qname}\tunmapped\tunassigned\tnone\t{rlen}\n")
            return

        ranked = sorted(per_org.items(), key=lambda kv: kv[1]["as"], reverse=True)

        if args.mode == "sequential":
            # Reproduce the classic subtraction chain -- "remove everything that
            # maps to the carrier, then everything that maps to the contaminant,
            # then map the survivors to the community" -- but evaluated against
            # the SAME alignments used by competitive mode. Deriving both rules
            # from one mapping pass means any difference between the two modes is
            # attributable to the decision rule alone, not to re-mapping the
            # survivors against a smaller index.
            #
            # The known consequence, and the reason competitive mode exists: a
            # read from the community's own E. coli that also aligns to the
            # carrier-derived E. coli K-12 is claimed by the subtraction step and
            # lost, even when it matches the community strain better.
            for role_to_subtract in [r.strip() for r in args.subtract_order.split(",")]:
                hits = [(o, v) for o, v in ranked
                        if org2role.get(o) == role_to_subtract
                        and v["as"] >= args.min_subtract_as]
                if hits:
                    org, v = max(hits, key=lambda kv: kv[1]["as"])
                    counts[org]["reads"] += 1
                    counts[org]["read_bases"] += rlen
                    counts[org]["aligned_bases"] += v["aligned"]
                    assign_fh.write(f"{qname}\t{org}\tsubtracted\t{role_to_subtract}\t"
                                    f"{v['as']}\t0\t0\t{rlen}\t{v['aligned']}\n")
                    len_fh.write(f"{qname}\tsubtracted\t{org}\t{role_to_subtract}\t{rlen}\n")
                    return

            remaining = [(o, v) for o, v in ranked
                         if org2role.get(o) not in
                         [r.strip() for r in args.subtract_order.split(",")]]
            if not remaining:
                counts["unassigned"]["reads"] += 1
                counts["unassigned"]["read_bases"] += rlen
                assign_fh.write(f"{qname}\tunassigned\tno_survivor\tnone\t0\t0\t0\t{rlen}\t0\n")
                len_fh.write(f"{qname}\tno_survivor\tunassigned\tnone\t{rlen}\n")
                return
            ranked = remaining

        best_org, best = ranked[0]
        runner_as = ranked[1][1]["as"] if len(ranked) > 1 else 0
        margin = best["as"] - runner_as
        threshold = max(args.min_margin_abs, int(args.min_margin_frac * best["as"]))

        if len(ranked) == 1:
            # Only one organism in the reference set aligned this read at all,
            # so there is no competition to resolve. Ambiguity is a statement
            # about two organisms being indistinguishable, not about a single
            # alignment being weak -- weak-but-unique alignments are still
            # unique, and calling them "ambiguous" would both mislabel them and
            # silently drain reads out of the per-organism counts.
            call, organism = "assigned", best_org
        elif margin >= threshold:
            call, organism = "assigned", best_org
        else:
            # Tied between two or more organisms: name them so the ambiguity is
            # auditable rather than hidden (e.g. "ambiguous:E. coli|E. coli K-12").
            tied = sorted(o for o, v in ranked if best["as"] - v["as"] < threshold)
            call, organism = "ambiguous", "ambiguous:" + "|".join(tied)

        role = org2role.get(best_org, "unknown") if call == "assigned" else "ambiguous"
        counts[organism]["reads"] += 1
        counts[organism]["read_bases"] += rlen
        counts[organism]["aligned_bases"] += best["aligned"]

        assign_fh.write(f"{qname}\t{organism}\t{call}\t{role}\t{best['as']}\t"
                        f"{runner_as}\t{margin}\t{rlen}\t{best['aligned']}\n")
        len_fh.write(f"{qname}\t{call}\t{organism}\t{role}\t{rlen}\n")

    mode = "rb" if args.bam.endswith(".bam") else "r"
    with pysam.AlignmentFile(args.bam, mode, check_sq=False) as bam:
        cur_q, per_org, cur_len = None, {}, 0
        unknown_contigs = set()

        for aln in bam.fetch(until_eof=True):
            qname = aln.query_name
            if qname != cur_q:
                if cur_q is not None:
                    flush(cur_q, per_org, cur_len)
                cur_q, per_org, cur_len = qname, {}, 0

            rlen = read_length(aln)
            if rlen > cur_len:
                cur_len = rlen

            if aln.is_unmapped:
                continue
            if aln.mapping_quality < args.min_mapq:
                continue

            org = contig2org.get(aln.reference_name)
            if org is None:
                unknown_contigs.add(aln.reference_name)
                continue

            score = alignment_score(aln)
            aligned = aln.query_alignment_length or 0
            prev = per_org.get(org)
            # Keep the single best-scoring alignment per organism. Supplementary
            # alignments of one read to the same organism are not summed: doing so
            # would let a fragmented reference outscore a contiguous one.
            if prev is None or score > prev["as"]:
                per_org[org] = {"as": score, "aligned": aligned}

        if cur_q is not None:
            flush(cur_q, per_org, cur_len)

    assign_fh.close()
    len_fh.close()

    if unknown_contigs:
        print(f"warning: {len(unknown_contigs)} contig(s) absent from the contig map, "
              f"e.g. {sorted(unknown_contigs)[:3]}", file=sys.stderr)

    with open(f"{args.prefix}.counts.tsv", "w") as fh:
        fh.write("organism\trole\treads\tread_bases\taligned_bases\n")
        for organism in sorted(counts):
            c = counts[organism]
            if organism == "unassigned":
                role = "none"
            elif organism.startswith("ambiguous:"):
                role = "ambiguous"
            else:
                role = org2role.get(organism, "unknown")
            fh.write(f"{organism}\t{role}\t{c['reads']}\t{c['read_bases']}\t{c['aligned_bases']}\n")

    total = sum(c["reads"] for c in counts.values())
    print(f"[assign_reads] {total:,} reads across {len(counts)} classes -> {args.prefix}.counts.tsv")


if __name__ == "__main__":
    main()
