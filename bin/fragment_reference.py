#!/usr/bin/env python3
"""Cut a complete genome where a contig-level assembly of the same species breaks.

Used by `make fragmentation` (docs/reference-fragmentation.md) to ask whether a
fragmented reference changes read assignment: the cut genome has exactly the
sequence of the complete one, so any difference between the two is the effect
of the breaks alone.

Where to cut comes from aligning the assembly's contigs to the complete genome
(minimap2 -x asm20 -c --secondary=no; PAF on --paf). Walking along each
chromosome of the complete genome through the alignments of at least --min-aln
bp with MAPQ > 0, a cut goes wherever the aligned contig changes, midway
between the two alignments. Cuts within --edge bp of a chromosome end or of the
previous cut are dropped. Real assemblies break in repeats (rDNA, Ty elements,
subtelomeres), so the cuts fall where an assembly actually would.

Writes the cut genome (pieces named <chromosome>_frag<N>) and, with --cuts, the
cut positions as TSV. Prints the number of pieces and the N50.

Usage:
    fragment_reference.py --complete B.fna --paf A_on_B.paf --out C.fna [--cuts cuts.tsv]
"""
import argparse


def read_fasta(path):
    seqs, order, name = {}, [], None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                order.append(name)
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}, order


def cut_positions(paf, lengths, min_aln, edge):
    alns = {}
    with open(paf) as f:
        for line in f:
            p = line.split("\t")
            if int(p[11]) == 0 or int(p[10]) < min_aln:
                continue
            alns.setdefault(p[5], []).append((int(p[7]), int(p[8]), p[0]))
    cuts = {}
    for chrom, a in alns.items():
        a.sort()
        c = []
        for (s0, e0, q0), (s1, _e1, q1) in zip(a, a[1:]):
            if q1 != q0:
                x = (e0 + s1) // 2 if s1 > e0 else (s1 + e0) // 2
                if edge < x < lengths[chrom] - edge and (not c or x - c[-1] > edge):
                    c.append(x)
        cuts[chrom] = c
    return cuts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--complete", required=True, help="complete genome FASTA")
    ap.add_argument("--paf", required=True, help="assembly contigs aligned to it")
    ap.add_argument("--out", required=True, help="cut genome FASTA to write")
    ap.add_argument("--cuts", help="also write the cut positions here (TSV)")
    ap.add_argument("--min-aln", type=int, default=10_000)
    ap.add_argument("--edge", type=int, default=5_000)
    a = ap.parse_args()

    seqs, order = read_fasta(a.complete)
    cuts = cut_positions(a.paf, {k: len(v) for k, v in seqs.items()}, a.min_aln, a.edge)
    pieces = []
    with open(a.out, "w") as f:
        for chrom in order:
            s = seqs[chrom]
            pos = [0] + cuts.get(chrom, []) + [len(s)]
            for i, (p0, p1) in enumerate(zip(pos, pos[1:]), 1):
                f.write(f">{chrom}_frag{i}\n")
                for j in range(p0, p1, 80):
                    f.write(s[j:min(j + 80, p1)] + "\n")
                pieces.append(p1 - p0)
    if a.cuts:
        with open(a.cuts, "w") as f:
            f.write("chrom\tpos\n")
            for chrom in order:
                for x in cuts.get(chrom, []):
                    f.write(f"{chrom}\t{x}\n")
    pieces.sort(reverse=True)
    total, acc, n50 = sum(pieces), 0, 0
    for length in pieces:
        acc += length
        if acc >= total / 2:
            n50 = length
            break
    print(f"{a.out}: {len(order)} sequences -> {len(pieces)} pieces "
          f"({sum(len(v) for v in cuts.values())} cuts), N50 {n50:,} bp, total {total:,} bp")


if __name__ == "__main__":
    main()
