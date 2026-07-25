#!/usr/bin/env python3
"""
Locate and annotate low-coverage regions in a per-base depth profile.

The coverage-artifact analysis reports Gini and coefficient of variation per
organism, which say *whether* coverage is uneven but not *where* or *why*. This
script answers those two: it finds contiguous runs of depth below a fraction of
the organism's median, merges runs separated by less than a gap tolerance, and --
given a GFF3 -- lists the annotated features each run overlaps.

Why a fraction of the median rather than an absolute floor: the organisms in
these datasets span two orders of magnitude in depth, so any fixed threshold
either misses depressions in the deep ones or drowns in noise on the shallow
ones. Poisson noise at 160x has a standard deviation of ~13x, so a run that
stays below 50% of the median for hundreds of bases is far outside sampling
variation and is structural.

Merging matters as much as thresholding. A broad depression is not flat -- it
dips, recovers slightly, dips again -- so a naive run-finder reports one
depression as several "discrete" dropouts, and the count depends entirely on the
threshold. Merging across --max-gap collapses those back into the single feature
they are.

Outputs a TSV: sample, organism, contig, start, end, length, min/mean depth,
depth relative to median, and (with --gff) the overlapping features.

Usage:
    python3 bin/coverage_dropouts.py results/<s>/coverage/<s>.depth.tsv.gz \\
        --contig-map results/references/<set>/contig_map.tsv \\
        --sample-id  <s> \\
        --organism   'Listeria monocytogenes' \\
        --gff        genomic.gff \\
        --out-tsv    dropouts.tsv
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict

import numpy as np

PRODUCT_RE = re.compile(r"product=([^;]*)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("depth", help="per-base depth TSV(.gz) from COVERAGE_PROFILE")
    p.add_argument("--contig-map", required=True, help="contig_map.tsv")
    p.add_argument("--sample-id", required=True)
    p.add_argument("--organism", default=None,
                   help="restrict to one organism (default: every sample-role organism)")
    p.add_argument("--gff", default=None,
                   help="GFF3 for feature annotation; coordinates must match the depth contigs")
    p.add_argument("--frac-of-median", type=float, default=0.5,
                   help="a position is low when depth < this fraction of the organism median")
    p.add_argument("--min-length", type=int, default=200,
                   help="discard merged regions shorter than this")
    p.add_argument("--max-gap", type=int, default=2000,
                   help="merge low runs separated by at most this many bases")
    p.add_argument("--out-tsv", required=True)
    return p.parse_args()


def load_contig_map(path):
    contig2org, org2role = {}, {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                contig2org[f[0]] = f[1]
                org2role[f[1]] = f[2]
    return contig2org, org2role


def load_depth(path, wanted_contigs):
    opener = gzip.open if path.endswith(".gz") else open
    pos = defaultdict(list)
    dep = defaultdict(list)
    with opener(path, "rt") as fh:
        for line in fh:
            c, s, d = line.rstrip("\n").split("\t")
            if c in wanted_contigs:
                pos[c].append(int(s))
                dep[c].append(int(d))
    return {c: (np.array(pos[c]), np.array(dep[c])) for c in pos}


def find_regions(pos, dep, threshold, max_gap, min_length):
    """Contiguous low runs, merged across gaps, as (start, end) in reference coords."""
    idx = np.flatnonzero(dep < threshold)
    if idx.size == 0:
        return []
    brk = np.flatnonzero(np.diff(idx) > 1)
    starts = np.r_[idx[0], idx[brk + 1]]
    ends = np.r_[idx[brk], idx[-1]]

    merged = []
    for s, e in zip(starts, ends):
        if merged and pos[s] - merged[-1][1] <= max_gap:
            merged[-1][1] = pos[e]
            merged[-1][3] = e
        else:
            merged.append([pos[s], pos[e], s, e])
    return [(a, b, si, ei) for a, b, si, ei in merged if b - a + 1 >= min_length]


def load_gff(path):
    feats = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("CDS", "rRNA", "tRNA"):
                continue
            m = PRODUCT_RE.search(f[8])
            feats[f[0]].append((int(f[3]), int(f[4]), f[2],
                                m.group(1) if m else "unannotated"))
    for c in feats:
        feats[c].sort()
    return feats


def overlapping(feats, contig, start, end):
    return [f for f in feats.get(contig, []) if f[1] >= start and f[0] <= end]


def main():
    args = parse_args()
    contig2org, org2role = load_contig_map(args.contig_map)

    if args.organism:
        wanted = {c for c, o in contig2org.items() if o == args.organism}
    else:
        wanted = {c for c, o in contig2org.items() if org2role.get(o) == "sample"}
    if not wanted:
        sys.exit(f"error: no contigs matched "
                 f"{'organism ' + args.organism if args.organism else 'role=sample'}")

    depth = load_depth(args.depth, wanted)
    if not depth:
        sys.exit(f"error: {args.depth} holds no rows for the requested contigs")

    feats = load_gff(args.gff) if args.gff else {}

    rows = []
    for contig, (pos, dep) in sorted(depth.items()):
        organism = contig2org[contig]
        median = float(np.median(dep))
        if median <= 0:
            print(f"skipping {contig}: median depth is 0", file=sys.stderr)
            continue
        threshold = args.frac_of_median * median

        for start, end, si, ei in find_regions(pos, dep, threshold,
                                               args.max_gap, args.min_length):
            window = dep[si:ei + 1]
            products = [f"{p} ({t})" if t != "CDS" else p
                        for _, _, t, p in overlapping(feats, contig, start, end)]
            rows.append({
                "sample_id": args.sample_id,
                "organism": organism,
                "contig": contig,
                "start": start,
                "end": end,
                "length": end - start + 1,
                "median_depth": f"{median:.1f}",
                "min_depth": int(window.min()),
                "mean_depth": f"{window.mean():.1f}",
                "frac_of_median": f"{window.mean() / median:.3f}",
                "n_features": len(products),
                "features": "; ".join(products) if products else "",
            })

    cols = ["sample_id", "organism", "contig", "start", "end", "length",
            "median_depth", "min_depth", "mean_depth", "frac_of_median",
            "n_features", "features"]
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"[coverage_dropouts] {len(rows)} regions below "
          f"{args.frac_of_median:.0%} of median -> {args.out_tsv}")
    for r in rows:
        print(f"    {r['organism']} {r['start']:,}-{r['end']:,} "
              f"({r['length']:,} bp, {r['frac_of_median']}x median, "
              f"{r['n_features']} features)")


if __name__ == "__main__":
    main()
