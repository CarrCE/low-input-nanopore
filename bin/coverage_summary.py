#!/usr/bin/env python3
"""
Characterise coverage uniformity per community organism.

Uneven coverage across mock-community members was one of the motivating
observations for this work, so this step quantifies it rather than leaving it as
an impression from a genome browser.

For each organism it reports:

    mean_depth, median_depth
    breadth_1x / _5x / _10x   fraction of the genome at >= that depth
    cv                        SD/mean of per-base depth -- scale-free spread
    gini                      inequality of depth across positions; 0 = perfectly
                              even, ->1 = all reads piled in one place
    dropout_fraction          fraction of positions at zero depth

CV and Gini are both reported because they fail differently: CV is inflated by a
few very deep positions, while Gini responds to the whole distribution. An
organism that is evenly covered except for one amplified locus shows high CV and
modest Gini; one that is covered in patches shows both.

IMPORTANT -- do not read these as evenness when depth is low. Below roughly 1x
mean depth most positions are zero simply because the genome was not sampled,
not because coverage is non-uniform, and Gini then approaches 1 for any organism
regardless of its true behaviour. At 0.1x mean depth a perfectly uniform genome
still yields Gini ~0.9. These statistics only carry information about coverage
artifacts once mean depth is comfortably above 1x, which for a log-distributed
community means only the most abundant members qualify. `mean_depth` is reported
first in the output for exactly this reason: check it before interpreting
anything to its right.

Also emits a binned depth profile (default 1 kb windows) so coverage can be
plotted along each genome without shipping per-base depth.

Input is `samtools depth -a` output (contig, position, depth), gzipped, which is
streamed rather than loaded.
"""
from __future__ import annotations

import argparse
import gzip
from collections import defaultdict

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("depth", help="samtools depth -a output, gzipped")
    p.add_argument("--contig-map", required=True, help="contig<TAB>organism<TAB>role")
    p.add_argument("--sample-id", required=True)
    p.add_argument("--window", type=int, default=1000, help="bin size in bp")
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-profile", required=True)
    return p.parse_args()


def load_contig_map(path):
    m = {}
    with open(path) as fh:
        header = True
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if header:
                header = False
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                m[parts[0]] = (parts[1], parts[2])
    return m


def gini(depths: np.ndarray) -> float:
    """
    Gini coefficient of the per-base depth distribution.

    Zero-depth positions are included deliberately: dropout is exactly the kind
    of unevenness this is meant to detect, and excluding it would make a patchy
    genome look uniform.
    """
    if depths.size == 0:
        return float("nan")
    total = depths.sum(dtype=np.float64)
    if total <= 0:
        return float("nan")
    x = np.sort(depths.astype(np.float64))
    n = x.size
    idx = np.arange(1, n + 1, dtype=np.float64)
    return float((2.0 * (idx * x).sum()) / (n * total) - (n + 1.0) / n)


def main():
    args = parse_args()
    contig2org = load_contig_map(args.contig_map)

    per_org_depths = defaultdict(list)
    bins = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # org -> bin -> [sum, n]

    opener = gzip.open if args.depth.endswith(".gz") else open
    with opener(args.depth, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig, pos, depth = parts[0], int(parts[1]), int(parts[2])
            entry = contig2org.get(contig)
            if entry is None:
                continue
            organism, role = entry
            if role != "sample":
                continue
            per_org_depths[organism].append(depth)
            b = (pos - 1) // args.window
            slot = bins[organism][(contig, b)]
            slot[0] += depth
            slot[1] += 1

    with open(args.out_summary, "w") as fh:
        fh.write("sample_id\torganism\tpositions\tmean_depth\tmedian_depth\t"
                 "breadth_1x\tbreadth_5x\tbreadth_10x\tcv\tgini\tdropout_fraction\n")
        for organism in sorted(per_org_depths):
            d = np.asarray(per_org_depths[organism], dtype=np.int64)
            n = d.size
            mean = float(d.mean()) if n else float("nan")
            sd = float(d.std()) if n else float("nan")
            fh.write(
                f"{args.sample_id}\t{organism}\t{n}\t{mean:.6g}\t{float(np.median(d)):.6g}\t"
                f"{float((d >= 1).mean()):.6g}\t{float((d >= 5).mean()):.6g}\t"
                f"{float((d >= 10).mean()):.6g}\t"
                f"{(sd / mean) if mean else float('nan'):.6g}\t"
                f"{gini(d):.6g}\t{float((d == 0).mean()):.6g}\n")
            print(f"[coverage] {organism:34s} mean {mean:9.3f}x  "
                  f"breadth1x {float((d >= 1).mean()):6.4f}  "
                  f"CV {(sd / mean) if mean else float('nan'):7.3f}  "
                  f"Gini {gini(d):6.4f}")

    with open(args.out_profile, "w") as fh:
        fh.write("sample_id\torganism\tcontig\tbin_start\tmean_depth\tn_positions\n")
        for organism in sorted(bins):
            for (contig, b), (total, count) in sorted(bins[organism].items()):
                fh.write(f"{args.sample_id}\t{organism}\t{contig}\t{b * args.window}\t"
                         f"{total / count:.6g}\t{count}\n")

    if not per_org_depths:
        print("[coverage] warning: no sample-role positions found; "
              "check that the contig map matches the depth file")


if __name__ == "__main__":
    main()
