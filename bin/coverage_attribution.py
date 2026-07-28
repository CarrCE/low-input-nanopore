#!/usr/bin/env python3
"""
Per-replicate coverage depth, alignment vs attributable, and the 1x threshold.

`samtools depth` measures primary alignments to a community genome. That is not
the same as the depth of reads competitive assignment awarded to the organism,
and for a member sharing sequence with an abundant relative the two differ by
orders of magnitude (see docs/TODO.md item 3 and bin/pool_coverage.py).

Which depth the 1x interpretability threshold is applied to therefore decides
which organism-replicate pairs are reported as characterisable. It has to be the
attributable one: a pair whose alignment depth is 65x but whose awarded reads
cover 0.83x has not been sequenced deeply enough to say anything about its
uniformity, and reporting it implies otherwise.

    attributable_depth = aligned bases of reads awarded to the organism
                         / genome size

Ambiguous classes are not split across their tied organisms. A tie is a
statement that the read cannot be awarded; apportioning it would manufacture
attribution the assignment step declined to make, in exactly the direction that
flatters the result.

The smoke-test sample is excluded by default. It is a 40,000-read synthetic
subsample of one replicate, not an experiment, and counting its rows inflates
the denominator of every "N of M pairs" statement made from this table.

Usage:
    python3 bin/coverage_attribution.py --out results/summary/coverage_attribution.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd

THRESHOLD = 1.0
EXCLUDE_DEFAULT = ("test_s2",)
COLS = ["sample_id", "organism", "genome_size", "alignment_depth",
        "attributable_depth", "attributable_fraction", "breadth_1x", "cv", "gini",
        "above_1x_alignment", "above_1x_attributable"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", default="results")
    p.add_argument("--mode", default="competitive")
    p.add_argument("--exclude", nargs="*", default=list(EXCLUDE_DEFAULT),
                   help="sample_ids that are not experiments (default: %(default)s)")
    p.add_argument("--out", default=None)
    return p.parse_args()


def genome_sizes(results: Path) -> dict:
    sizes = {}
    for f in sorted(results.glob("references/*/genome_sizes.tsv")):
        for _, r in pd.read_csv(f, sep="\t").iterrows():
            sizes[r["organism"]] = int(r["genome_size"])
    return sizes


def main():
    args = parse_args()
    root = Path(__file__).resolve().parent.parent
    results = root / args.results

    sizes = genome_sizes(results)
    if not sizes:
        sys.exit(f"error: no genome_sizes.tsv under {results}/references")

    rows = []
    for d in sorted(results.iterdir()):
        sid = d.name
        if not d.is_dir() or sid in args.exclude:
            continue
        summ = d / "coverage" / f"{sid}.coverage_summary.tsv"
        counts = d / args.mode / f"{sid}.counts.tsv"
        if not (summ.is_file() and counts.is_file()):
            continue

        c = pd.read_csv(counts, sep="\t")
        awarded = dict(zip(c.loc[c["role"] == "sample", "organism"],
                           c.loc[c["role"] == "sample", "aligned_bases"]))

        for _, r in pd.read_csv(summ, sep="\t").iterrows():
            org = r["organism"]
            g = sizes.get(org)
            if not g:
                sys.exit(f"error: no genome size for {org!r}; reference sets and "
                         f"coverage summaries disagree")
            align = float(r["mean_depth"])
            att = awarded.get(org, 0.0) / g
            rows.append({
                "sample_id": sid, "organism": org, "genome_size": g,
                "alignment_depth": f"{align:.6g}",
                "attributable_depth": f"{att:.6g}",
                "attributable_fraction": f"{(att / align) if align else float('nan'):.6g}",
                "breadth_1x": r["breadth_1x"], "cv": r["cv"], "gini": r["gini"],
                "above_1x_alignment": int(align >= THRESHOLD),
                "above_1x_attributable": int(att >= THRESHOLD),
            })

    if not rows:
        sys.exit("error: no replicate has both a coverage summary and counts")

    df = pd.DataFrame(rows)
    n = len(df)
    n_align = int(df["above_1x_alignment"].sum())
    n_att = int(df["above_1x_attributable"].sum())
    kept = df[df["above_1x_attributable"] == 1]
    dropped = df[(df["above_1x_alignment"] == 1) & (df["above_1x_attributable"] == 0)]

    print(f"[attribution] {n} organism-replicate pairs "
          f"(excluding {', '.join(args.exclude) or 'nothing'})")
    print(f"[attribution]   >= {THRESHOLD:g}x on alignment depth    : {n_align}")
    print(f"[attribution]   >= {THRESHOLD:g}x on attributable depth : {n_att}"
          f"  ({kept['organism'].nunique()} organisms)")
    if len(dropped):
        print("[attribution] above threshold on alignment depth only "
              "-- not characterisable:")
        for _, r in dropped.sort_values("alignment_depth", key=lambda s: s.astype(float),
                                        ascending=False).iterrows():
            print(f"[attribution]   {r['sample_id']:16s} {r['organism']:26s} "
                  f"{float(r['alignment_depth']):8.1f}x -> "
                  f"{float(r['attributable_depth']):6.3f}x attributable")
    worst = kept["attributable_fraction"].astype(float)
    if len(worst):
        print(f"[attribution] of the {n_att} characterisable pairs, attributable/alignment "
              f"spans {worst.min():.3f}-{worst.max():.3f}")

    text = df[COLS].to_csv(sep="\t", index=False)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text)
        print(f"[attribution] -> {args.out}")
    else:
        print(text)


if __name__ == "__main__":
    main()
