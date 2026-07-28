#!/usr/bin/env python3
"""
Pool per-base depth across the replicates of an experiment, per organism.

Per-replicate coverage answers "was this organism evenly covered in this
library". Pooled coverage answers a different and, for a log-distributed
community, more useful question: "with every read this study collected, what can
be said about this organism at all". Replicates are independent libraries of the
same material, so summing per-base depth is the same measurement a single deeper
run would have produced, and it is the only way the rare members reach a depth at
which uniformity means anything.

Pooling is per EXPERIMENT, never across both: lowinput_s1 and lowinput_s2 are
different communities and share no sample organism.

    ATTRIBUTION -- read this before using mean_depth
    ------------------------------------------------
    The depth files this reads are `samtools depth` over PRIMARY ALIGNMENTS to
    community contigs. They are not filtered to the reads that competitive
    assignment awarded to the organism, so for an organism with an abundant
    relative in the reference they measure mostly other organisms' reads. In
    lowinput_s1_r1, Enterococcus faecalis is awarded 4 reads (1,483 bases) yet
    shows 0.88x alignment depth, because 1,942 reads tie between it and
    Listeria monocytogenes and roughly half of those place their primary
    alignment on the E. faecalis genome.

    So this script reports BOTH:

        mean_depth            alignment depth, as above
        assigned_depth        aligned bases of reads actually awarded to this
                              organism, divided by genome size
        attributable          assigned_depth / mean_depth

    and never presents mean_depth alone. Where `attributable` is near 1 the two
    agree and the coverage statistics describe the organism. Where it is small
    the profile is a picture of whatever else aligns there, and Gini/CV describe
    that instead. Nine of the thirteen community members in this study sit at
    0.97-1.03; four do not.

Outputs mirror coverage_summary.py so the same plotting code can consume either:
a pooled summary TSV and a pooled binned profile TSV.

Usage:
    python3 bin/pool_coverage.py --out-summary results/summary/pooled_coverage_summary.tsv \\
                                 --out-profile results/summary/pooled_coverage_profile.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", default="results")
    p.add_argument("--samplesheets", default="assets/samplesheets")
    p.add_argument("--mode", default="competitive")
    p.add_argument("--window", type=int, default=1000)
    p.add_argument("--exclude", nargs="*", default=[],
                   help="sample_ids to leave out of the pool")
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-profile", required=True)
    return p.parse_args()


def gini(d: np.ndarray) -> float:
    """Gini of per-base depth. Zero-depth positions are included deliberately --
    dropout is the unevenness this is meant to detect."""
    if d.size == 0:
        return float("nan")
    total = d.sum(dtype=np.float64)
    if total <= 0:
        return float("nan")
    x = np.sort(d.astype(np.float64))
    n = x.size
    idx = np.arange(1, n + 1, dtype=np.float64)
    return float((2.0 * (idx * x).sum()) / (n * total) - (n + 1.0) / n)


def read_samplesheets(d):
    """sample_id -> (experiment, reference_set stem)."""
    out = {}
    for p in sorted(Path(d).glob("*.csv")):
        if p.name in ("test.csv", "raghavendra_2023.csv", "all.csv"):
            continue
        rows = [l for l in p.read_text().splitlines()
                if l.strip() and not l.startswith("#")]
        for r in csv.DictReader(rows):
            if r.get("sample_id"):
                out[r["sample_id"]] = (r["experiment"], Path(r["reference_set"]).stem)
    return out


def load_contig_map(path):
    """contig -> organism, for role=sample contigs only."""
    m = {}
    for line in Path(path).read_text().splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 3 and f[2] == "sample":
            m[f[0]] = f[1]
    return m


def assigned_bases(results, sid, mode):
    """Aligned bases of reads competitive assignment actually awarded, per organism.

    Read from counts.tsv, which books each read to exactly one class. Ambiguous
    classes are NOT split across their tied organisms: a tie is a statement that
    the read cannot be awarded, and dividing it up would manufacture attribution
    that the assignment step declined to make.
    """
    f = Path(results) / sid / mode / f"{sid}.counts.tsv"
    if not f.is_file():
        return {}
    c = pd.read_csv(f, sep="\t")
    c = c[c["role"] == "sample"]
    return dict(zip(c["organism"], c["aligned_bases"]))


def main():
    args = parse_args()
    root = Path(__file__).resolve().parent.parent
    results = root / args.results
    sheets = read_samplesheets(root / args.samplesheets)

    by_experiment = defaultdict(list)
    for sid, (exp, refset) in sorted(sheets.items()):
        if sid in args.exclude:
            print(f"[pool] excluding {sid}", file=sys.stderr)
            continue
        if (results / sid / "coverage" / f"{sid}.depth.tsv.gz").is_file():
            by_experiment[exp].append((sid, refset))
    if not by_experiment:
        sys.exit("error: no replicate has a depth file; run the pipeline first")

    summary_rows, profile_rows = [], []

    for exp, members in sorted(by_experiment.items()):
        refsets = {rs for _, rs in members}
        if len(refsets) != 1:
            sys.exit(f"error: {exp} mixes reference sets {refsets}; cannot pool")
        contig2org = load_contig_map(results / "references" / refsets.pop() / "contig_map.tsv")

        print(f"[pool] {exp}: {len(members)} replicates "
              f"({', '.join(s for s, _ in members)})", file=sys.stderr)

        depth = {}            # contig -> summed per-base depth
        assigned = defaultdict(float)
        for sid, _ in members:
            for org, b in assigned_bases(results, sid, args.mode).items():
                assigned[org] += b

            path = results / sid / "coverage" / f"{sid}.depth.tsv.gz"
            print(f"[pool]   reading {path.name}", file=sys.stderr)
            reader = pd.read_csv(path, sep="\t", header=None,
                                 names=["contig", "pos", "depth"],
                                 dtype={"contig": str, "pos": np.int64, "depth": np.int32},
                                 chunksize=8_000_000, compression="gzip")
            for chunk in reader:
                for contig, sub in chunk.groupby("contig", sort=False):
                    if contig not in contig2org:
                        continue
                    pos = sub["pos"].to_numpy()
                    d = sub["depth"].to_numpy().astype(np.int64)
                    arr = depth.get(contig)
                    need = int(pos.max())
                    if arr is None:
                        arr = np.zeros(need, dtype=np.int64)
                        depth[contig] = arr
                    elif arr.size < need:
                        arr = np.concatenate([arr, np.zeros(need - arr.size, dtype=np.int64)])
                        depth[contig] = arr
                    np.add.at(arr, pos - 1, d)

        # Regroup contigs into organisms.
        per_org = defaultdict(list)
        for contig, arr in depth.items():
            per_org[contig2org[contig]].append((contig, arr))

        for org in sorted(per_org):
            pieces = sorted(per_org[org])
            d = np.concatenate([a for _, a in pieces])
            n = d.size
            mean = float(d.mean()) if n else float("nan")
            sd = float(d.std()) if n else float("nan")
            adep = assigned.get(org, 0.0) / n if n else float("nan")
            summary_rows.append({
                "experiment": exp, "n_replicates": len(members), "organism": org,
                "positions": n, "mean_depth": f"{mean:.6g}",
                "median_depth": f"{float(np.median(d)):.6g}",
                "assigned_depth": f"{adep:.6g}",
                "attributable": f"{(adep / mean) if mean else float('nan'):.6g}",
                "breadth_1x": f"{float((d >= 1).mean()):.6g}",
                "breadth_5x": f"{float((d >= 5).mean()):.6g}",
                "breadth_10x": f"{float((d >= 10).mean()):.6g}",
                "cv": f"{(sd / mean) if mean else float('nan'):.6g}",
                "gini": f"{gini(d):.6g}",
                "dropout_fraction": f"{float((d == 0).mean()):.6g}",
            })
            print(f"[pool]   {org:34s} mean {mean:10.3f}x  assigned {adep:9.3f}x  "
                  f"attributable {(adep / mean) if mean else float('nan'):6.3f}  "
                  f"Gini {gini(d):6.4f}", file=sys.stderr)

            offset = 0
            for contig, arr in pieces:
                nb = int(np.ceil(arr.size / args.window))
                pad = nb * args.window - arr.size
                binned = np.concatenate([arr, np.zeros(pad, dtype=np.int64)]).reshape(nb, -1)
                counts = np.full(nb, args.window, dtype=np.int64)
                counts[-1] -= pad
                means = binned.sum(axis=1) / counts
                for i, (m, c) in enumerate(zip(means, counts)):
                    profile_rows.append({
                        "experiment": exp, "organism": org, "contig": contig,
                        "bin_start": i * args.window,
                        "global_bin_start": offset + i * args.window,
                        "mean_depth": f"{m:.6g}", "n_positions": int(c),
                    })
                offset += arr.size

    for path, rows in ((args.out_summary, summary_rows), (args.out_profile, profile_rows)):
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"[pool] -> {path}  ({len(rows)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
