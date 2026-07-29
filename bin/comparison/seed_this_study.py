#!/usr/bin/env python3
"""
Regenerate `this_study.tsv` from pipeline output.

`this_study.tsv` is the seed the comparison figure falls back on when
`--results-dir` is not supplied. That fallback exists because the reads are
large: someone cloning the repository can redraw the published figure without
first acquiring ~100 GB of FASTQ and running the pipeline.

A fallback is only worth having if it agrees with the thing it stands in for.
The previous seed did not. It was extracted from a hand-maintained spreadsheet
before the carrier-derived *E. coli* correction, so running the documented
command without `--results-dir` silently reproduced the figure the paper
explicitly corrects:

    seeded      Round 1  0.4247 reads/fg    Round 2  0.0325 reads/fg (n=1)
    published   Round 1  0.1427 reads/fg    Round 2  0.1730 reads/fg (n=3)

It also carried a third copy of the input masses -- `dna_pg` of 1000 for
lowinput_s1_r2, the nominal 1 ng, against the 1.21 ng recorded in
assets/measurements.tsv -- and cited cells of a workbook that is not in the
repository and never will be.

So the seed is now a *snapshot* rather than a parallel record: this script asks
comparison_data for exactly the rows `--results-dir` would produce and writes
them out. The two paths cannot disagree, because they are the same code. Re-run
it whenever the pipeline output changes, and commit the result.

Usage:
    python3 bin/comparison/seed_this_study.py --results-dir results
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comparison_data as cd  # noqa: E402

COLUMNS = ["study", "study_short", "condition", "organism", "replicate_idx",
           "replicate_label", "row_type", "classifier", "reads", "bases",
           "dna_pg", "dna_basis", "reads_per_fg", "bases_per_fg", "verified",
           "source", "source_detail", "provenance_note", "round", "experiment",
           "sample_id", "mode"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", type=Path, default=cd.DEFAULT_RESULTS_DIR)
    p.add_argument("--out", type=Path, default=cd.THIS_STUDY_TSV)
    p.add_argument("--mode", default="competitive")
    return p.parse_args()


def main():
    args = parse_args()
    if not args.results_dir.exists():
        sys.exit(f"error: {args.results_dir} does not exist; run the pipeline first")

    live = cd._live_rows(args.results_dir, mode=args.mode)
    if live.empty:
        sys.exit(f"error: no headline rows under {args.results_dir}")

    # The provenance of a committed seed is not the provenance of a live row:
    # the path it was read from is local to whoever generated it and means
    # nothing in a clone. Record what the value IS and how to reproduce it.
    live = live.copy()
    live["source"] = "pipeline_run_snapshot"
    live["source_detail"] = live["sample_id"].map(
        lambda s: f"results/{s}/{args.mode}/{s}.metrics.tsv, "
                  f"'{cd.HEADLINE_ORGANISM}' row")
    live["provenance_note"] = (
        "Snapshot of pipeline output, written by bin/comparison/seed_this_study.py "
        "so the figure can be redrawn without the reads. Identical to what "
        "--results-dir produces; input masses come from assets/measurements.tsv "
        "via the pipeline, not from this file.")

    out = live.reindex(columns=COLUMNS)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False,
               float_format=lambda v: repr(float(v)))

    print(f"[seed] {len(out)} row(s) -> {args.out}")
    for rnd, grp in out.groupby("round"):
        print(f"[seed]   {rnd}: n={len(grp)}, "
              f"mean {grp['reads_per_fg'].astype(float).mean():.6g} reads/fg, "
              f"{grp['bases_per_fg'].astype(float).mean():.6g} bases/fg")


if __name__ == "__main__":
    main()
