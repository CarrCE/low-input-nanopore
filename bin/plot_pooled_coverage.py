#!/usr/bin/env python3
"""
Coverage uniformity across EVERY community member, pooled across replicates.

The per-replicate coverage figure can only show the handful of organisms that
reach interpretable depth in a single library. Pooling the replicates of an
experiment -- summing per-base depth, which is what a single deeper run would
have produced -- lets every member be shown on one axis, and moves two more
organisms above 1x.

  A  Pooled depth along each genome, normalised to that organism's own mean,
     one panel per community member, ordered by depth. Reading down the grid is
     the figure's argument: flat at high depth, progressively spikier as depth
     falls, until at the bottom the profile is a handful of isolated reads.

  B  Gini against pooled depth for all thirteen members.

ATTRIBUTION. Depth here is of PRIMARY ALIGNMENTS, not of reads competitive
assignment awarded (see pool_coverage.py). For most members the two agree to
within 1%, but for three they do not, and for those the profile is a picture of
whichever abundant relative aligns there. Rather than drop them or show them
unmarked, both quantities are drawn: a filled marker at the depth actually
attributable to the organism, a hollow marker at the raw alignment depth, and a
connector between. A long connector is the point -- it says the alignment depth
is borrowed.

Colour separates the two experiments (different communities, no shared member).
Organism identity is carried by panel titles and direct labels rather than by a
thirteen-hue legend, which no reader could decode and no palette could keep
colourblind-safe.

Usage:
    python3 bin/plot_pooled_coverage.py \\
        --summary results/summary/pooled_coverage_summary.tsv \\
        --profile results/summary/pooled_coverage_profile.tsv \\
        --outdir results/summary
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "dejavusans",
})

# Two hues only, one per experiment. Validated for CVD separation
# (worst adjacent pair dE 21.9 protan, 31.2 normal) against a light surface.
EXP_COLOR = {"lowinput_s1": "#0072B2", "lowinput_s2": "#D55E00"}

MIN_DEPTH_INTERPRETABLE = 1.0
# Below this, alignment depth is mostly reads assigned elsewhere and the profile
# describes another organism. Chosen an order of magnitude clear of the ~1%
# measurement noise on the 1.00 cases.
MIN_ATTRIBUTABLE = 0.9
SMOOTH_BINS = 25
NCOL = 5


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", required=True,
                   help="pooled summary; normally the assignment-filtered one")
    p.add_argument("--alignment-summary", default=None,
                   help="the same pooling over raw alignment depth. When given, each "
                        "point also shows where it would have sat had depth not been "
                        "filtered to awarded reads -- which for three of the thirteen "
                        "members is a different number by orders of magnitude")
    p.add_argument("--profile", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--basename", default="pooled_coverage")
    return p.parse_args()


def italic(name):
    return "$\\mathit{" + name.replace(" ", r"\ ") + "}$" if " " in name else name


def abbrev(name):
    """'Listeria monocytogenes' -> 'L. monocytogenes', for direct labels."""
    parts = name.split(" ", 1)
    return f"{parts[0][0]}. {parts[1]}" if len(parts) == 2 else name


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summ = pd.read_csv(args.summary, sep="\t")
    prof = pd.read_csv(args.profile, sep="\t")

    # `attributable` inside the summary compares two routes to the same quantity
    # once depth is already assignment-filtered, so it is ~1 by construction and
    # says nothing. What a reader needs is how far the filtering moved each
    # point, which is this pooling against the alignment-depth pooling.
    if args.alignment_summary:
        al = pd.read_csv(args.alignment_summary, sep="\t")
        amap = dict(zip(al["organism"], al["mean_depth"]))
        summ["alignment_depth"] = summ["organism"].map(amap)
    else:
        summ["alignment_depth"] = summ["mean_depth"]
    summ["attributable"] = summ["mean_depth"] / summ["alignment_depth"]

    # Ordered by depth so the grid reads as a depth series, not an alphabet.
    summ = summ.sort_values("mean_depth", ascending=False).reset_index(drop=True)
    organisms = summ["organism"].tolist()
    nrow = int(np.ceil(len(organisms) / NCOL))

    fig = plt.figure(figsize=(2.65 * NCOL, 1.8 * nrow + 5.0))
    gs = fig.add_gridspec(nrow + 1, NCOL, height_ratios=[1] * nrow + [2.9],
                          hspace=0.95, wspace=0.32, top=0.95)

    plotted = []

    # ---- Panel A: one profile per member ----------------------------------
    for i, org in enumerate(organisms):
        ax = fig.add_subplot(gs[i // NCOL, i % NCOL])
        row = summ[summ["organism"] == org].iloc[0]
        color = EXP_COLOR.get(row["experiment"], "0.4")
        weak = row["attributable"] < MIN_ATTRIBUTABLE

        sub = prof[prof["organism"] == org].sort_values(["contig", "bin_start"])
        mean = float(row["mean_depth"])
        if sub.empty or mean <= 0:
            ax.set_axis_off()
            continue
        y = sub["mean_depth"].to_numpy() / mean
        x = np.arange(len(y)) * 1000 / 1e6
        smooth = pd.Series(y).rolling(SMOOTH_BINS, center=True, min_periods=1).median()

        # A genome covered at 0.03x has its reads in features far narrower than a
        # pixel; a translucent fill alone renders them invisible and the panel
        # reads as broken rather than as empty. Draw the raw series as a hairline
        # too, so isolated reads survive rasterisation.
        ax.fill_between(x, 0, np.clip(y, 0, 2.5), color=color, alpha=0.30, lw=0, zorder=2)
        ax.plot(x, np.clip(y, 0, 2.5), color=color, lw=0.25, alpha=0.75, zorder=3)
        ax.plot(x, smooth, color=color, lw=1.0, zorder=4)
        ax.axhline(1.0, color="0.45", lw=0.7, ls=(0, (4, 3)), zorder=4)
        ax.set_ylim(0, 2.5)
        ax.set_xlim(0, x.max() if len(x) else 1)
        ax.tick_params(labelsize=7)
        ax.grid(True, color="0.93", lw=0.5, zorder=0)

        # Depth is the thing that decides whether the panel means anything, so
        # it goes in the title next to the name.
        note = f"{mean:,.2f}$\\times$" if mean < 10 else f"{mean:,.0f}$\\times$"
        if weak:
            note += (f" of {row['alignment_depth']:,.0f}$\\times$ aligned")
        ax.set_title(f"{italic(org)}\n{note}", fontsize=7.4, loc="left",
                     color="0.25" if not weak else "#8a3800")
        if weak:
            for s in ax.spines.values():
                s.set_color("#8a3800"); s.set_linewidth(1.1); s.set_linestyle((0, (3, 2)))
        if i % NCOL == 0:
            ax.set_ylabel("depth / mean", fontsize=7.6)
        ax.set_xlabel("Mb", fontsize=7.2, labelpad=1)

        for xi, yi, si in zip(x, y, smooth):
            plotted.append({"panel": "A", "experiment": row["experiment"],
                            "organism": org, "position_mb": xi,
                            "depth_over_mean_raw": yi,
                            "depth_over_mean_smoothed": si,
                            "mean_depth": mean, "assigned_depth": row["assigned_depth"],
                            "attributable": row["attributable"], "gini": ""})

    for j in range(len(organisms), nrow * NCOL):          # blank out the tail
        fig.add_subplot(gs[j // NCOL, j % NCOL]).set_axis_off()

    fig.text(0.005, 0.995,
             "A  Pooled depth along each genome, every community member, ordered by depth "
             f"({SMOOTH_BINS}-bin rolling median over shaded raw 1 kb bins; clipped at 2.5)",
             fontsize=9.5, va="top", ha="left")

    # ---- Panel B: Gini vs depth -------------------------------------------
    axB = fig.add_subplot(gs[nrow, :])
    lo = float(summ[["mean_depth", "alignment_depth"]].values.min()) * 0.4
    axB.axvspan(lo, MIN_DEPTH_INTERPRETABLE, color="0.90", zorder=0)
    axB.axvline(MIN_DEPTH_INTERPRETABLE, color="0.4", lw=0.8, ls=(0, (4, 3)), zorder=2)

    labels = []
    for _, r in summ.iterrows():
        c = EXP_COLOR.get(r["experiment"], "0.4")
        weak = r["attributable"] < MIN_ATTRIBUTABLE
        if weak:
            # The gap between what aligns and what is attributable, drawn.
            axB.plot([r["mean_depth"], r["alignment_depth"]], [r["gini"]] * 2,
                     color=c, lw=1.0, ls=":", zorder=3, alpha=0.85)
            axB.scatter(r["alignment_depth"], r["gini"], s=44, facecolors="none",
                        edgecolors=c, linewidths=1.2, zorder=4)
        axB.scatter(r["mean_depth"], r["gini"],
                    s=40, color=c, edgecolors="black", linewidths=0.4, zorder=5)
        labels.append((r["mean_depth"],
                       float(r["gini"]), abbrev(r["organism"])))
        plotted.append({"panel": "B", "experiment": r["experiment"],
                        "organism": r["organism"], "position_mb": "",
                        "depth_over_mean_raw": "", "depth_over_mean_smoothed": "",
                        "mean_depth": r["mean_depth"], "assigned_depth": r["assigned_depth"],
                        "attributable": r["attributable"], "gini": r["gini"]})

    axB.set_xscale("log")
    axB.set_xlim(lo, float(summ["alignment_depth"].max()) * 2.6)
    axB.set_ylim(-0.03, 1.16)

    # Six of the thirteen sit within 0.03 Gini of each other along the top of the
    # panel, so a fixed label offset overplots them. Alternate above/below in
    # order of depth: with the points effectively collinear there, alternating is
    # enough to separate every pair without a layout solver.
    for k, (lx, ly, txt) in enumerate(sorted(labels)):
        above = (k % 2 == 0)
        axB.annotate(txt, (lx, ly), textcoords="offset points",
                     xytext=(0, 8 if above else -14), ha="center",
                     fontsize=6.4, color="0.25", style="italic", zorder=6)
    axB.set_xlabel("Pooled depth (×)")
    axB.set_ylabel("Gini coefficient of per-base depth")
    axB.grid(True, which="major", color="0.93", lw=0.5, zorder=0)
    axB.text(0.006, 0.05, "shaded: below 1× pooled depth,\nGini reflects sparse sampling,"
             "\nnot true unevenness", transform=axB.transAxes, fontsize=7,
             color="0.35", va="bottom")
    axB.set_title("B  Uniformity vs pooled sequencing depth", loc="left")

    handles = [Line2D([], [], marker="o", ls="none", color=EXP_COLOR["lowinput_s1"],
                      markeredgecolor="black", markeredgewidth=0.4, markersize=7,
                      label="lowinput_s1 (D6311, 3 replicates)"),
               Line2D([], [], marker="o", ls="none", color=EXP_COLOR["lowinput_s2"],
                      markeredgecolor="black", markeredgewidth=0.4, markersize=7,
                      label="lowinput_s2 (D6321, 4 replicates)"),
               Line2D([], [], marker="o", ls=":", color="0.35", markerfacecolor="none",
                      markeredgecolor="0.35", markersize=7,
                      label="filled: depth of awarded reads\nhollow: raw alignment depth")]
    # Parked mid-right: the only region of this panel with no marks in it, and
    # clear of the shaded-region note at lower left.
    axB.legend(handles=handles, fontsize=6.8, frameon=True, framealpha=0.95,
               edgecolor="0.7", loc="center right", bbox_to_anchor=(0.995, 0.60),
               handlelength=2.2, borderpad=0.5)

    fig.savefig(outdir / f"{args.basename}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{args.basename}.png", bbox_inches="tight", dpi=600)
    plt.close(fig)

    pd.DataFrame(plotted).to_csv(outdir / f"{args.basename}.csv", index=False)

    weak = summ[summ["attributable"] < MIN_ATTRIBUTABLE]
    payload = {
        "id": args.basename,
        "title": "Coverage uniformity across every community member, pooled across replicates",
        "caption": (
            "(A) Pooled per-base depth along each reference genome, normalised to that "
            "organism's own pooled mean, for every community member of both experiments, "
            "ordered by depth; 25-bin rolling median over shaded raw 1 kb bins, clipped at "
            "2.5x the mean. Replicates are pooled by summing per-base depth within an "
            "experiment, which is the measurement a single deeper run would have produced. "
            "(B) Gini coefficient of pooled per-base depth against pooled depth. Below 1x "
            "most positions are zero because the genome was not sampled, so Gini there "
            "approaches 1 regardless of true uniformity. Depth is measured over primary "
            "alignments, which for an organism with an abundant relative in the reference "
            "is not the same as the depth of reads competitive assignment awarded to it; "
            "filled markers give the attributable depth, hollow markers the raw alignment "
            "depth, and the dotted connector the difference. Organisms whose panels are "
            "outlined in the same way carry alignment depth that is mostly other "
            "organisms' reads."),
        "source_files": [Path(args.summary).name, Path(args.profile).name],
        "software": ["python 3.12", f"pandas {pd.__version__}",
                     f"numpy {np.__version__}", f"matplotlib {matplotlib.__version__}"],
        "metrics": {
            "n_members": int(len(summ)),
            "n_above_1x_alignment_depth": int((summ["mean_depth"] >= 1).sum()),
            "n_above_1x_attributable_depth": int((summ["assigned_depth"] >= 1).sum()),
            "low_attribution": weak[["organism", "mean_depth", "assigned_depth",
                                     "attributable"]].to_dict("records"),
            "rows": summ.to_dict("records"),
        },
    }
    (outdir / f"{args.basename}.json").write_text(json.dumps(payload, indent=2))

    print(f"[pooled-coverage] {len(summ)} members; "
          f"{int((summ['mean_depth'] >= 1).sum())} above 1x by alignment depth, "
          f"{int((summ['assigned_depth'] >= 1).sum())} by attributable depth")
    if len(weak):
        print("[pooled-coverage] low attribution: "
              + ", ".join(f"{r.organism} ({r.attributable:.3f})"
                          for r in weak.itertuples()))
    print(f"[pooled-coverage] wrote {args.basename}.pdf/.png/.csv/.json in {outdir}")


if __name__ == "__main__":
    main()
