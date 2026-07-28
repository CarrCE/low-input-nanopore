#!/usr/bin/env python3
"""
Characterise coverage uniformity across community members.

Two panels:

  A  Binned depth along each genome, normalised to that organism's own mean, for
     organisms sequenced deeply enough for the profile to mean anything. A
     uniformly covered genome sits flat at 1.0; excursions are the coverage
     artifacts this panel exists to expose.

  B  Gini coefficient against mean depth for every organism and replicate, with
     the shallow-coverage region shaded. This panel is the honesty check on
     panel A: below ~1x mean depth almost every position is zero simply because
     the genome was not sampled, so Gini approaches 1 for any organism
     regardless of how uniform it truly is. Points in the shaded region say
     nothing about uniformity.

Inputs are the per-sample coverage_summary.tsv and coverage_profile.tsv written
by coverage_summary.py.
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

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "dejavusans",
})

# Okabe-Ito
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#F0E442", "#000000"]

MIN_DEPTH_INTERPRETABLE = 1.0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summaries", nargs="+", required=True,
                   help="*.coverage_summary.tsv")
    p.add_argument("--profiles", nargs="+", required=True,
                   help="*.coverage_profile.tsv")
    p.add_argument("--attribution", default=None,
                   help="coverage_attribution.tsv; marks panels whose alignment "
                        "depth is largely not attributable to awarded reads")
    p.add_argument("--exclude", nargs="*", default=["test_s2"],
                   help="sample_ids that are not experiments. The smoke test is a "
                        "40,000-read synthetic subsample; the input globs match it, "
                        "and counting it inflates every 'N of M pairs' statement "
                        "made from this figure (default: %(default)s)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--basename", default="coverage")
    return p.parse_args()


def italicize(name):
    if not name or " " not in name:
        return name
    return "$\\mathit{" + name.replace(" ", r"\ ") + "}$"


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summ = pd.concat([pd.read_csv(f, sep="\t") for f in args.summaries],
                     ignore_index=True)
    prof = pd.concat([pd.read_csv(f, sep="\t") for f in args.profiles],
                     ignore_index=True)

    if args.exclude:
        summ = summ[~summ["sample_id"].isin(args.exclude)]
        prof = prof[~prof["sample_id"].isin(args.exclude)]

    # Attributable depth decides what may be *interpreted*; alignment depth still
    # decides what is *drawn*, because the panel that is mostly another
    # organism's reads is the one Section S4 exists to explain. Drawn and marked
    # beats quietly dropped.
    attribution = {}
    if args.attribution:
        att = pd.read_csv(args.attribution, sep="\t")
        attribution = {(r["sample_id"], r["organism"]): float(r["attributable_fraction"])
                       for _, r in att.iterrows()}

    summ = summ.sort_values(["organism", "sample_id"])
    summ.to_csv(outdir / f"{args.basename}_summary.csv", index=False)

    deep = summ[summ["mean_depth"] >= MIN_DEPTH_INTERPRETABLE]
    organisms = sorted(deep["organism"].unique())

    # Small multiples for the profiles: overlaying six genomes of different
    # lengths on one shared axis is unreadable, and raw 1 kb bins are mostly
    # sampling noise at these depths. Each organism gets its own panel, its own
    # genome length, and a rolling median.
    SMOOTH_BINS = 25
    ncol = min(3, max(1, len(organisms)))
    nrow_prof = int(np.ceil(len(organisms) / ncol)) if organisms else 1

    fig = plt.figure(figsize=(4.0 * ncol, 2.35 * nrow_prof + 4.2))
    gs = fig.add_gridspec(nrow_prof + 1, ncol, height_ratios=[1] * nrow_prof + [1.9],
                          hspace=0.65, wspace=0.28)

    # Accumulate exactly what gets drawn, so the accompanying CSV lets a reader
    # redraw the figure without rerunning the pipeline.
    plotted_rows = []

    colors = {o: PALETTE[i % len(PALETTE)] for i, o in enumerate(organisms)}
    for i, org in enumerate(organisms):
        axp = fig.add_subplot(gs[i // ncol, i % ncol])
        cand = deep[deep["organism"] == org]
        sid = cand.loc[cand["mean_depth"].idxmax(), "sample_id"]   # deepest replicate
        mean_depth = float(cand.loc[cand["mean_depth"].idxmax(), "mean_depth"])
        sub = prof[(prof["organism"] == org) & (prof["sample_id"] == sid)] \
                .sort_values(["contig", "bin_start"])
        if sub.empty or mean_depth <= 0:
            continue
        y = sub["mean_depth"].to_numpy() / mean_depth
        binw = sub["bin_start"].diff().median() or 1000
        x = np.arange(len(y)) * binw / 1e6
        smooth = pd.Series(y).rolling(SMOOTH_BINS, center=True, min_periods=1).median()
        for xi, yi, si_ in zip(x, y, smooth):
            plotted_rows.append({
                "panel": "A", "organism": org, "sample_id": sid,
                "mean_depth": mean_depth, "position_mb": xi,
                "depth_over_mean_raw": yi, "depth_over_mean_smoothed": si_,
                "gini": "",
            })
        axp.fill_between(x, 0, y, color=colors[org], alpha=0.22, lw=0, zorder=2)
        axp.plot(x, smooth, color=colors[org], lw=1.1, zorder=3)
        axp.axhline(1.0, color="0.4", lw=0.8, ls=(0, (4, 3)), zorder=4)
        axp.set_ylim(0, 2.5)
        axp.set_xlim(0, x.max() if len(x) else 1)
        frac = attribution.get((sid, org))
        weak = frac is not None and frac < 0.9
        title = f"{italicize(org)}\n{sid}, {mean_depth:.1f}$\\times$"
        if weak:
            title += f" ({mean_depth * frac:.2f}$\\times$ attributable)"
        axp.set_title(title, fontsize=8, loc="left",
                      color="#8a3800" if weak else "black")
        if weak:
            for s in axp.spines.values():
                s.set_color("#8a3800"); s.set_linewidth(1.1); s.set_linestyle((0, (3, 2)))
        axp.tick_params(labelsize=7.5)
        axp.grid(True, color="0.93", lw=0.5, zorder=0)
        if i % ncol == 0:
            axp.set_ylabel("depth / mean", fontsize=8)
        if i // ncol == nrow_prof - 1:
            axp.set_xlabel("Genome position (Mb)", fontsize=8)

    if not organisms:
        axp = fig.add_subplot(gs[0, :])
        axp.text(0.5, 0.5, "no organism reached 1x mean depth",
                 ha="center", va="center", transform=axp.transAxes, color="0.4")
        axp.set_axis_off()

    fig.text(0.005, 0.985,
             f"A  Depth uniformity along each genome (deepest replicate; "
             f"{int(SMOOTH_BINS)}-bin rolling median over shading of raw 1 kb bins)",
             fontsize=9.5, va="top", ha="left")

    # ---- Panel B: Gini vs mean depth --------------------------------------
    axB = fig.add_subplot(gs[nrow_prof, :])
    axB.axvspan(summ["mean_depth"][summ["mean_depth"] > 0].min() * 0.5
                if (summ["mean_depth"] > 0).any() else 1e-4,
                MIN_DEPTH_INTERPRETABLE, color="0.88", zorder=0)
    axB.text(0.02, 0.04,
             "shaded: below 1× mean depth,\nGini reflects sparse sampling,\nnot true unevenness",
             transform=axB.transAxes, fontsize=7.5, color="0.35", va="bottom")

    all_orgs = sorted(summ["organism"].unique())
    colorsB = {o: PALETTE[i % len(PALETTE)] for i, o in enumerate(all_orgs)}
    for org in all_orgs:
        sub = summ[summ["organism"] == org]
        axB.scatter(sub["mean_depth"], sub["gini"], s=34, alpha=0.9,
                    color=colorsB[org], edgecolors="black", linewidths=0.4,
                    label=italicize(org), zorder=3)
        for _, r in sub.iterrows():
            plotted_rows.append({
                "panel": "B", "organism": org, "sample_id": r["sample_id"],
                "mean_depth": r["mean_depth"], "position_mb": "",
                "depth_over_mean_raw": "", "depth_over_mean_smoothed": "",
                "gini": r["gini"],
            })
    axB.set_xscale("log")
    axB.set_xlabel("Mean depth (×)")
    axB.set_ylabel("Gini coefficient of per-base depth")
    axB.set_ylim(0, 1.02)
    axB.axvline(MIN_DEPTH_INTERPRETABLE, color="0.4", lw=0.8, ls=(0, (4, 3)), zorder=2)
    axB.legend(fontsize=6.5, frameon=True, framealpha=0.95, edgecolor="0.7",
               loc="upper right", ncol=1)
    axB.set_title("B  Uniformity vs sequencing depth", loc="left")
    axB.grid(True, which="major", color="0.93", lw=0.5, zorder=0)

    # tight_layout fights the explicit gridspec spacing set above; the figure is
    # already laid out, so only trim the surrounding whitespace on save.
    fig.savefig(outdir / f"{args.basename}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{args.basename}.png", bbox_inches="tight", dpi=600)
    plt.close(fig)

    pd.DataFrame(plotted_rows).to_csv(outdir / f"{args.basename}.csv", index=False)

    # "Interpretable" is now attributable-depth based. An organism-replicate pair
    # whose awarded reads cover less than 1x has not been sequenced deeply enough
    # to characterise, whatever its alignment depth says.
    if attribution:
        att_depth = deep.apply(
            lambda r: r["mean_depth"] * attribution.get((r["sample_id"], r["organism"]), 1.0),
            axis=1)
        interp_rows = deep[att_depth >= MIN_DEPTH_INTERPRETABLE]
    else:
        interp_rows = deep
    interpretable = interp_rows[["sample_id", "organism", "mean_depth", "breadth_1x",
                                 "cv", "gini"]].to_dict("records")
    payload = {
        "id": args.basename,
        "title": "Coverage uniformity across community members",
        "caption": (
            "(A) Depth in 1 kb bins along each reference genome, normalised to that "
            "organism's own mean depth, for organisms whose mean depth reached 1x; the "
            "deepest replicate is shown, with a 25-bin rolling median over the raw bins. "
            "A uniformly covered genome sits at 1.0. Depth is over primary alignments "
            "and is not restricted to reads competitive assignment awarded to the "
            "organism; panels outlined in dashed rule are those where the two differ, "
            "and carry the attributable depth in the panel title. The Escherichia coli "
            "panel is the extreme case at 1.3% attributable, and does NOT show a "
            "biological coverage artifact: its depth is concentrated at sequence "
            "B-1109 shares with the lambda carrier and with the K-12 contaminant, so it "
            "is a picture of the reference set rather than of the sample. It is drawn "
            "rather than dropped because it is the cautionary example the accompanying "
            "text explains, but it is excluded from the characterisable pairs on "
            "attributable depth. Listeria monocytogenes, by contrast, is flat "
            "at 162x apart from two discrete dropouts, which are genuine artifacts. "
            "(B) Gini "
            "coefficient of per-base depth against mean depth for every organism and "
            "replicate. The shaded region marks coverage below 1x, where most positions "
            "are zero because the genome was not sampled rather than because coverage is "
            "uneven; Gini there approaches 1 for any organism and carries no information "
            "about uniformity. In a log-distributed community only the most abundant "
            "members clear that threshold, which is itself the central limitation on "
            "assessing coverage artifacts at these input masses."),
        "source_files": sorted(Path(f).name for f in args.summaries + args.profiles),
        "software": ["python 3.12", f"pandas {pd.__version__}",
                     f"numpy {np.__version__}", f"matplotlib {matplotlib.__version__}"],
        "metrics": {
            "n_plotted_points": int(len(plotted_rows)),
            "n_organism_replicate_pairs": int(len(summ)),
            "excluded_samples": list(args.exclude),
            "n_above_1x_alignment_depth": int(len(deep)),
            "n_above_1x_attributable_depth": int(len(interp_rows)),
            "organisms_above_1x_alignment": organisms,
            "organisms_characterisable": sorted(interp_rows["organism"].unique()),
            "interpretable_rows": interpretable,
        },
    }
    (outdir / f"{args.basename}.json").write_text(json.dumps(payload, indent=2))

    print(f"[coverage] {len(summ)} organism-replicate pairs "
          f"(excluded: {', '.join(args.exclude) or 'nothing'}); "
          f"{len(deep)} above {MIN_DEPTH_INTERPRETABLE:g}x alignment depth, "
          f"{len(interp_rows)} above {MIN_DEPTH_INTERPRETABLE:g}x attributable depth "
          f"({interp_rows['organism'].nunique()} organisms)")
    if organisms:
        print(f"[coverage] interpretable: {', '.join(organisms)}")
    print(f"[coverage] wrote {args.basename}.pdf/.png/.csv/.json in {outdir}")


if __name__ == "__main__":
    main()
