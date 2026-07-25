#!/usr/bin/env python3
"""
Quantify what sequential subtraction costs, per organism.

Sequential subtraction removes every read that aligns to the carrier, then every
read that aligns to the carrier-derived contaminant, and assigns whatever
survives. When a community member shares sequence with either of those --
as the lowinput_s1 community's own *E. coli* shares a core genome with the
lambda-prep *E. coli* K-12 -- its reads are removed along with the contaminant
and the loss leaves no trace in the output.

Both rules are evaluated over the *same* alignments, so every difference here is
attributable to the decision rule alone rather than to re-mapping.

The figure plots, per organism, the fraction of competitive-mode reads that
sequential mode retains. An organism with no homology to carrier or contaminant
sits at 1.0; anything below that is signal destroyed by subtraction.
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


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--per-organism", required=True,
                   help="results/summary/per_organism.tsv (must contain both modes)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--basename", default="mode_delta")
    return p.parse_args()


def italicize(name):
    if not name or " " not in name or name.startswith(("ambiguous", "All ", "unassigned")):
        return name
    return "$\\mathit{" + name.replace(" ", r"\ ") + "}$"


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.per_organism, sep="\t")
    modes = set(df["mode"].unique())
    if not {"competitive", "sequential"} <= modes:
        raise SystemExit(f"error: need both modes in {args.per_organism}; found {sorted(modes)}. "
                         "Re-run the pipeline with --mode both.")

    keep = df[df["role"].isin(["sample", "contaminant", "carrier"])]
    piv = (keep.pivot_table(index=["experiment", "organism", "role"],
                            columns="mode", values=["reads", "bases"],
                            aggfunc="sum")
             .reset_index())
    piv.columns = [c[0] if not c[1] else f"{c[0]}_{c[1]}" for c in piv.columns]

    piv["reads_retained"] = np.where(piv["reads_competitive"] > 0,
                                     piv["reads_sequential"] / piv["reads_competitive"],
                                     np.nan)
    piv["bases_retained"] = np.where(piv["bases_competitive"] > 0,
                                     piv["bases_sequential"] / piv["bases_competitive"],
                                     np.nan)
    piv["reads_lost"] = piv["reads_competitive"] - piv["reads_sequential"]
    piv = piv.sort_values(["experiment", "role", "reads_retained"])
    piv.to_csv(outdir / f"{args.basename}.csv", index=False)

    plot = piv[(piv["role"] == "sample") & piv["reads_retained"].notna()
               & (piv["reads_competitive"] > 0)].copy()
    plot = plot.sort_values("reads_retained")

    fig, ax = plt.subplots(figsize=(7.2, max(3.2, 0.30 * len(plot) + 1.4)))
    ypos = np.arange(len(plot))
    # Anything that loses more than a couple of percent is worth flagging.
    colors = ["#D55E00" if v < 0.98 else "#009E73" for v in plot["reads_retained"]]
    ax.barh(ypos, plot["reads_retained"], color=colors, edgecolor="black",
            linewidth=0.4, height=0.68, zorder=3)
    ax.axvline(1.0, color="0.35", lw=0.9, ls=(0, (4, 3)), zorder=4)

    labels = [f"{italicize(o)}  ({e.replace('lowinput_', '')})"
              for o, e in zip(plot["organism"], plot["experiment"])]
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Reads retained by sequential subtraction "
                  "(fraction of competitive assignment)")
    ax.set_xlim(0, 1.08)
    ax.set_title("What sequential subtraction removes", loc="left")

    for y, (_, r) in zip(ypos, plot.iterrows()):
        if r["reads_retained"] < 0.98:
            ax.text(r["reads_retained"] + 0.015, y,
                    f"{r['reads_retained']*100:.1f}%  "
                    f"({int(r['reads_competitive']):,} → {int(r['reads_sequential']):,})",
                    va="center", fontsize=7.5, color="#D55E00")
    ax.grid(True, axis="x", color="0.93", lw=0.5, zorder=0)
    fig.tight_layout()
    fig.savefig(outdir / f"{args.basename}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{args.basename}.png", bbox_inches="tight", dpi=600)
    plt.close(fig)

    worst = plot.nsmallest(1, "reads_retained").iloc[0] if len(plot) else None
    payload = {
        "id": args.basename,
        "title": "Community reads destroyed by sequential subtraction",
        "caption": (
            "Fraction of each community organism's reads that survive sequential "
            "carrier-then-contaminant subtraction, relative to competitive assignment "
            "over the identical alignments. Organisms sharing no sequence with the "
            "lambda carrier or the carrier-derived Escherichia coli K-12 are unaffected "
            "and sit at 1.0. The exception is the mock community's own E. coli, whose "
            "core genome is shared with the contaminant: subtraction cannot tell the two "
            "strains apart and removes both. Because subtraction leaves no record of "
            "what it discarded, this loss is invisible in a sequential analysis."),
        "source_files": [Path(args.per_organism).name],
        "software": ["python 3.12", f"pandas {pd.__version__}",
                     f"numpy {np.__version__}", f"matplotlib {matplotlib.__version__}"],
        "metrics": {
            "n_organisms": int(len(plot)),
            "n_materially_affected": int((plot["reads_retained"] < 0.98).sum()),
            "worst_case": ({
                "organism": worst["organism"],
                "experiment": worst["experiment"],
                "reads_competitive": int(worst["reads_competitive"]),
                "reads_sequential": int(worst["reads_sequential"]),
                "fraction_retained": float(worst["reads_retained"]),
            } if worst is not None else None),
            "per_organism": plot[["experiment", "organism", "reads_competitive",
                                  "reads_sequential", "reads_retained",
                                  "bases_retained"]].to_dict("records"),
        },
    }
    (outdir / f"{args.basename}.json").write_text(json.dumps(payload, indent=2))

    print(plot[["experiment", "organism", "reads_competitive", "reads_sequential",
                "reads_retained"]].to_string(index=False))
    print(f"\n[mode-delta] wrote {args.basename}.pdf/.png/.csv/.json in {outdir}")


if __name__ == "__main__":
    main()
