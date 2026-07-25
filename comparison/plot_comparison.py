#!/usr/bin/env python3
"""
Display item: low-input nanopore performance, this study vs prior work.

Log-log scatter of reads / fg DNA into library prep (x) against bases / fg DNA
into library prep (y). Every prior-study sample is one point, coloured and
marker-shaped by study (Okabe-Ito, colourblind-safe). This study's rounds are
large bold black markers, read from `this_study.tsv` or from live pipeline
output. Dashed circles are iso-improvement contours around the Round 1 mean:
the geometric-mean fold improvement, expressed as Euclidean distance in
log10-log10 space, at 10x / 100x / 1000x.

Inputs are the versioned TSVs in this directory, NOT the legacy .xlsx:
    prior_studies.tsv   one row per prior-study sample per classifier variant
    this_study.tsv      seeded Round 1 / Round 2 values (source=legacy_spreadsheet)
    --results-dir       optional; live results/<sample_id>/<mode>/*.metrics.tsv
                        headline rows supersede the seeded values

Outputs (per the repo's display-item convention), written to --outdir:
    <basename>.pdf      vector, editable text (Type 42 fonts)
    <basename>.png      600 dpi raster
    <basename>.csv      every plotted point, with provenance
    <basename>.json     sidecar: id / title / caption / source_files / software / metrics

Known data defects handled here (full discussion in README.md):
  (a) Zorzano et al. 2025 axis inconsistency -> --zorzano-classifier
  (b) Mojarro et al. 2019 unsourced counts   -> loud warning when plotted
  (c) Raghavendra published vs reanalysed    -> --raghavendra-classifier
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import platform
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import comparison_data as cd  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_OUTDIR = HERE / "figures"
DEFAULT_BASENAME = "low_input_comparison"
DISPLAY_ITEM_ID = "fig_low_input_comparison"

# ---------------------------------------------------------------------------
# Style - colourblind-safe Okabe-Ito (unchanged from the legacy figure script)
# ---------------------------------------------------------------------------
COLORS = {
    "Mojarro et al. 2019":        "#0072B2",   # blue
    "B. Raghavendra et al. 2023": "#D55E00",   # vermilion
    "Zorzano et al. 2025":        "#009E73",   # bluish green
}
MARKERS = {
    "Mojarro et al. 2019":        "o",
    "B. Raghavendra et al. 2023": "s",
    "Zorzano et al. 2025":        "^",
}
FALLBACK_COLOR, FALLBACK_MARKER = "#CC79A7", "v"

THIS_STUDY_COLOR = "#000000"
ROUND_STYLE = {
    "Round 1": {"marker": "*", "s": 260, "legend_ms": 14, "edge_lw": 0.9},
    "Round 2": {"marker": "D", "s": 130, "legend_ms": 8,  "edge_lw": 0.8},
}
FALLBACK_ROUND_STYLE = {"marker": "P", "s": 150, "legend_ms": 9, "edge_lw": 0.8}

# Per-study callout offsets in log10 units, hand-tuned to avoid overlaps.
CALLOUT_OFFSET = {
    "Mojarro et al. 2019":        (+0.20, -0.75),
    "B. Raghavendra et al. 2023": (-0.10, -0.65),
    "Zorzano et al. 2025":        (-0.35, +1.25),
}

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":         9,
    "axes.titlesize":    10,
    "axes.labelsize":    10,
    "xtick.labelsize":   8.5,
    "ytick.labelsize":   8.5,
    "legend.fontsize":   8,
    "axes.linewidth":    0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype":      42,   # editable text in PDF
    "ps.fonttype":       42,
    "mathtext.fontset":  "dejavusans",
})

SPECIES_NAMES = (
    "Bacillus spizizenii", "Bacillus subtilis", "Escherichia coli",
    "Saccharomyces cerevisiae", "E. coli", "S. cerevisiae",
)


def italicize_species(text: str) -> str:
    """Wrap any species name in mathtext italic for matplotlib rendering."""
    for sp in SPECIES_NAMES:
        if sp in text:
            text = text.replace(sp, "$\\mathit{" + sp.replace(" ", r"\ ") + "}$")
    return text


def color_of(study):
    return COLORS.get(study, FALLBACK_COLOR)


def marker_of(study):
    return MARKERS.get(study, FALLBACK_MARKER)


def fmt_rate(v: float) -> str:
    """Readable fixed-point for a per-fg rate spanning many orders of magnitude."""
    if v >= 100:
        return f"{v:,.0f}"
    if v >= 1:
        return f"{v:,.2f}"
    if v >= 0.01:
        return f"{v:.4f}"
    return f"{v:.2e}"


# ---------------------------------------------------------------------------
# Figure pieces
# ---------------------------------------------------------------------------
def add_iso_improvement_contour(ax, cx, cy, log_radius, label,
                                color="0.4", lw=0.7, ls=(0, (3, 3))):
    """Circle of `log_radius` log10 units around (cx, cy), drawn on log-log axes."""
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(10 ** (np.log10(cx) + log_radius * np.cos(theta)),
            10 ** (np.log10(cy) + log_radius * np.sin(theta)),
            color=color, linewidth=lw, linestyle=ls, alpha=0.55, zorder=1)
    angle = np.deg2rad(225)   # label lower-left, where the prior studies live
    ax.text(10 ** (np.log10(cx) + log_radius * np.cos(angle)),
            10 ** (np.log10(cy) + log_radius * np.sin(angle)), label,
            color=color, fontsize=7.5, ha="center", va="center",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.8, alpha=0.9),
            zorder=2)


def pick_callouts(prior: pd.DataFrame, ref: dict) -> list[dict]:
    """Per study, the sample furthest below the reference point in log-log space."""
    callouts = []
    for study, grp in prior.groupby("study"):
        d = np.sqrt((np.log10(ref["reads"]) - np.log10(grp["reads_per_fg"])) ** 2 +
                    (np.log10(ref["bases"]) - np.log10(grp["bases_per_fg"])) ** 2)
        r = grp.loc[d.idxmax()]
        f_reads = ref["reads"] / r["reads_per_fg"]
        f_bases = ref["bases"] / r["bases_per_fg"]
        text = (f"{r['study_short']}\n{italicize_species(str(r['condition']))}\n"
                f"{f_reads:,.0f}x reads, {f_bases:,.0f}x bases")
        callouts.append({"x": float(r["reads_per_fg"]), "y": float(r["bases_per_fg"]),
                         "study": study, "condition": r["condition"],
                         "fold_reads": float(f_reads), "fold_bases": float(f_bases),
                         "text": text.replace("x reads", "× reads")
                                     .replace("x bases", "× bases")})
    return callouts


def build_figure(prior: pd.DataFrame, rounds: dict, title: str):
    fig, ax = plt.subplots(figsize=(7.0, 5.2))

    # ---- prior-study points ------------------------------------------------
    for study, grp in prior.groupby("study"):
        ax.scatter(grp["reads_per_fg"], grp["bases_per_fg"],
                   marker=marker_of(study), color=color_of(study),
                   s=42, alpha=0.85, edgecolors="black", linewidths=0.5,
                   label=study, zorder=3)

    # ---- this study, one series per round ----------------------------------
    for name, info in rounds.items():
        st = ROUND_STYLE.get(name, FALLBACK_ROUND_STYLE)
        ax.scatter(info["reads_values"], info["bases_values"],
                   marker=st["marker"], s=st["s"], c=THIS_STUDY_COLOR,
                   edgecolors="white", linewidths=st["edge_lw"], zorder=6)

    # ---- iso-improvement contours around the reference round ---------------
    ref_name = next(iter(rounds))
    ref = {"reads": rounds[ref_name]["reads_mean"], "bases": rounds[ref_name]["bases_mean"]}
    for log_r, label in [(1, "10× less"), (2, "100× less"), (3, "1000× less")]:
        add_iso_improvement_contour(ax, ref["reads"], ref["bases"], log_r, label)

    # ---- axis limits -------------------------------------------------------
    all_x = list(prior["reads_per_fg"]) + [v for r in rounds.values() for v in r["reads_values"]]
    all_y = list(prior["bases_per_fg"]) + [v for r in rounds.values() for v in r["bases_values"]]
    pad, pad_y_lo = 0.45, 0.95   # extra room at the bottom for callouts
    ax.set_xlim(10 ** (np.log10(min(all_x)) - pad), 10 ** (np.log10(max(all_x)) + pad))
    ax.set_ylim(10 ** (np.log10(min(all_y)) - pad_y_lo), 10 ** (np.log10(max(all_y)) + pad))
    ax.set_xscale("log")
    ax.set_yscale("log")

    # ---- constant bases-per-read reference ---------------------------------
    bpr_median = float(np.median(prior["bases_per_fg"] / prior["reads_per_fg"]))
    xs = np.array(ax.get_xlim())
    ax.plot(xs, bpr_median * xs, color="0.7", linewidth=0.7,
            linestyle=(0, (1, 2)), zorder=0)
    ax.text(xs[1] * 0.92, bpr_median * xs[1] * 0.55, f"~{bpr_median:,.0f} bp/read",
            color="0.5", fontsize=7.5, ha="right", va="center")

    # ---- per-study callouts ------------------------------------------------
    callouts = pick_callouts(prior, ref)
    for c in callouts:
        dx, dy = CALLOUT_OFFSET.get(c["study"], (0.20, 0.45))
        ax.annotate(
            c["text"], xy=(c["x"], c["y"]),
            xytext=(10 ** (np.log10(c["x"]) + dx), 10 ** (np.log10(c["y"]) + dy)),
            fontsize=7.5, ha="left", va="bottom", color=color_of(c["study"]),
            bbox=dict(facecolor="white", edgecolor=color_of(c["study"]),
                      linewidth=0.6, pad=2.5, alpha=0.92),
            arrowprops=dict(arrowstyle="-", color=color_of(c["study"]),
                            lw=0.6, alpha=0.7, connectionstyle="arc3,rad=0.12"),
            zorder=5)

    # ---- this-study labels -------------------------------------------------
    # Both rounds now sit in the same corner: correcting the Round 1 analysis
    # (which had been counting carrier-derived E. coli as community) moved it
    # down into Round 2's range, so the two clusters interleave. Push both
    # labels left and down, away from the points and away from the title -- the
    # old "Round 2 up and right" offset ran the box straight through the title.
    label_offsets = {"Round 1": (-0.55, -0.55, "right", "top"),
                     "Round 2": (-0.90, +0.18, "right", "center")}
    for name, info in rounds.items():
        dx, dy, ha, va = label_offsets.get(name, (+0.25, +0.40, "left", "bottom"))
        n = len(info["reads_values"])
        if n > 1:
            # anchor on the median replicate so the leader lands on a real marker
            order = np.argsort(info["reads_values"])
            j = int(order[n // 2])
            anchor = (info["reads_values"][j], info["bases_values"][j])
            head = f"{name}\n(n={n}, this study)\nmean:"
        else:
            anchor = (info["reads_values"][0], info["bases_values"][0])
            head = f"{name}\n(this study)"
        ax.annotate(
            f"{head}\n{fmt_rate(info['reads_mean'])} reads/fg\n"
            f"{fmt_rate(info['bases_mean'])} bases/fg",
            xy=anchor,
            xytext=(10 ** (np.log10(info["reads_mean"]) + dx),
                    10 ** (np.log10(info["bases_mean"]) + dy)),
            fontsize=8, fontweight="bold", ha=ha, va=va,
            bbox=dict(facecolor="white", edgecolor=THIS_STUDY_COLOR,
                      linewidth=0.8, pad=2.5, alpha=0.95),
            arrowprops=dict(arrowstyle="-", color=THIS_STUDY_COLOR, lw=0.6),
            zorder=7)

    # ---- marginal rugs -----------------------------------------------------
    y_hi, x_hi = np.log10(ax.get_ylim()[1]), np.log10(ax.get_xlim()[1])
    for _, r in prior.iterrows():
        c = color_of(r["study"])
        ax.plot([r["reads_per_fg"]] * 2, [10 ** (y_hi - 0.10), 10 ** (y_hi - 0.04)],
                color=c, linewidth=0.8, alpha=0.7, zorder=2)
        ax.plot([10 ** (x_hi - 0.10), 10 ** (x_hi - 0.04)], [r["bases_per_fg"]] * 2,
                color=c, linewidth=0.8, alpha=0.7, zorder=2)

    # ---- axes / legend -----------------------------------------------------
    ax.set_xlabel("Reads / fg DNA into library prep")
    ax.set_ylabel("Bases / fg DNA into library prep")
    ax.set_title(title)

    handles = [Line2D([0], [0], marker=marker_of(s), color=color_of(s),
                      markeredgecolor="black", markeredgewidth=0.5,
                      linestyle="None", markersize=6.5, label=s)
               for s in prior["study"].drop_duplicates()]
    for name in rounds:
        st = ROUND_STYLE.get(name, FALLBACK_ROUND_STYLE)
        handles.append(Line2D([0], [0], marker=st["marker"], color=THIS_STUDY_COLOR,
                              markeredgecolor="white", markeredgewidth=0.7,
                              linestyle="None", markersize=st["legend_ms"],
                              label=f"{name} (this study)"))
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.965, 0.02),
              frameon=True, framealpha=0.95, edgecolor="0.7")

    ax.grid(True, which="major", color="0.93", linewidth=0.5, zorder=0)
    ax.grid(True, which="minor", color="0.97", linewidth=0.4, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_color("0.6")

    fig.tight_layout()
    return fig, callouts, bpr_median, ref_name


# ---------------------------------------------------------------------------
# Sidecar helpers
# ---------------------------------------------------------------------------
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def summarise_round(df: pd.DataFrame) -> dict:
    reads = [float(v) for v in df["reads_per_fg"]]
    bases = [float(v) for v in df["bases_per_fg"]]
    return {
        "n": len(reads),
        "reads_values": reads,
        "bases_values": bases,
        "reads_mean": float(np.mean(reads)),
        "bases_mean": float(np.mean(bases)),
        "reads_min": float(np.min(reads)),
        "reads_max": float(np.max(reads)),
        "bases_min": float(np.min(bases)),
        "bases_max": float(np.max(bases)),
        "sources": sorted(set(df["source"])),
    }


def build_caption(prior, rounds, ref_name, args, callouts) -> str:
    ref = rounds[ref_name]
    studies = ", ".join(sorted(set(prior["study"])))
    folds = "; ".join(
        f"{c['study']} ({c['condition']}): {c['fold_reads']:,.0f}× reads, "
        f"{c['fold_bases']:,.0f}× bases" for c in callouts)
    return (
        "Low-input nanopore sequencing performance expressed per femtogram of DNA "
        "entering library preparation. Each prior-study point is one sample "
        f"({len(prior)} samples from {studies}); reads/fg = reads/(pg DNA x 1000) "
        "and bases/fg = bases/(pg DNA x 1000). This study is shown as "
        + " and ".join(f"{k} (n={v['n']})" for k, v in rounds.items())
        + f". Dashed contours are iso-improvement circles around the {ref_name} mean "
        f"({fmt_rate(ref['reads_mean'])} reads/fg, {fmt_rate(ref['bases_mean'])} "
        "bases/fg) at 10x, 100x and 1000x geometric-mean fold improvement in "
        "log10-log10 space. Callouts give the largest per-study improvement: "
        f"{folds}. Zorzano et al. 2025 counts use the '{args.zorzano_classifier}' "
        f"variant and Basapathi Raghavendra et al. 2023 the "
        f"'{args.raghavendra_classifier}' variant; see prior_studies.tsv for the "
        "per-row provenance of every point.")


# ---------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prior-tsv", type=Path, default=cd.PRIOR_TSV,
                    help="prior-study table (default: %(default)s)")
    ap.add_argument("--this-study-tsv", type=Path, default=cd.THIS_STUDY_TSV,
                    help="seeded this-study table (default: %(default)s)")
    ap.add_argument("--results-dir", type=Path, default=None,
                    help="pipeline results/ directory; headline rows of "
                         "<sample_id>/<mode>/<sample_id>.metrics.tsv supersede the "
                         "seeded this-study values (default: use the seeded values)")
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                    help="output directory (default: %(default)s)")
    ap.add_argument("--basename", default=DEFAULT_BASENAME,
                    help="output file stem (default: %(default)s)")
    ap.add_argument("--zorzano-classifier", choices=cd.ZORZANO_CHOICES,
                    default="kraken2_q1",
                    help="DEFECT (a): which Zorzano et al. 2025 read/base assignment "
                         "to plot. 'kraken2_q1' (default) is internally consistent on "
                         "both axes; 'published_squeezemeta' uses the paper's own hits "
                         "on both axes (bases are a derived estimate); "
                         "'legacy_hybrid_workbook' reproduces the old spreadsheet "
                         "figure, which mixed classifiers between axes")
    ap.add_argument("--raghavendra-classifier", choices=cd.RAGHAVENDRA_CHOICES,
                    default="kraken2_q10",
                    help="DEFECT (c): 'kraken2_q10' (default) is our reanalysis of the "
                         "deposited reads; 'published' uses the paper's pass reads, "
                         "which have no base counts and so cannot be plotted")
    ap.add_argument("--title",
                    default="Low-input nanopore performance — this study vs prior work",
                    help="figure title")
    ap.add_argument("--dpi", type=int, default=600, help="raster DPI (default: %(default)s)")
    ap.add_argument("--allow-unverified", action="store_true", default=True,
                    help=argparse.SUPPRESS)
    ap.add_argument("--drop-unverified", dest="allow_unverified", action="store_false",
                    help="exclude rows whose counts have no recorded citation "
                         "(DEFECT (b): the Mojarro et al. 2019 point)")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    prior, prior_notes = cd.load_prior(args.prior_tsv,
                                       zorzano_classifier=args.zorzano_classifier,
                                       raghavendra_classifier=args.raghavendra_classifier)
    mine, mine_notes = cd.load_this_study(args.this_study_tsv, args.results_dir)

    for n in prior_notes + mine_notes:
        print(f"[data] {n}")

    # ---- DEFECT (b): shout about unverified rows ---------------------------
    unverified = prior[~prior["verified"]]
    if len(unverified) and not args.allow_unverified:
        print(f"[data] --drop-unverified: excluding {len(unverified)} unverified row(s).")
        prior = prior[prior["verified"]].copy()
        unverified = unverified.iloc[0:0]
    if len(unverified):
        bar = "!" * 78
        print(f"\n{bar}", file=sys.stderr)
        print(f"!! WARNING: plotting {len(unverified)} UNVERIFIED point(s) "
              "(no citation on record)", file=sys.stderr)
        for _, r in unverified.iterrows():
            print(f"!!   {r['study']} | {r['condition']} | classifier={r['classifier']}",
                  file=sys.stderr)
            print(f"!!   reads={r['reads']:.0f} bases={r['bases']:.0f} "
                  f"dna_pg={r['dna_pg']:g} -> {r['reads_per_fg']:g} reads/fg, "
                  f"{r['bases_per_fg']:g} bases/fg", file=sys.stderr)
            print(f"!!   {r['provenance_note']}", file=sys.stderr)
        print("!! Trace these counts to a published table before submission, or "
              "re-run with --drop-unverified.", file=sys.stderr)
        print(f"{bar}\n", file=sys.stderr)

    derived = prior[prior["source"] == "derived_estimate"]
    if len(derived):
        print(f"[data] NOTE: {len(derived)} plotted point(s) use a derived "
              "(not published) base-count estimate; see provenance_note.")

    if prior.empty:
        sys.exit("[error] no plottable prior-study rows for the selected classifiers.")
    if mine.empty:
        sys.exit("[error] no this-study rows to plot.")

    # ---- rounds, data-driven ----------------------------------------------
    mine = mine.dropna(subset=["reads_per_fg", "bases_per_fg"])
    order = [r for r in ("Round 1", "Round 2") if r in set(mine["round"])]
    order += [r for r in dict.fromkeys(mine["round"]) if r not in order]
    rounds = {name: summarise_round(mine[mine["round"] == name].sort_values("replicate_idx"))
              for name in order}
    for name, info in rounds.items():
        print(f"[data] {name}: n={info['n']}, mean {info['reads_mean']:.6g} reads/fg, "
              f"{info['bases_mean']:.6g} bases/fg, source(s)={','.join(info['sources'])}")

    # ---- figure ------------------------------------------------------------
    fig, callouts, bpr_median, ref_name = build_figure(prior, rounds, args.title)

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    pdf = outdir / f"{args.basename}.pdf"
    png = outdir / f"{args.basename}.png"
    csv = outdir / f"{args.basename}.csv"
    js = outdir / f"{args.basename}.json"

    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=args.dpi)
    plt.close(fig)

    # ---- CSV of every plotted point ---------------------------------------
    ref = rounds[ref_name]
    cols = ["series", "study", "study_short", "condition", "organism",
            "replicate_idx", "replicate_label", "classifier", "reads", "bases",
            "dna_pg", "reads_per_fg", "bases_per_fg", "verified", "source",
            "source_detail", "provenance_note"]
    plotted = pd.concat([
        prior.assign(series="prior study"),
        mine.assign(series=mine["round"]),
    ], ignore_index=True)
    plotted = plotted.reindex(columns=cols)
    plotted[f"fold_reads_vs_{ref_name.lower().replace(' ', '')}_mean"] = \
        ref["reads_mean"] / plotted["reads_per_fg"]
    plotted[f"fold_bases_vs_{ref_name.lower().replace(' ', '')}_mean"] = \
        ref["bases_mean"] / plotted["bases_per_fg"]
    # repr() is the shortest text that round-trips a float64 exactly; pandas'
    # default serialiser drops the last significant digit on mixed-dtype frames.
    plotted.to_csv(csv, index=False, float_format=lambda v: repr(float(v)))

    # ---- JSON sidecar ------------------------------------------------------
    def rel(p):
        """Repo-relative path when possible, absolute otherwise."""
        p = Path(p).resolve()
        try:
            return str(p.relative_to(HERE.parent))
        except ValueError:
            return str(p)

    source_files = [
        {"path": rel(args.prior_tsv),
         "sha256": sha256(args.prior_tsv), "role": "prior-study data"},
        {"path": rel(args.this_study_tsv),
         "sha256": sha256(args.this_study_tsv), "role": "this-study data"},
        {"path": f"comparison/{Path(__file__).name}",
         "sha256": sha256(Path(__file__)), "role": "figure script"},
        {"path": "comparison/comparison_data.py",
         "sha256": sha256(HERE / "comparison_data.py"), "role": "loader"},
    ]
    if args.results_dir is not None:
        source_files.append({"path": str(args.results_dir), "sha256": None,
                             "role": "live pipeline results directory"})

    sidecar = {
        "id": DISPLAY_ITEM_ID,
        "title": args.title,
        "caption": build_caption(prior, rounds, ref_name, args, callouts),
        "source_files": source_files,
        "software": {
            "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
            "command": " ".join([Path(__file__).name] + (list(argv) if argv else sys.argv[1:])),
        },
        "metrics": {
            "convention": {
                "reads_per_fg": "reads / (dna_pg * 1000)",
                "bases_per_fg": "bases / (dna_pg * 1000)",
                "dna_basis": "DNA into library prep (post-extraction)",
            },
            "classifier_selection": {
                "Zorzano et al. 2025": args.zorzano_classifier,
                "B. Raghavendra et al. 2023": args.raghavendra_classifier,
                "Mojarro et al. 2019": "published_table1",
            },
            "n_points_plotted": int(len(prior) + len(mine)),
            "n_prior_points": int(len(prior)),
            "n_prior_points_by_study": {k: int(v) for k, v in
                                        prior["study"].value_counts().items()},
            "rounds": {k: {kk: vv for kk, vv in v.items()} for k, v in rounds.items()},
            "reference_round": ref_name,
            "median_bases_per_read_prior": bpr_median,
            "largest_improvement_per_study": [
                {"study": c["study"], "condition": c["condition"],
                 "fold_reads": c["fold_reads"], "fold_bases": c["fold_bases"]}
                for c in callouts],
            "unverified_points_plotted": [
                {"study": r["study"], "condition": r["condition"],
                 "classifier": r["classifier"], "reads_per_fg": float(r["reads_per_fg"]),
                 "bases_per_fg": float(r["bases_per_fg"])}
                for _, r in unverified.iterrows()],
        },
    }
    js.write_text(json.dumps(sidecar, indent=2) + "\n", encoding="utf-8")

    for p in (pdf, png, csv, js):
        print(f"[saved] {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
