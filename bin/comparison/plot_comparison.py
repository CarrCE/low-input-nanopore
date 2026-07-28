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

Inputs are the versioned TSVs under assets/comparison/:
    prior_studies.tsv   one row per prior-study sample per classifier variant,
                        each citing the published table or reanalysis it rests on
    this_study.tsv      Round 1 / Round 2 values, a committed snapshot of
                        pipeline output (see seed_this_study.py)
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

HERE = Path(__file__).resolve().parent          # bin/comparison
REPO = HERE.parent.parent                       # repository root
# Outputs go under results/ with everything else the pipeline produces;
# results/ is gitignored, so figures are regenerated rather than tracked.
DEFAULT_OUTDIR = REPO / "results" / "comparison"
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

# Label placement is solved, not tabulated. The previous version carried
# hand-tuned per-study offsets in log10 units; they were correct for the data as
# it stood when they were tuned, and silently wrong afterwards -- correcting the
# Round 1 analysis moved that cluster on top of the "bp/read" guide, and the
# Round 1 box came to sit over two Basapathi Raghavendra points. Offsets that
# have to be re-tuned whenever a number changes are a defect in a figure that is
# regenerated from live pipeline output. `Placer` below chooses positions by
# measuring overlap instead, so the figure stays legible when the data move.
PLACER_CANDIDATE_ANGLES = 24
PLACER_CANDIDATE_RADII = (0.45, 0.70, 1.00, 1.35, 1.75)

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
# Label placement
# ---------------------------------------------------------------------------
class Placer:
    """Chooses label positions that do not cover data.

    Everything is done in log10 data units, which on these log-log axes are
    linear and so let a text box be treated as a plain rectangle. Obstacles are
    the plotted points, the marginal rugs, and every label already placed --
    labels are registered as they are positioned, so later ones route around
    earlier ones rather than stacking on them.

    The scoring is deliberately blunt: overlapping a data point is a hard cost,
    overlapping another label a slightly softer one, and among the candidates
    that clear both, the one nearest its anchor wins so leader lines stay short.
    """

    #: A point marker is about this many log10 units across; used to give
    #: scatter points a footprint rather than treating them as infinitesimal.
    POINT_PAD = 0.055

    def __init__(self, fig, ax, points_xy):
        self.fig, self.ax = fig, ax
        (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
        self.xlim = (np.log10(x0), np.log10(x1))
        self.ylim = (np.log10(y0), np.log10(y1))

        # Font metrics are in points; obstacles are in log10 data units. The
        # bridge is the axes box in inches, so this must run after the layout is
        # settled or every box is sized against the wrong axes width.
        bbox = ax.get_position()
        w_in = fig.get_size_inches()[0] * bbox.width
        h_in = fig.get_size_inches()[1] * bbox.height
        self.logx_per_in = (self.xlim[1] - self.xlim[0]) / w_in
        self.logy_per_in = (self.ylim[1] - self.ylim[0]) / h_in

        self.obstacles = [self._point_box(x, y) for x, y in points_xy]
        # The rug ticks live in the top and right margins; nothing may sit there.
        self.obstacles.append((self.xlim[0], self.ylim[1] - 0.12,
                               self.xlim[1], self.ylim[1]))
        self.obstacles.append((self.xlim[1] - 0.12, self.ylim[0],
                               self.xlim[1], self.ylim[1]))
        self.labels = []

    def _point_box(self, x, y):
        lx, ly = np.log10(x), np.log10(y)
        return (lx - self.POINT_PAD, ly - self.POINT_PAD,
                lx + self.POINT_PAD, ly + self.POINT_PAD)

    def text_halfsize(self, text, fontsize, pad_pt=3.0):
        """Half-width and half-height of a rendered text box, in log10 units."""
        lines = text.split("\n")
        # 0.58 em average advance for this sans stack; mathtext italics run a
        # little wider, which the padding absorbs.
        w_in = (max(len(l) for l in lines) * 0.58 * fontsize + 2 * pad_pt) / 72.0
        h_in = (len(lines) * 1.30 * fontsize + 2 * pad_pt) / 72.0
        return (w_in * self.logx_per_in / 2.0, h_in * self.logy_per_in / 2.0)

    @staticmethod
    def _overlap(a, b):
        """Area of intersection of two (x0, y0, x1, y1) rectangles."""
        dx = min(a[2], b[2]) - max(a[0], b[0])
        dy = min(a[3], b[3]) - max(a[1], b[1])
        return dx * dy if (dx > 0 and dy > 0) else 0.0

    @staticmethod
    def _segment_hits(p0, p1, rect, n=32):
        """Does the leader from p0 to p1 pass through rect?

        Sampled rather than clipped: at this resolution a box narrower than the
        sample spacing would have to be smaller than the text it contains, and
        sampling keeps the test to three lines.
        """
        for t in np.linspace(0.0, 1.0, n):
            x = p0[0] + t * (p1[0] - p0[0])
            y = p0[1] + t * (p1[1] - p0[1])
            if rect[0] <= x <= rect[2] and rect[1] <= y <= rect[3]:
                return True
        return False

    def _cost(self, box, anchor_log):
        cost = 0.0
        for o in self.obstacles:
            cost += 60.0 * self._overlap(box, o)
        for l in self.labels:
            cost += 25.0 * self._overlap(box, l)

        # A box can sit in clear space and still ruin a neighbour by dragging its
        # leader line across it -- which is how the Round 2 leader came to strike
        # through the Round 1 box. Charge for that too.
        centre = ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2)
        for l in self.labels:
            if self._segment_hits(anchor_log, centre, l):
                cost += 3.0
        for o in self.obstacles[:-2]:            # points only, not the rug strips
            if self._segment_hits(anchor_log, centre, o):
                cost += 0.4
        # Staying inside the axes matters more than anything except covering a
        # point: a box clipped by the frame is unreadable.
        outside = (max(0.0, self.xlim[0] - box[0]) + max(0.0, box[2] - self.xlim[1])
                   + max(0.0, self.ylim[0] - box[1]) + max(0.0, box[3] - self.ylim[1]))
        cost += 120.0 * outside
        cx, cy = (box[0] + box[2]) / 2, (box[1] + box[3]) / 2
        cost += 0.35 * np.hypot(cx - anchor_log[0], cy - anchor_log[1])
        return cost

    def place(self, anchor_xy, text, fontsize, radii=PLACER_CANDIDATE_RADII,
              angles=PLACER_CANDIDATE_ANGLES, register=True):
        """Best (x, y) in data coordinates for a label anchored at anchor_xy."""
        hw, hh = self.text_halfsize(text, fontsize)
        ax_, ay_ = np.log10(anchor_xy[0]), np.log10(anchor_xy[1])
        best, best_cost = None, np.inf
        for r in radii:
            for t in np.linspace(0, 2 * np.pi, angles, endpoint=False):
                cx, cy = ax_ + r * np.cos(t) * 1.35, ay_ + r * np.sin(t)
                box = (cx - hw, cy - hh, cx + hw, cy + hh)
                c = self._cost(box, (ax_, ay_))
                if c < best_cost:
                    best, best_cost = (cx, cy, box), c
        cx, cy, box = best
        if register:
            self.labels.append(box)
        return 10 ** cx, 10 ** cy

    def place_on_path(self, xs_log, ys_log, text, fontsize, register=True):
        """Best position from a set of candidate points along a path."""
        hw, hh = self.text_halfsize(text, fontsize)
        best, best_cost = None, np.inf
        for lx, ly in zip(xs_log, ys_log):
            box = (lx - hw, ly - hh, lx + hw, ly + hh)
            c = self._cost(box, (lx, ly))
            if c < best_cost:
                best, best_cost = (lx, ly, box), c
        lx, ly, box = best
        if register:
            self.labels.append(box)
        return 10 ** lx, 10 ** ly

    def register_rect(self, x0, y0, x1, y1):
        """Reserve a region in data coordinates (used for the legend)."""
        self.labels.append((np.log10(x0), np.log10(y0),
                            np.log10(x1), np.log10(y1)))


# ---------------------------------------------------------------------------
# Figure pieces
# ---------------------------------------------------------------------------
def draw_iso_improvement_contour(ax, cx, cy, log_radius,
                                 color="0.4", lw=0.7, ls=(0, (3, 3))):
    """Circle of `log_radius` log10 units around (cx, cy), drawn on log-log axes."""
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(10 ** (np.log10(cx) + log_radius * np.cos(theta)),
            10 ** (np.log10(cy) + log_radius * np.sin(theta)),
            color=color, linewidth=lw, linestyle=ls, alpha=0.55, zorder=1)


def label_iso_contour(ax, placer, cx, cy, log_radius, label, color="0.4"):
    """Put the contour's label on the emptiest stretch of its own circle."""
    theta = np.linspace(0, 2 * np.pi, 180, endpoint=False)
    xs = np.log10(cx) + log_radius * np.cos(theta)
    ys = np.log10(cy) + log_radius * np.sin(theta)
    x, y = placer.place_on_path(xs, ys, label, 7.5)
    ax.text(x, y, label, color=color, fontsize=7.5, ha="center", va="center",
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


def verify_placement(fig, ax, artists, points_xy, pad_px=1.5):
    """Report any label that ends up covering a plotted point.

    The placer works from estimated text metrics; this checks the rendered
    result, so a bad estimate surfaces as a warning at generation time instead
    of as a reviewer's comment. Returns a list of (label, x, y) collisions.
    """
    fig.canvas.draw()
    hits = []
    for art, name in artists:
        # The drawn box, not Annotation.get_window_extent(): that returns the
        # union of the text and its leader, so a label with a long leader claims
        # a rectangle spanning everything the leader flies over and reports
        # collisions with points it does not touch. The bbox patch is the white
        # rectangle the reader actually sees.
        patch = art.get_bbox_patch()
        try:
            bb = patch.get_window_extent() if patch is not None else art.get_window_extent()
        except Exception:                       # artist never rendered
            continue
        for x, y in points_xy:
            px, py = ax.transData.transform((x, y))
            if (bb.x0 - pad_px <= px <= bb.x1 + pad_px
                    and bb.y0 - pad_px <= py <= bb.y1 + pad_px):
                hits.append((name, float(x), float(y)))
    return hits


def build_figure(prior: pd.DataFrame, rounds: dict, title: str | None):
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

    ref_name = next(iter(rounds))
    ref = {"reads": rounds[ref_name]["reads_mean"], "bases": rounds[ref_name]["bases_mean"]}

    # ---- axis limits -------------------------------------------------------
    all_x = list(prior["reads_per_fg"]) + [v for r in rounds.values() for v in r["reads_values"]]
    all_y = list(prior["bases_per_fg"]) + [v for r in rounds.values() for v in r["bases_values"]]
    pad, pad_y_lo = 0.45, 0.95   # extra room at the bottom for callouts
    ax.set_xlim(10 ** (np.log10(min(all_x)) - pad), 10 ** (np.log10(max(all_x)) + pad))
    ax.set_ylim(10 ** (np.log10(min(all_y)) - pad_y_lo), 10 ** (np.log10(max(all_y)) + pad))
    ax.set_xscale("log")
    ax.set_yscale("log")

    # ---- axes, legend, layout ----------------------------------------------
    # These come before any label is positioned. The legend occupies a real
    # patch of the axes and the layout fixes the axes size in inches, and the
    # placer needs both: without the layout it sizes every text box against the
    # wrong axes width, and without the legend it will happily park a callout
    # underneath it.
    ax.set_xlabel("Reads / fg DNA into library prep")
    ax.set_ylabel("Bases / fg DNA into library prep")
    if title:
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
    legend = ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.965, 0.02),
                       frameon=True, framealpha=0.95, edgecolor="0.7")

    ax.grid(True, which="major", color="0.93", linewidth=0.5, zorder=0)
    ax.grid(True, which="minor", color="0.97", linewidth=0.4, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_color("0.6")

    fig.tight_layout()
    fig.canvas.draw()          # settle the layout so text metrics are meaningful

    # ---- placement ---------------------------------------------------------
    points = ([(x, y) for x, y in zip(prior["reads_per_fg"], prior["bases_per_fg"])]
              + [(x, y) for r in rounds.values()
                 for x, y in zip(r["reads_values"], r["bases_values"])])
    placer = Placer(fig, ax, points)
    placed_artists = []

    lb = legend.get_window_extent().transformed(ax.transData.inverted())
    placer.register_rect(lb.x0, lb.y0, lb.x1, lb.y1)

    # ---- iso-improvement contours around the reference round ---------------
    for log_r, label in [(1, "10× less"), (2, "100× less"), (3, "1000× less")]:
        draw_iso_improvement_contour(ax, ref["reads"], ref["bases"], log_r)
        label_iso_contour(ax, placer, ref["reads"], ref["bases"], log_r, label)

    # ---- constant bases-per-read reference ---------------------------------
    # The guide line runs corner to corner through the densest part of the plot,
    # so its label is placed along the line itself rather than pinned to one end
    # -- pinning it to the right end put it under this study's markers.
    bpr_median = float(np.median(prior["bases_per_fg"] / prior["reads_per_fg"]))
    xs = np.array(ax.get_xlim())
    ax.plot(xs, bpr_median * xs, color="0.7", linewidth=0.7,
            linestyle=(0, (1, 2)), zorder=0)
    bpr_text = f"~{bpr_median:,.0f} bp/read"
    lxs = np.linspace(np.log10(xs[0]) + 0.25, np.log10(xs[1]) - 0.25, 60)
    lys = np.log10(bpr_median) + lxs
    bx, by = placer.place_on_path(lxs, lys, bpr_text, 7.5)
    placed_artists.append((
        ax.text(bx, by, bpr_text, color="0.5", fontsize=7.5, ha="center", va="center",
                bbox=dict(facecolor="white", edgecolor="none", pad=0.8, alpha=0.85),
                zorder=2), bpr_text))

    # ---- per-study callouts ------------------------------------------------
    callouts = pick_callouts(prior, ref)
    for c in callouts:
        tx, ty = placer.place((c["x"], c["y"]), c["text"], 7.5)
        placed_artists.append((ax.annotate(
            c["text"], xy=(c["x"], c["y"]), xytext=(tx, ty),
            fontsize=7.5, ha="center", va="center", color=color_of(c["study"]),
            bbox=dict(facecolor="white", edgecolor=color_of(c["study"]),
                      linewidth=0.6, pad=2.5, alpha=0.92),
            arrowprops=dict(arrowstyle="-", color=color_of(c["study"]),
                            lw=0.6, alpha=0.7, connectionstyle="arc3,rad=0.12"),
            zorder=5), c["study"]))

    # ---- this-study labels -------------------------------------------------
    for name, info in rounds.items():
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
        text = (f"{head}\n{fmt_rate(info['reads_mean'])} reads/fg\n"
                f"{fmt_rate(info['bases_mean'])} bases/fg")
        tx, ty = placer.place(anchor, text, 8)
        placed_artists.append((ax.annotate(
            text, xy=anchor, xytext=(tx, ty),
            fontsize=8, fontweight="bold", ha="center", va="center",
            bbox=dict(facecolor="white", edgecolor=THIS_STUDY_COLOR,
                      linewidth=0.8, pad=2.5, alpha=0.95),
            arrowprops=dict(arrowstyle="-", color=THIS_STUDY_COLOR, lw=0.6),
            zorder=7), name))

    # ---- marginal rugs -----------------------------------------------------
    y_hi, x_hi = np.log10(ax.get_ylim()[1]), np.log10(ax.get_xlim()[1])
    for _, r in prior.iterrows():
        c = color_of(r["study"])
        ax.plot([r["reads_per_fg"]] * 2, [10 ** (y_hi - 0.10), 10 ** (y_hi - 0.04)],
                color=c, linewidth=0.8, alpha=0.7, zorder=2)
        ax.plot([10 ** (x_hi - 0.10), 10 ** (x_hi - 0.04)], [r["bases_per_fg"]] * 2,
                color=c, linewidth=0.8, alpha=0.7, zorder=2)

    collisions = verify_placement(fig, ax, placed_artists, points)
    if collisions:
        print("[layout] WARNING: label(s) overlap plotted points:", file=sys.stderr)
        for name, x, y in collisions:
            print(f"[layout]   {name!r} covers point ({x:.6g}, {y:.6g})", file=sys.stderr)
    else:
        print(f"[layout] {len(placed_artists)} labels placed, none covering a data point")

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
                         "on both axes (bases are a derived estimate)")
    ap.add_argument("--raghavendra-classifier", choices=cd.RAGHAVENDRA_CHOICES,
                    default="kraken2_q10",
                    help="DEFECT (c): 'kraken2_q10' (default) is our reanalysis of the "
                         "deposited reads; 'published' uses the paper's pass reads, "
                         "which have no base counts and so cannot be plotted")
    # The figure carries no drawn title: in the manuscript it sits above a
    # caption that says the same thing, and a title there is redundant and eats
    # plot area. The string is still the display item's title in the JSON
    # sidecar, which is what the display-item registry reads.
    ap.add_argument("--title",
                    default="Low-input nanopore performance — this study vs prior work",
                    help="display-item title recorded in the JSON sidecar")
    ap.add_argument("--draw-title", action="store_true",
                    help="also render the title above the axes (off by default; "
                         "the manuscript caption already carries it)")
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
    fig, callouts, bpr_median, ref_name = build_figure(
        prior, rounds, args.title if args.draw_title else None)

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
            return str(p.relative_to(REPO))
        except ValueError:
            return str(p)

    source_files = [
        {"path": rel(args.prior_tsv),
         "sha256": sha256(args.prior_tsv), "role": "prior-study data"},
        {"path": rel(args.this_study_tsv),
         "sha256": sha256(args.this_study_tsv), "role": "this-study data"},
        {"path": f"bin/comparison/{Path(__file__).name}",
         "sha256": sha256(Path(__file__)), "role": "figure script"},
        {"path": "bin/comparison/comparison_data.py",
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
            "n_plotted_points": int(len(prior) + len(mine)),
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
