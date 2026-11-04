#!/usr/bin/env python3
"""
Combine per-sample metrics into the study-level tables and display items.

Consumes every `<sample>.metrics.tsv` and `<sample>.summary.json` produced by
compute_metrics.py and emits:

    per_organism.tsv        every organism x sample row, all modes
    per_sample.tsv          one headline row per sample ("All organisms")
    experiment_summary.tsv  mean/SD across replicates, headline samples only
    abundance.pdf/.png      theoretical vs measured abundance, per experiment
    readlengths.pdf/.png    read-length distribution by role -- the
                            adaptive-sampling ejection signature
    <item>.csv/.json        data + provenance sidecar for each display item

Replicates flagged include_in_headline=0 in assets/measurements.tsv (e.g.
lowinput_s2_r0, whose mass is an extrapolation below the Qubit HS reporting
range rather than a reading within it) are carried in the per-sample table but
excluded from the experiment-level statistics. The exclusion rests on the
recorded basis of the mass, not on the mass being missing -- r0 has one.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Okabe-Ito, colourblind-safe; matches the comparison figure's palette.
ROLE_COLORS = {
    "sample":      "#009E73",
    "carrier":     "#0072B2",
    "contaminant": "#D55E00",
    "ambiguous":   "#CC79A7",
    "none":        "#999999",   # unmapped; assign_reads.py labels these "none"
    "unassigned":  "#999999",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "dejavusans",
})


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metrics", nargs="+", required=True,
                   help="per-sample *.metrics.tsv files")
    p.add_argument("--summaries", nargs="+", default=[],
                   help="per-sample *.summary.json files")
    p.add_argument("--readlengths", nargs="+", default=[],
                   help="per-sample *.readlengths.tsv.gz files")
    p.add_argument("--measurements", default=None,
                   help="assets/measurements.tsv; authoritative for include_in_headline")
    p.add_argument("--samplesheet", default=None,
                   help="fallback include_in_headline for samples measurements.tsv "
                        "does not name (the smoke test and the reanalysed prior study)")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_headline_flags(path, delimiter="\t"):
    """Which replicates contribute to the headline statistics.

    Read from measurements.tsv rather than from a samplesheet: the flag belongs
    with the mass whose provenance justifies it. An exclusion must be an
    explicit 0, never inferred from a value being blank -- inferring it is how a
    replicate silently rejoined a published mean the moment its mass was filled
    in.
    """
    if not path:
        return {}
    lines = [l for l in Path(path).read_text().splitlines()
             if l.strip() and not l.strip().startswith("#")]
    if not lines:
        return {}
    import csv
    import io
    rows = list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter=delimiter))
    flags = {}
    for r in rows:
        sid = r.get("sample_id")
        if not sid:
            continue
        # `== "1"` was the old test, which excluded a replicate whose flag was
        # blank -- the opposite of the rule stated above, and the opposite of
        # what comparison_data.py does with the same column, so the two could
        # select different replicate sets from one file. A blank means "not
        # stated", which is an inclusion; only an explicit 0 excludes.
        raw = str(r.get("include_in_headline") or "").strip()
        if raw not in ("", "0", "1"):
            sys.exit(f"error: {path}: {sid} has include_in_headline={raw!r}; "
                     f"expected 0, 1, or blank.")
        flags[sid] = raw != "0"
    return flags


def read_experiment_names(path):
    """Display name for each experiment: its standard's catalogue number.

    Figures say "D6311", not the internal id "lowinput_s1". The name comes from
    the standard_catalogue column of measurements.tsv, so it lives beside the
    other facts about each replicate instead of in a mapping in plotting code.
    An experiment the file does not name (the smoke test, the reanalysed prior
    study) keeps its id. A pack-size suffix is dropped: lowinput_s2 used D6321
    and D6321-10, "the same product in a ten-pack" (measurements.tsv header),
    and the product is what the figure names.

    A copy of bin/display_names.py:read_experiment_names(), inlined because
    the AGGREGATE process stages this script alone; keep the two in step.
    """
    if not path:
        return {}
    import csv
    lines = [l for l in Path(path).read_text().splitlines()
             if l.strip() and not l.lstrip().startswith("#")]
    names = {}
    for r in csv.DictReader(lines, delimiter="\t"):
        exp = (r.get("experiment") or "").strip()
        cat = (r.get("standard_catalogue") or "").strip().split("-")[0]
        if exp and cat:
            if names.get(exp, cat) != cat:
                sys.exit(f"error: {path}: experiment {exp} lists two standards "
                         f"({names[exp]} and {cat}); its display name is ambiguous.")
            names[exp] = cat
    return names


def place_labels(ax, points, texts, fontsize=6.5, color="0.25", obstacles=()):
    """Label each point at the first nearby position that overlaps nothing.

    Points along the diagonal of a log-distributed community crowd together, so
    a fixed offset (or alternating sides) still lets neighbours overprint. Each
    label tries eight positions around its point, nearest first, then the same
    eight further out, and keeps the first whose box clears every label already
    placed, every marker and error bar, and the axes edge. If none clears, it
    takes the position with the least overlap.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_px = ax.transData.transform
    axbox = ax.get_window_extent(renderer)
    taken = list(obstacles)
    for (x, y) in points:                       # every marker is an obstacle
        px, py = to_px((x, y))
        taken.append((px - 5, py - 5, px + 5, py + 5))

    def overlap(a, b):
        w = min(a[2], b[2]) - max(a[0], b[0])
        h = min(a[3], b[3]) - max(a[1], b[1])
        return w * h if w > 0 and h > 0 else 0.0

    base = [(6, 4, "left", "bottom"), (6, -4, "left", "top"),
            (-6, 4, "right", "bottom"), (-6, -4, "right", "top"),
            (8, 0, "left", "center"), (-8, 0, "right", "center"),
            (0, 8, "center", "bottom"), (0, -8, "center", "top"),
            # steep diagonals, to clear a neighbour's tall error bar
            (4, 10, "left", "bottom"), (4, -10, "left", "top"),
            (-4, 10, "right", "bottom"), (-4, -10, "right", "top")]
    # Three rings. Beyond the first, a thin leader line ties the label to its
    # point, because a label that far out could otherwise be read as naming a
    # neighbour.
    cands = [(dx * k, dy * k, ha, va, k > 1)
             for k in (1, 2.2, 3.5) for dx, dy, ha, va in base]
    dpi_scale = fig.dpi / 72.0
    for (x, y), text in zip(points, texts):
        best = None
        for dx, dy, ha, va, leader in cands:
            # Measure without the leader line: an Annotation's box includes its
            # arrow, which always reaches back to the point and so would count
            # the label's own marker as an overlap.
            ann = ax.annotate(text, (x, y), textcoords="offset points",
                              xytext=(dx, dy), ha=ha, va=va,
                              fontsize=fontsize, color=color)
            bb = ann.get_window_extent(renderer)
            ann.remove()
            box = (bb.x0 - 1.5 * dpi_scale, bb.y0 - 1.5 * dpi_scale,
                   bb.x1 + 1.5 * dpi_scale, bb.y1 + 1.5 * dpi_scale)
            outside = (box[0] < axbox.x0 or box[2] > axbox.x1
                       or box[1] < axbox.y0 or box[3] > axbox.y1)
            cost = sum(overlap(box, t) for t in taken) + (1e9 if outside else 0)
            if best is None or cost < best[0]:
                best = (cost, (dx, dy, ha, va, leader), box)
            if cost == 0:
                break
        dx, dy, ha, va, leader = best[1]
        ax.annotate(text, (x, y), textcoords="offset points", xytext=(dx, dy),
                    ha=ha, va=va, fontsize=fontsize, color=color,
                    arrowprops=(dict(arrowstyle="-", color="0.6", lw=0.5,
                                     shrinkA=1, shrinkB=4) if leader else None))
        taken.append(best[2])


def italicize(name):
    """Species names render italic; everything else stays upright."""
    if not name or name[0].islower() or " " not in name:
        return name
    if name.startswith(("ambiguous", "All ", "unassigned")):
        return name
    return "$\\mathit{" + name.replace(" ", r"\ ") + "}$"


def write_sidecar(path, item_id, title, caption, sources, metrics=None):
    payload = {
        "id": item_id,
        "title": title,
        "caption": caption,
        "source_files": sorted(sources),
        "software": ["python 3.12", f"pandas {pd.__version__}",
                     f"numpy {np.__version__}", f"matplotlib {matplotlib.__version__}"],
    }
    if metrics:
        payload["metrics"] = metrics
    Path(path).write_text(json.dumps(payload, indent=2))


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # measurements.tsv wins wherever it speaks; the samplesheet answers only for
    # samples it does not name.
    headline_flags = read_headline_flags(args.samplesheet, delimiter=",")
    headline_flags.update(read_headline_flags(args.measurements))

    exp_names = read_experiment_names(args.measurements)

    per_org = pd.concat([pd.read_csv(f, sep="\t") for f in args.metrics],
                        ignore_index=True)
    per_org["include_in_headline"] = per_org["sample_id"].map(
        lambda s: headline_flags.get(s, True))
    per_org.to_csv(outdir / "per_organism.tsv", sep="\t", index=False)

    # ---- headline per-sample table ----------------------------------------
    per_sample = per_org[per_org["organism"] == "All organisms"].copy()

    enrich = {}
    for f in args.summaries:
        d = json.loads(Path(f).read_text())
        enrich[(d["sample_id"], d["mode"])] = {
            "total_reads": d["totals"]["reads"],
            "total_bases": d["totals"]["bases"],
            "carrier_reads": d["by_role"]["carrier"]["reads"],
            "carrier_bases": d["by_role"]["carrier"]["bases"],
            "contaminant_reads": d["by_role"]["contaminant"]["reads"],
            "contaminant_bases": d["by_role"]["contaminant"]["bases"],
            "ambiguous_reads": d["by_role"]["ambiguous"]["reads"],
            "unassigned_reads": d["by_role"]["unassigned"]["reads"],
            "input_sample_fraction": d["input"]["input_sample_fraction"],
            "output_sample_fraction_bases": d["output_sample_fraction"]["bases"],
            "enrichment_bases": d["enrichment"]["bases"],
            "enrichment_reads": d["enrichment"]["reads"],
        }
    for col in ["total_reads", "total_bases", "carrier_reads", "carrier_bases",
                "contaminant_reads", "contaminant_bases", "ambiguous_reads",
                "unassigned_reads", "input_sample_fraction",
                "output_sample_fraction_bases", "enrichment_bases", "enrichment_reads"]:
        per_sample[col] = per_sample.apply(
            lambda r: enrich.get((r["sample_id"], r["mode"]), {}).get(col), axis=1)

    per_sample = per_sample.sort_values(["experiment", "mode", "replicate"])
    per_sample.to_csv(outdir / "per_sample.tsv", sep="\t", index=False)

    # ---- experiment-level statistics --------------------------------------
    hl = per_sample[per_sample["include_in_headline"]]
    agg_cols = ["reads_per_fg", "bases_per_fg", "enrichment_bases",
                "enrichment_reads", "output_sample_fraction_bases"]
    if len(hl):
        summary = (hl.groupby(["experiment", "mode"])[agg_cols]
                     .agg(["count", "mean", "std", "min", "max"]))
        summary.columns = ["_".join(c) for c in summary.columns]
        summary = summary.reset_index()
        summary.to_csv(outdir / "experiment_summary.tsv", sep="\t", index=False)
        print(summary.to_string(index=False))
    else:
        print("warning: no headline samples; experiment_summary.tsv not written")
        summary = pd.DataFrame()

    excluded = per_sample.loc[~per_sample["include_in_headline"], "sample_id"].unique()
    if len(excluded):
        print(f"[aggregate] excluded from experiment statistics (no quantified "
              f"input mass): {', '.join(sorted(excluded))}")

    # ---- Figure: theoretical vs measured abundance ------------------------
    samples = per_org[(per_org["role"] == "sample")
                      & per_org["theoretical_dna_fraction"].notna()
                      & (per_org["theoretical_dna_fraction"] > 0)].copy()
    samples = samples[samples["measured_base_fraction"].notna()]

    # With --mode both, per_org holds one row per (sample, mode, organism), so
    # averaging without first selecting a mode silently mixes the two
    # assignment rules: the "mean across replicates" becomes the mean of the
    # competitive and the sequential value, and the error bar becomes
    # between-mode variance rather than between-replicate. The panel title
    # counts distinct sample_ids and so would still read n=3 while the
    # statistic behind it was over 6 rows. Pick one mode -- competitive, which
    # every headline number uses -- and fall back to whichever single mode is
    # present when only one was run. MAKE_READLENGTHS applies the same
    # deduplication upstream (main.nf) for the same reason.
    modes = sorted(samples["mode"].unique())
    abundance_mode = "competitive" if "competitive" in modes else (
        modes[0] if modes else None)
    if abundance_mode is not None:
        samples = samples[samples["mode"] == abundance_mode].copy()
        if len(modes) > 1:
            print(f"[aggregate] abundance figure uses mode={abundance_mode} "
                  f"(available: {', '.join(modes)})")

    if len(samples):
        experiments = sorted(samples["experiment"].unique())
        fig, axes = plt.subplots(1, len(experiments),
                                 figsize=(5.2 * len(experiments), 4.6), squeeze=False)
        axes_by_exp = {}
        for ax, exp in zip(axes[0], experiments):
            sub = samples[samples["experiment"] == exp]
            grp = sub.groupby("organism").agg(
                theory=("theoretical_dna_fraction", "first"),
                mean=("measured_base_fraction", "mean"),
                sd=("measured_base_fraction", "std"),
                n=("measured_base_fraction", "size")).reset_index()
            grp["sd"] = grp["sd"].fillna(0.0)

            ax.errorbar(grp["theory"], grp["mean"], yerr=grp["sd"], fmt="o",
                        color="#D55E00", ecolor="#D55E00", elinewidth=0.9,
                        capsize=2.5, markersize=5, markeredgecolor="black",
                        markeredgewidth=0.4, linestyle="none", zorder=3)

            lo = min(grp["theory"].min(), grp["mean"][grp["mean"] > 0].min()
                     if (grp["mean"] > 0).any() else grp["theory"].min())
            hi = max(grp["theory"].max(), grp["mean"].max())
            span = [lo * 0.3, hi * 3]
            ax.plot(span, span, color="0.6", lw=0.8, ls=(0, (4, 3)), zorder=1)

            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlim(*span); ax.set_ylim(*span)
            ax.set_xlabel("Abundance (theoretical, genomic DNA fraction)")
            ax.set_ylabel("Abundance (measured, fraction of sample bases)")
            ax.set_title(f"{exp_names.get(exp, exp)}  "
                         f"(n={sub['sample_id'].nunique()} replicates)")
            ax.grid(True, which="major", color="0.93", lw=0.5, zorder=0)
            ax.grid(True, which="minor", color="0.97", lw=0.4, zorder=0)
            axes_by_exp[exp] = (ax, grp)

        # Labels go on last, once tight_layout has fixed where every axes sits,
        # so the overlap test measures the boxes that will actually be drawn.
        fig.tight_layout()
        for exp, (ax, grp) in axes_by_exp.items():
            shown = grp[grp["mean"] > 0].sort_values("theory")
            to_px = ax.transData.transform
            bars = []
            for _, r in shown.iterrows():           # error bars are obstacles too
                lo_y = max(r["mean"] - r["sd"], ax.get_ylim()[0])
                (x0, y0), (_, y1) = to_px((r["theory"], lo_y)), to_px((r["theory"], r["mean"] + r["sd"]))
                bars.append((x0 - 3, y0, x0 + 3, y1))
            place_labels(ax, list(zip(shown["theory"], shown["mean"])),
                         [italicize(o) for o in shown["organism"]], obstacles=bars)
        fig.savefig(outdir / "abundance.pdf", bbox_inches="tight")
        fig.savefig(outdir / "abundance.png", bbox_inches="tight", dpi=600)
        plt.close(fig)

        grp_all = (samples.groupby(["experiment", "mode", "organism"])
                   .agg(theoretical=("theoretical_dna_fraction", "first"),
                        measured_mean=("measured_base_fraction", "mean"),
                        measured_sd=("measured_base_fraction", "std"),
                        n_replicates=("measured_base_fraction", "size"))
                   .reset_index())
        grp_all.to_csv(outdir / "abundance.csv", index=False)
        write_sidecar(
            outdir / "abundance.json", "abundance",
            "Theoretical versus measured community abundance",
            "Theoretical abundance is the manufacturer's genomic-DNA fraction for each "
            "mock-community member. Measured abundance is that organism's share of the "
            "bases assigned to the community after competitive assignment against a "
            "combined reference containing the lambda carrier, the E. coli K-12 "
            "contaminant carried over from carrier production, and every community "
            "member. Points are the mean across replicates; error bars are the "
            "standard deviation. The dashed line is 1:1.",
            [Path(f).name for f in args.metrics],
            {"n_organisms": int(len(grp_all)),
             "n_plotted_points": int(len(grp_all)),
             "assignment_mode": abundance_mode,
             "n_replicates_by_experiment": {
                 e: int(g["sample_id"].nunique())
                 for e, g in samples.groupby("experiment")}})
        print(f"[aggregate] wrote abundance.pdf/.png/.csv/.json")

    # ---- Figure: read length by role (the ejection signature) -------------
    if args.readlengths:
        # Every mode writes <sample_id>.readlengths.tsv.gz, so a --mode both run
        # offers two files per sample holding the same reads under different
        # labels. Summing both double-counts every read: doing exactly that by
        # hand put carrier n at 105,559,752 instead of 50,513,219. main.nf
        # filters to one mode before staging, but nothing stopped a direct
        # invocation, and the guard belongs with the code that can be misused.
        by_name = defaultdict(list)
        for f in args.readlengths:
            by_name[Path(f).name].append(f)
        dupes = {k: v for k, v in by_name.items() if len(v) > 1}
        if dupes:
            raise SystemExit(
                "error: read lengths supplied more than once for the same "
                "sample, which double-counts every read: "
                + "; ".join(f"{k} x{len(v)}" for k, v in sorted(dupes.items()))
                + ". Pass a single assignment mode.")

        def int_edges(per_decade, top=2_000_000):
            """Log-spaced bin edges placed on half-integers.

            Read lengths are integers. A log bin narrower than a few bp holds 0,
            1 or 2 integer lengths depending only on where its edges fall, which
            draws a comb over the short end of every distribution. Putting edges
            on half-integers makes each bin hold a whole number of lengths, and
            dividing each count by that number (see density()) removes the comb.
            """
            n = int(round(per_decade * np.log10(top))) + 1
            e = np.unique(np.round(np.logspace(0, np.log10(top), n))) - 0.5
            return e.astype(float)

        # FINE: for quantiles (2000 bins, as before). SHOW: what is drawn.
        FINE = int_edges(2000 / np.log10(2_000_000))
        SHOW = int_edges(50)
        sample_exp = dict(zip(per_org["sample_id"], per_org["experiment"]))

        hist = defaultdict(lambda: np.zeros(len(FINE) - 1, dtype=np.int64))
        shown = defaultdict(lambda: np.zeros(len(SHOW) - 1, dtype=np.int64))
        exact = defaultdict(lambda: {"n": 0, "sum": 0, "max": 0})
        # Panel b: carrier lengths per experiment, by read and by base.
        car_reads = defaultdict(lambda: np.zeros(len(SHOW) - 1, dtype=np.int64))
        car_bases = defaultdict(lambda: np.zeros(len(SHOW) - 1, dtype=np.float64))
        car_exact = defaultdict(lambda: {"n": 0, "sum": 0})
        car_fine = defaultdict(lambda: np.zeros(len(FINE) - 1, dtype=np.int64))

        # These files hold one row per read -- ~120 million rows across this
        # study's replicates. Concatenating them into a single DataFrame needs
        # well over 10 GB and gets the process OOM-killed, so accumulate
        # per-role histograms in fixed bins instead and never hold more than one
        # chunk at a time. Counts, sums and maxima stay exact; quantiles are
        # derived from the histogram and are therefore approximate, to within
        # the bin width (under 1%).
        for f in args.readlengths:
            sid = Path(f).name.split(".readlengths")[0]
            exp = sample_exp.get(sid, sid.rsplit("_r", 1)[0])
            for chunk in pd.read_csv(f, sep="\t", compression="gzip",
                                     usecols=["role", "read_length"],
                                     dtype={"role": str, "read_length": np.int64},
                                     chunksize=4_000_000):
                chunk = chunk[chunk["read_length"] > 0]
                for role, sub in chunk.groupby("role", sort=False):
                    v = sub["read_length"].to_numpy()
                    hist[role] += np.histogram(v, bins=FINE)[0]
                    shown[role] += np.histogram(v, bins=SHOW)[0]
                    e = exact[role]
                    e["n"] += int(v.size)
                    e["sum"] += int(v.sum())
                    e["max"] = max(e["max"], int(v.max()))
                    if role == "carrier":
                        car_reads[exp] += np.histogram(v, bins=SHOW)[0]
                        car_bases[exp] += np.histogram(v, bins=SHOW, weights=v)[0]
                        car_fine[exp] += np.histogram(v, bins=FINE)[0]
                        car_exact[exp]["n"] += int(v.size)
                        car_exact[exp]["sum"] += int(v.sum())

        def hist_quantile(counts, q, edges=FINE):
            """Quantile from binned counts, interpolated within the bin."""
            total = counts.sum()
            if total == 0:
                return float("nan")
            cum = np.cumsum(counts)
            i = int(np.searchsorted(cum, q * total))
            i = min(i, len(counts) - 1)
            lo, hi = edges[i], edges[i + 1]
            before = cum[i - 1] if i > 0 else 0
            frac = ((q * total - before) / counts[i]) if counts[i] else 0.0
            return float(lo + frac * (hi - lo))

        SHOW_N_INT = np.diff(SHOW)                       # integer lengths per bin
        SHOW_MID = np.sqrt(np.maximum(SHOW[:-1], 0.5) * SHOW[1:])

        def density(counts, total):
            """Share of the total per drawn bin, free of the integer comb.

            For a bin wide enough to hold many lengths this is simply
            counts/total, the fraction per bin. For a narrow bin it is the
            per-length density scaled to the bin's nominal log width, so every
            bin is on the same footing however many integers it happens to hold.
            """
            nominal = SHOW_MID * np.log(10) / 50          # bp spanned by 1/50 decade
            return counts / total / SHOW_N_INT * nominal

        # assign_reads.py writes role="none" for reads no reference aligned, so
        # "none" must be listed or unmapped reads silently vanish from the figure.
        roles = [r for r in ["carrier", "contaminant", "sample", "ambiguous",
                             "none", "unassigned"] if exact.get(r, {}).get("n")]
        exps = sorted(car_exact)
        LAMBDA_BP = 48_502

        fig = plt.figure(figsize=(7.2, 6.6))
        gs = fig.add_gridspec(2, max(1, len(exps)), height_ratios=[1.25, 1],
                              hspace=0.42, wspace=0.08)
        ax = fig.add_subplot(gs[0, :])
        plotted = []
        for role in roles:
            n = exact[role]["n"]
            y = density(shown[role], n)
            med = hist_quantile(hist[role], 0.5)
            ax.stairs(y, SHOW, color=ROLE_COLORS.get(role, "0.5"), linewidth=1.3,
                      label=f"{role} (n={n:,}, median {med:,.0f} bp)")
            for i in np.nonzero(shown[role])[0]:
                plotted.append({"panel": "a", "series": role, "weighting": "reads",
                                "bin_lower_bp": SHOW[i], "bin_upper_bp": SHOW[i + 1],
                                "count": int(shown[role][i]), "plotted_value": y[i]})
        ax.set_xscale("log")
        # Reads below ~20 bp are adapter-length artefacts, not molecules; keep
        # them out of the view so the axis is set by the real distributions.
        ax.set_xlim(20, max(1e4, max(exact[r]["max"] for r in roles) * 1.1))
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Read length (bp)")
        ax.set_ylabel("Fraction of reads per bin")
        ax.legend(frameon=True, framealpha=0.95, edgecolor="0.7", fontsize=7.5)
        ax.grid(True, which="major", color="0.93", lw=0.5, zorder=0)
        ax.text(-0.09, 1.02, "a", transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom")

        # Panel b: carrier reads only, by read and by base. Weighting by bases
        # shows where the carrier SEQUENCE is: the few intact lambda molecules
        # that escaped ejection (the 48.5 kb mode) hold a share of carrier bases
        # far larger than their share of reads.
        c = ROLE_COLORS["carrier"]
        bax = None
        car_stats = {}
        for j, exp in enumerate(exps):
            bax = fig.add_subplot(gs[1, j], sharey=bax)
            ce = car_exact[exp]
            yr = density(car_reads[exp], ce["n"])
            yb = density(car_bases[exp], ce["sum"])
            bax.stairs(yb, SHOW, fill=True, color=c, alpha=0.30, lw=0,
                       label="of carrier bases")
            bax.stairs(yr, SHOW, color=c, lw=1.3, label="of carrier reads")
            bax.axvline(LAMBDA_BP, color="0.35", lw=0.8, ls=(0, (3, 3)), zorder=1)
            bax.set_xscale("log")
            bax.set_xlim(20, 2e5)
            bax.set_ylim(bottom=0)
            med = hist_quantile(car_fine[exp], 0.5)
            intact = car_bases[exp][SHOW[:-1] >= 40_000].sum() / ce["sum"]
            intact_r = car_reads[exp][SHOW[:-1] >= 40_000].sum() / ce["n"]
            car_stats[exp] = {"n": ce["n"], "median_bp": med,
                              "share_of_bases_at_or_above_40kb": float(intact),
                              "share_of_reads_at_or_above_40kb": float(intact_r)}
            bax.set_title(f"{exp_names.get(exp, exp)} carrier reads\n"
                          f"n={ce['n']:,}, median {med:,.0f} bp", fontsize=8.5)
            bax.set_xlabel("Read length (bp)")
            bax.grid(True, which="major", color="0.93", lw=0.5, zorder=0)
            if j == 0:
                bax.set_ylabel("Fraction per bin")
                bax.text(-0.2, 1.06, "b", transform=bax.transAxes, fontsize=11,
                         fontweight="bold", va="bottom")
            else:
                plt.setp(bax.get_yticklabels(), visible=False)
            for weighting, arr, y in (("reads", car_reads[exp], yr),
                                      ("bases", car_bases[exp], yb)):
                for i in np.nonzero(arr)[0]:
                    plotted.append({"panel": "b", "series": f"carrier:{exp}",
                                    "weighting": weighting,
                                    "bin_lower_bp": SHOW[i], "bin_upper_bp": SHOW[i + 1],
                                    "count": float(arr[i]), "plotted_value": y[i]})
        if bax is not None:
            ymax = bax.get_ylim()[1]
            bax.text(LAMBDA_BP * 0.75, ymax * 0.95, "λ genome\n48.5 kb", ha="right",
                     va="top", fontsize=7, color="0.3")
            bax.legend(frameon=False, fontsize=7, loc="upper left")

        fig.savefig(outdir / "readlengths.pdf", bbox_inches="tight")
        fig.savefig(outdir / "readlengths.png", bbox_inches="tight", dpi=600)
        plt.close(fig)

        stats = pd.DataFrame([{
            "role": role,
            "n": exact[role]["n"],
            "mean": exact[role]["sum"] / exact[role]["n"],
            "median": hist_quantile(hist[role], 0.50),
            "p25": hist_quantile(hist[role], 0.25),
            "p75": hist_quantile(hist[role], 0.75),
            "max": exact[role]["max"],
        } for role in roles])

        # The CSV that accompanies a figure must hold the points the figure
        # actually draws, so a reader can redraw it without rerunning anything.
        # Summary statistics alone are not that; they go in the JSON sidecar.
        pd.DataFrame(plotted).to_csv(outdir / "readlengths.csv", index=False)
        stats.to_csv(outdir / "readlengths_summary.csv", index=False)
        write_sidecar(
            outdir / "readlengths.json", "readlengths",
            "Read-length distribution by assignment class, and of carrier reads",
            "(a) Read length by assignment class, all replicates pooled. "
            "Depletion-mode adaptive sampling ejects a molecule once it is recognised "
            "as carrier, truncating that read, while molecules that are not rejected "
            "are sequenced to full length, so carrier reads are expected to be "
            "systematically shorter than community reads. Distributions are normalised "
            "per class. (b) Carrier (lambda) reads per experiment, as a fraction of "
            "carrier reads (line) and of carrier bases (shaded). The dashed line marks "
            "the 48,502 bp lambda genome: intact molecules that escaped ejection are a "
            "small share of carrier reads but a large share of carrier bases. "
            "Bins are log-spaced, 50 per decade, with edges on half-integers; each "
            "count is scaled by the number of integer lengths its bin spans relative "
            "to the bin's nominal width, which removes the comb that integer read "
            "lengths otherwise draw in narrow bins. Counts, means and maxima are "
            "exact; medians and quartiles come from 2000 log-spaced bins over 1 bp "
            "to 2 Mb (accurate to well under 1%), because the per-read tables are too "
            "large to hold in memory at once.",
            [Path(f).name for f in args.readlengths],
            {**{r["role"]: {"n": int(r["n"]), "median_bp": float(r["median"])}
                for _, r in stats.iterrows()},
             "carrier_by_experiment": {exp_names.get(e, e): v for e, v in car_stats.items()},
             "n_plotted_points": int(len(plotted))})
        print(f"[aggregate] wrote readlengths.pdf/.png/.csv/.json")
        print(stats.to_string(index=False))
        for e, v in car_stats.items():
            print(f"[aggregate] carrier {exp_names.get(e, e)}: {v}")

    print(f"[aggregate] outputs in {outdir}")


if __name__ == "__main__":
    main()
