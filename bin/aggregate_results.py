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

Replicates flagged include_in_headline=0 in the samplesheet (e.g. lowinput_s2_r0,
whose input mass was below Qubit quantification) are carried in the per-sample
table but excluded from the experiment-level statistics, because every per-fg
metric divides by an input mass those samples do not have.
"""
from __future__ import annotations

import argparse
import json
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
    p.add_argument("--samplesheet", default=None,
                   help="samplesheet CSV, for include_in_headline flags")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_samplesheet(path):
    if not path:
        return {}
    lines = [l for l in Path(path).read_text().splitlines()
             if l.strip() and not l.strip().startswith("#")]
    if not lines:
        return {}
    import csv
    import io
    rows = list(csv.DictReader(io.StringIO("\n".join(lines))))
    return {r["sample_id"]: (str(r.get("include_in_headline", "1")).strip() == "1")
            for r in rows}


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

    headline_flags = read_samplesheet(args.samplesheet)

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

    if len(samples):
        experiments = sorted(samples["experiment"].unique())
        fig, axes = plt.subplots(1, len(experiments),
                                 figsize=(5.2 * len(experiments), 4.6), squeeze=False)
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

            # Alternate the label side by rank so neighbouring points in a
            # log-distributed community (which crowd along the diagonal) do not
            # overprint each other.
            ranked = grp[grp["mean"] > 0].sort_values("theory").reset_index(drop=True)
            for i, r in ranked.iterrows():
                right = (i % 2 == 0)
                ax.annotate(italicize(r["organism"]), (r["theory"], r["mean"]),
                            textcoords="offset points",
                            xytext=(6, 5) if right else (-6, -11),
                            ha="left" if right else "right",
                            fontsize=6.5, color="0.25")

            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlim(*span); ax.set_ylim(*span)
            ax.set_xlabel("Abundance (theoretical, genomic DNA fraction)")
            ax.set_ylabel("Abundance (measured, fraction of sample bases)")
            ax.set_title(f"{exp}  (n={sub['sample_id'].nunique()} replicates)")
            ax.grid(True, which="major", color="0.93", lw=0.5, zorder=0)
            ax.grid(True, which="minor", color="0.97", lw=0.4, zorder=0)

        fig.tight_layout()
        fig.savefig(outdir / "abundance.pdf", bbox_inches="tight")
        fig.savefig(outdir / "abundance.png", bbox_inches="tight", dpi=600)
        plt.close(fig)

        grp_all = (samples.groupby(["experiment", "organism"])
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
             "n_plotted_points": int(len(grp_all))})
        print(f"[aggregate] wrote abundance.pdf/.png/.csv/.json")

    # ---- Figure: read length by role (the ejection signature) -------------
    if args.readlengths:
        # These files hold one row per read -- ~120 million rows across this
        # study's replicates. Concatenating them into a single DataFrame needs
        # well over 10 GB and gets the process OOM-killed, so accumulate
        # per-role histograms in fixed bins instead and never hold more than one
        # chunk at a time. Counts, sums and maxima stay exact; quantiles are
        # derived from the histogram and are therefore approximate, to within
        # the bin width (2000 log-spaced bins over 1 bp - 2 Mb, so well under 1%).
        BINS = np.logspace(0, np.log10(2_000_000), 2001)
        widths = np.diff(BINS)
        hist = defaultdict(lambda: np.zeros(len(BINS) - 1, dtype=np.int64))
        exact = defaultdict(lambda: {"n": 0, "sum": 0, "max": 0})

        for f in args.readlengths:
            for chunk in pd.read_csv(f, sep="\t", compression="gzip",
                                     usecols=["role", "read_length"],
                                     dtype={"role": str, "read_length": np.int64},
                                     chunksize=4_000_000):
                chunk = chunk[chunk["read_length"] > 0]
                for role, sub in chunk.groupby("role", sort=False):
                    v = sub["read_length"].to_numpy()
                    hist[role] += np.histogram(v, bins=BINS)[0]
                    e = exact[role]
                    e["n"] += int(v.size)
                    e["sum"] += int(v.sum())
                    e["max"] = max(e["max"], int(v.max()))

        def hist_quantile(counts, q):
            """Quantile from binned counts, interpolated within the bin."""
            total = counts.sum()
            if total == 0:
                return float("nan")
            cum = np.cumsum(counts)
            i = int(np.searchsorted(cum, q * total))
            i = min(i, len(counts) - 1)
            lo, hi = BINS[i], BINS[i + 1]
            before = cum[i - 1] if i > 0 else 0
            frac = ((q * total - before) / counts[i]) if counts[i] else 0.0
            return float(lo + frac * (hi - lo))

        # assign_reads.py writes role="none" for reads no reference aligned, so
        # "none" must be listed or unmapped reads silently vanish from the figure.
        roles = [r for r in ["carrier", "contaminant", "sample", "ambiguous",
                             "none", "unassigned"] if exact.get(r, {}).get("n")]

        fig, ax = plt.subplots(figsize=(6.8, 4.4))
        for role in roles:
            counts = hist[role]
            n = exact[role]["n"]
            med = hist_quantile(counts, 0.5)
            # Plot the fraction of reads per bin, NOT counts/(n*width).
            # The bins are equal-width in log10, so counts/n already is a
            # density with respect to log(length) -- which is the right density
            # for a log x axis. Dividing by the linear bin width instead makes
            # narrow low-length bins explode: a spike of ~10 bp unmapped reads
            # then dominates the axis and flattens every distribution of
            # interest into the baseline.
            ax.stairs(counts / n, BINS, color=ROLE_COLORS.get(role, "0.5"),
                      linewidth=1.3,
                      label=f"{role} (n={n:,}, median {med:,.0f} bp)")
        ax.set_xscale("log")
        # Reads below ~20 bp are adapter-length artefacts, not molecules; keep
        # them out of the view so the axis is set by the real distributions.
        ax.set_xlim(20, max(1e4, max(exact[r]["max"] for r in roles) * 1.1))
        ax.set_xlabel("Read length (bp)")
        ax.set_ylabel("Fraction of reads per bin")
        ax.set_title("Read length by assignment class")
        ax.legend(frameon=True, framealpha=0.95, edgecolor="0.7", fontsize=7.5)
        ax.grid(True, which="major", color="0.93", lw=0.5, zorder=0)
        fig.tight_layout()
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
        plotted = []
        for role in roles:
            counts = hist[role]
            n = exact[role]["n"]
            for i, c in enumerate(counts):
                if c == 0:
                    continue
                plotted.append({
                    "role": role,
                    "bin_lower_bp": BINS[i],
                    "bin_upper_bp": BINS[i + 1],
                    "count": int(c),
                    "fraction_of_reads": c / n,
                })
        pd.DataFrame(plotted).to_csv(outdir / "readlengths.csv", index=False)
        stats.to_csv(outdir / "readlengths_summary.csv", index=False)
        write_sidecar(
            outdir / "readlengths.json", "readlengths",
            "Read-length distribution by assignment class",
            "Depletion-mode adaptive sampling ejects a molecule once it is recognised "
            "as carrier, truncating that read, while molecules that are not rejected "
            "are sequenced to full length. Carrier reads are therefore expected to be "
            "systematically shorter than community reads; this figure tests that "
            "expectation directly. Densities are normalised per class. Counts, means "
            "and maxima are exact; medians and quartiles are derived from 2000 "
            "log-spaced bins spanning 1 bp to 2 Mb (accurate to well under 1%), "
            "because the per-read tables are too large to hold in memory at once.",
            [Path(f).name for f in args.readlengths],
            {**{r["role"]: {"n": int(r["n"]), "median_bp": float(r["median"])}
                for _, r in stats.iterrows()},
             "n_plotted_points": int(len(plotted))})
        print(f"[aggregate] wrote readlengths.pdf/.png/.csv/.json")
        print(stats.to_string(index=False))

    print(f"[aggregate] outputs in {outdir}")


if __name__ == "__main__":
    main()
