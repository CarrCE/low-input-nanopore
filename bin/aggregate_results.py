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
import gzip
import json
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

            for _, r in grp.iterrows():
                if r["mean"] > 0:
                    ax.annotate(italicize(r["organism"]), (r["theory"], r["mean"]),
                                textcoords="offset points", xytext=(5, 4),
                                fontsize=7, color="0.25")

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
            {"n_organisms": int(len(grp_all))})
        print(f"[aggregate] wrote abundance.pdf/.png/.csv/.json")

    # ---- Figure: read length by role (the ejection signature) -------------
    if args.readlengths:
        frames = []
        for f in args.readlengths:
            sid = Path(f).name.split(".")[0]
            with gzip.open(f, "rt") as fh:
                df = pd.read_csv(fh, sep="\t")
            df["sample_id"] = sid
            frames.append(df[["sample_id", "role", "read_length"]])
        rl = pd.concat(frames, ignore_index=True)

        roles = [r for r in ["carrier", "contaminant", "sample", "ambiguous", "unassigned"]
                 if (rl["role"] == r).any()]
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        bins = np.logspace(np.log10(max(rl["read_length"].min(), 10)),
                           np.log10(max(rl["read_length"].max(), 100)), 70)
        for role in roles:
            v = rl.loc[rl["role"] == role, "read_length"]
            v = v[v > 0]
            if not len(v):
                continue
            ax.hist(v, bins=bins, histtype="step", density=True, linewidth=1.3,
                    color=ROLE_COLORS.get(role, "0.5"),
                    label=f"{role} (n={len(v):,}, median {int(v.median()):,} bp)")
        ax.set_xscale("log")
        ax.set_xlabel("Read length (bp)")
        ax.set_ylabel("Density")
        ax.set_title("Read length by assignment class")
        ax.legend(frameon=True, framealpha=0.95, edgecolor="0.7", fontsize=7.5)
        ax.grid(True, which="major", color="0.93", lw=0.5, zorder=0)
        fig.tight_layout()
        fig.savefig(outdir / "readlengths.pdf", bbox_inches="tight")
        fig.savefig(outdir / "readlengths.png", bbox_inches="tight", dpi=600)
        plt.close(fig)

        stats = (rl.groupby("role")["read_length"]
                   .agg(n="size", mean="mean", median="median",
                        p25=lambda s: s.quantile(0.25),
                        p75=lambda s: s.quantile(0.75), max="max")
                   .reset_index())
        stats.to_csv(outdir / "readlengths.csv", index=False)
        write_sidecar(
            outdir / "readlengths.json", "readlengths",
            "Read-length distribution by assignment class",
            "Depletion-mode adaptive sampling ejects a molecule once it is recognised "
            "as carrier, truncating that read, while molecules that are not rejected "
            "are sequenced to full length. Carrier reads are therefore expected to be "
            "systematically shorter than community reads; this figure tests that "
            "expectation directly. Densities are normalised per class.",
            [Path(f).name for f in args.readlengths],
            {r["role"]: {"n": int(r["n"]), "median_bp": float(r["median"])}
             for _, r in stats.iterrows()})
        print(f"[aggregate] wrote readlengths.pdf/.png/.csv/.json")
        print(stats.to_string(index=False))

    print(f"[aggregate] outputs in {outdir}")


if __name__ == "__main__":
    main()
