#!/usr/bin/env python3
"""Compare the carrier-depletion model with the deposited reads.

Reads <workdir>/<s>.reads.npz (extract_reads.py), <workdir>/<s>.instrument.npz
(join_instrument.py), <outdir>/cos_ends.tsv (cos_ends.py),
assets/measurements.tsv, and the pipeline's results/summary/per_sample.tsv and
per_organism.tsv.

Writes to <outdir>/ (analysis tables):
  decomposition.tsv    per replicate: read classes, lambda-derived read length
                       structure, observed enrichment, every model's prediction
  time_resolved.tsv    per replicate x 15-min bin: yield, sample share, escape
                       fraction, length of reads the instrument unblocked
  per_organism_s1.tsv  D6311 organisms: bases per fg against mean read length

and to <figdir>/ (display items, docs/display-items.md):
  carrier-model.{pdf,png,csv,json}         Supplementary Figure S4
  carrier-model-table.{csv,json}           Supplementary Table S15

Every model quantity counts lambda-derived reads (model.py). The paper's carrier
class (lambda-assigned only) appears only in the read-class accounting columns
and as the 465 bp median that the simple model assumes.

Usage: analyze.py <repo_root> <results> <workdir> <outdir> <figdir>
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model as cm  # noqa: E402
from common import (BIN_S, ESCAPE_BP, SAMPLES, classes, display, read_tsv,  # noqa: E402
                    write_tsv)

import matplotlib  # noqa: E402
import matplotlib.lines  # noqa: E402,F401
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
})
# Panel b colours the two experiments (Okabe-Ito blue and vermillion, as in the
# other figures). Panel a separates its two model series by marker and grey
# level instead, so that no colour means two things in one figure.
C_EXP = {"lowinput_s1": "#0072B2", "lowinput_s2": "#D55E00"}
C_SIMPLE, C_INTACT = "0.62", "0.12"
LS_REP = {"r0": (0, (5, 1.5, 1, 1.5)), "r1": "-", "r2": (0, (4, 2)), "r3": (0, (1, 1.4))}
EJECT_WINDOW = (450, 2000)
MIN_BIN_READS = 1000
FIG_ID, TABLE_ID = "carrier-model", "carrier-model-table"


def main(repo, ror, work, out, figdir):
    repo, ror, work, out, figdir = map(Path, (repo, ror, work, out, figdir))
    out.mkdir(parents=True, exist_ok=True); figdir.mkdir(parents=True, exist_ok=True)
    meas = {r["sample_id"]: r for r in read_tsv(repo / "assets/measurements.tsv")}
    per = {r["sample_id"]: r for r in read_tsv(ror / "summary/per_sample.tsv")
           if r["mode"] == "competitive" and r["role"] == "sample_total"}
    cosr = read_tsv(out / "cos_ends.tsv")
    cos = {r["sample_id"]: r for r in cosr if r["read_set"] == "hash"}
    esc_cos = {r["sample_id"]: r for r in cosr if r["read_set"] == "escaped"}

    rows, time_rows = [], []
    for s in SAMPLES:
        d = np.load(work / f"{s}.reads.npz")
        ins = np.load(work / f"{s}.instrument.npz")
        act = ins["action"]
        L = d["length"].astype(np.int64)
        assert len(act) == len(L), f"{s}: instrument join is not in deposit order"
        # Run time from the sequencing summary, not the FASTQ header tags: an
        # archive-fetched FASTQ has none. Every deposited read is in the summary.
        assert ins["in_summary"].all(), f"{s}: deposited reads missing from the summary"
        st = ins["start_time"]
        R = classes(d)
        lam, car, smp = R["lambda_derived"], R["carrier"], R["sample"]

        # consistency with the pipeline's per-sample table
        assert len(L) == int(per[s]["total_reads"]), s
        assert L.sum() == int(per[s]["total_bases"]), s
        assert smp.sum() == int(per[s]["reads"]) and L[smp].sum() == int(per[s]["bases"]), s

        m_s = float(meas[s]["sample_dna_ng"]); m_c = float(meas[s]["carrier_dna_ng"])
        Ll = L[lam]
        obs = dict(m_s=m_s, m_c=m_c, N_s=smp.sum(), B_s=L[smp].sum(),
                   N_c=lam.sum(), N_tot=len(L), B_tot=L.sum())
        obs_class = dict(obs, N_c=car.sum())
        E_b, E_r = float(per[s]["enrichment_bases"]), float(per[s]["enrichment_reads"])
        b_s = L[smp].mean()

        ej = Ll[(Ll >= EJECT_WINDOW[0]) & (Ll < EJECT_WINDOW[1])]
        counts, edges = np.histogram(Ll[Ll < EJECT_WINDOW[1]], bins=np.arange(100, 1201, 10))
        unb = Ll[act[lam] == 1]
        esc = Ll >= ESCAPE_BP
        f_cos = float(cos[s]["f_cos"])
        full = float(esc_cos[s]["full_length_frac"])
        si_b, si_r = cm.simple_enrichment(m_s, b_s, m_c)
        p_int = cm.predict(obs, cm.LAMBDA_BP)
        p_cos = cm.predict(obs, cm.LAMBDA_BP * f_cos)      # phi = 1: every cos end adaptable
        L_eff = cm.implied_L_eff(obs)
        phi, beta = cm.decompose(L_eff, f_cos)
        # early vs late in the run: first 6 h and last 12 h from the first deposited read
        t_h = (st - np.nanmin(st)) / 3600
        win = {}
        for tag, (a, z) in (("first6h", (0, 6)), ("last12h", (np.nanmax(t_h) - 12, np.inf))):
            w = (t_h >= a) & (t_h < z)
            lw = lam & w
            win[f"escape_frac_{tag}"] = (L[lw] >= ESCAPE_BP).mean()
            win[f"kept_frac_{tag}"] = (act[lw] == 2).mean()
            win[f"unblocked_median_bp_{tag}"] = float(np.median(L[lw & (act == 1)]))
            win[f"sample_per_lambda_derived_{tag}"] = (smp & w).sum() / lw.sum()
        rows.append({
            "sample_id": s, "library": display(s),
            "minknow": meas[s]["minknow_version"], "lambda_lot": meas[s]["carrier_lot"],
            "include_in_headline": meas[s]["include_in_headline"],
            "m_s_ng": m_s, "m_c_ng": m_c,
            "reads": len(L), "sample_reads": int(smp.sum()),
            "carrier_class_reads": int(car.sum()),
            "ambiguous_lambda_reads": int(R["ambiguous_lambda"].sum()),
            "lambda_derived_reads": int(lam.sum()),
            "contaminant_reads": int(R["contaminant"].sum()),
            "ambiguous_reads": int(R["ambiguous"].sum()),
            "unassigned_reads": int(R["none"].sum()),
            "unassigned_bases_frac": L[R["none"]].sum() / L.sum(),
            "sample_mean_len": b_s, "all_mean_len": L.mean(),
            "carrier_class_mean_len": L[car].mean(),
            "carrier_class_median_len": float(np.median(L[car])),
            "lambda_mean_len": Ll.mean(), "lambda_median_len": float(np.median(Ll)),
            "lambda_lt450_frac": (Ll < 450).mean(),
            "lambda_eject_window_frac": len(ej) / len(Ll),
            "eject_mode_bp": float(edges[np.argmax(counts)] + 5),
            "eject_window_median_bp": float(np.median(ej)),
            "unblocked_frac": len(unb) / len(Ll),
            "unblocked_median_bp": float(np.median(unb)),
            "escape_frac_reads": esc.mean(), "escape_frac_lambda_bases": Ll[esc].sum() / Ll.sum(),
            "escaped_mean_len": Ll[esc].mean(),
            "lambda_reads_ge40k": int((Ll >= 40000).sum()),
            "E_bases_obs": E_b, "E_reads_obs": E_r,
            "E_bases_simple": si_b, "E_reads_simple": si_r,
            "L_carrier_simple_inverted": E_b * cm.B_EJECT_SIMPLE,
            "E_bases_intact": p_int["E_bases"], "E_reads_intact": p_int["E_reads"],
            "f_cos": f_cos, "f_cos_lo": float(cos[s]["f_cos_lo"]), "f_cos_hi": float(cos[s]["f_cos_hi"]),
            "E_bases_cos": p_cos["E_bases"], "E_reads_cos": p_cos["E_reads"],
            "L_eff_implied": L_eff,
            "L_eff_implied_carrier_class": cm.implied_L_eff(obs_class),
            "phi_cos_adaptable": phi, "beta_breaks_per_genome": beta,
            "escaped_cos_full_length_frac": full,
            "beta_upper_bound": -math.log(full) if full > 0 else math.inf,
            "m_s_implied_intact_ng": cm.implied_m_s(obs, cm.LAMBDA_BP),
            **win,
        })

        # time-resolved, 15-min bins from the first deposited read
        t = st - np.nanmin(st)
        b = (t // BIN_S).astype(int)
        for k in range(b.max() + 1):
            m = b == k
            if m.sum() < MIN_BIN_READS:
                continue
            lm = lam & m
            ll = L[lm]
            lu = L[lm & (act == 1)]
            time_rows.append({
                "sample_id": s, "t_start_h": k * BIN_S / 3600, "t_mid_h": (k + 0.5) * BIN_S / 3600,
                "reads": int(m.sum()), "bases": int(L[m].sum()),
                "sample_read_frac": smp[m].mean(), "sample_base_frac": L[smp & m].sum() / L[m].sum(),
                "sample_per_lambda_derived": (smp & m).sum() / max(lm.sum(), 1),
                "lambda_escape_frac": (ll >= ESCAPE_BP).mean() if len(ll) else math.nan,
                "lambda_mean_len": ll.mean() if len(ll) else math.nan,
                "lambda_unblocked_frac": len(lu) / len(ll) if len(ll) else math.nan,
                "unblocked_median_bp": float(np.median(lu)) if len(lu) else math.nan,
                "eject_window_median_bp": float(np.median(ll[(ll >= 450) & (ll < 2000)]))
                if len(ll) else math.nan,
            })

    write_tsv(out / "decomposition.tsv", rows)
    write_tsv(out / "time_resolved.tsv", time_rows)
    slopes = per_organism(ror, out)
    figure(rows, time_rows, figdir, slopes)
    table(rows, figdir)
    for r in rows:
        print(f"[carrier-model] {r['library']:9} obs {r['E_bases_obs']:5.1f}  simple {r['E_bases_simple']:5.1f}  "
              f"intact {r['E_bases_intact']:5.1f}  f_cos {r['f_cos']:.3f} L_eff {r['L_eff_implied']:7.0f} "
              f"phi {r['phi_cos_adaptable']:.2f} beta {r['beta_breaks_per_genome']:.2f} "
              f"(<= {r['beta_upper_bound']:.2f})  escape {r['escape_frac_reads']:.4f}/"
              f"{r['escape_frac_lambda_bases']:.3f}  mean {r['lambda_mean_len']:.0f}")


def per_organism(ror, out):
    """A1/A2 check from natural variation. Under per-molecule capture with sample
    reads whole, bases per fg do not depend on an organism's molecule length
    (slope 0) and reads per fg scale as 1/length (slope -1). E. coli is excluded:
    its D6311 row is dominated by carrier-derived K-12 reads
    (docs/ecoli-partitioning.md)."""
    po = [r for r in read_tsv(ror / "summary/per_organism.tsv")
          if r["mode"] == "competitive" and r["role"] == "sample" and r["experiment"] == "lowinput_s1"]
    rows = []
    for r in po:
        reads = int(r["reads"])
        if reads < 100 or r["organism"].startswith("Escherichia coli"):
            continue
        rows.append({"sample_id": r["sample_id"], "organism": r["organism"], "reads": reads,
                     "mean_read_len": int(r["bases"]) / reads, "bases_per_fg": float(r["bases_per_fg"]),
                     "reads_per_fg": float(r["reads_per_fg"])})
    write_tsv(out / "per_organism_s1.tsv", rows)
    x = np.log([r["mean_read_len"] for r in rows])
    sid = [r["sample_id"] for r in rows]
    slopes = {"n_organism_replicates": len(rows)}
    for key, want in (("bases_per_fg", 0), ("reads_per_fg", -1)):
        y = np.log([r[key] for r in rows])
        # replicate fixed effect: demean within replicate, since masses differ between them
        xd, yd = x.copy(), y.copy()
        for s in set(sid):
            m = np.array([v == s for v in sid]); xd[m] -= x[m].mean(); yd[m] -= y[m].mean()
        b = float(np.polyfit(xd, yd, 1)[0])
        slopes[f"slope_log_{key}_vs_log_mean_read_len"] = b
        print(f"[carrier-model] D6311 per organism: log({key}) vs log(mean read length), "
              f"within-replicate slope {b:.2f} (model predicts {want}) over {len(rows)} organism-replicates")
    return slopes


def figure(rows, time_rows, figdir, slopes):
    fig = plt.figure(figsize=(7.2, 4.4))
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 1.3], wspace=0.32, hspace=0.18,
                          left=0.08, right=0.985, top=0.88, bottom=0.1)
    ax = fig.add_subplot(gs[:, 0])
    bx = [fig.add_subplot(gs[0, 1])]
    bx += [fig.add_subplot(gs[i, 1], sharex=bx[0]) for i in (1, 2)]
    pts = []

    # (a) predicted vs observed
    obs = np.array([r["E_bases_obs"] for r in rows])
    inc = np.array([r["include_in_headline"] == "1" for r in rows])
    hi = max(obs.max(), max(r["E_bases_simple"] for r in rows)) * 1.12
    ax.plot([0, hi], [0, hi], color="0.6", lw=0.8, ls=(0, (4, 3)), zorder=1)
    for key, c, mk, lab in [("E_bases_simple", C_SIMPLE, "s", "465 bp per carrier capture"),
                            ("E_bases_intact", C_INTACT, "o", "Intact λ, measured bases\nper λ-derived read")]:
        y = np.array([r[key] for r in rows])
        ax.scatter(obs[inc], y[inc], s=34, marker=mk, color=c, edgecolor="white",
                   linewidth=1.0, zorder=3, label=lab)
        ax.scatter(obs[~inc], y[~inc], s=34, marker=mk, facecolor="white", edgecolor=c,
                   linewidth=1.1, zorder=3)
        for r, yy in zip(rows, y):
            pts.append({"panel": "a", "series": key, "sample_id": r["sample_id"], "library": r["library"],
                        "x": r["E_bases_obs"], "y": yy, "hollow": r["include_in_headline"] != "1"})
    # (dx, dy, ha) per label, beside its intact-lambda point
    offsets = {"lowinput_s1_r1": (6, 0, "left"), "lowinput_s1_r2": (6, 3, "left"),
               "lowinput_s1_r3": (-7, -3, "right"), "lowinput_s2_r0": (0, -14, "center"),
               "lowinput_s2_r1": (6, -6, "left"), "lowinput_s2_r2": (6, -3, "left"),
               "lowinput_s2_r3": (-7, 4, "right")}
    for r in rows:
        dx, dy, ha = offsets.get(r["sample_id"], (5, 0, "left"))
        ax.annotate(r["library"] + ("\n(excluded)" if r["include_in_headline"] != "1" else ""),
                    (r["E_bases_obs"], r["E_bases_intact"]), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, va="center", fontsize=6.5, color="0.3")
    ax.set_xlim(0, hi); ax.set_ylim(0, hi); ax.set_aspect("equal")
    ax.set_xlabel("Observed enrichment (bases)"); ax.set_ylabel("Predicted enrichment (bases)")
    ax.grid(True, color="0.93", lw=0.5, zorder=0)
    # lower right is the one empty region: every point sits near or above 1:1
    ax.legend(frameon=False, fontsize=6.5, loc="lower right", handletextpad=0.3)

    # (b) over each run
    series = [("sample_per_lambda_derived", "Community reads\nper λ-derived read"),
              ("lambda_escape_frac", "λ-derived reads\n≥2 kb"),
              ("unblocked_median_bp", "Median unblocked\nread (bp)")]
    for axb, (key, lab) in zip(bx, series):
        for s in SAMPLES:
            tr = [r for r in time_rows if r["sample_id"] == s]
            rep = s.rsplit("_", 1)[1]
            axb.plot([r["t_mid_h"] for r in tr], [r[key] for r in tr], color=C_EXP[s[:11]], lw=0.9,
                     ls=LS_REP[rep], label=display(s))
            for r in tr:
                pts.append({"panel": "b", "series": key, "sample_id": s, "library": display(s),
                            "x": r["t_mid_h"], "y": r[key], "hollow": False})
        axb.set_ylabel(lab, fontsize=7)
        axb.grid(True, color="0.93", lw=0.5)
        axb.tick_params(labelsize=7)
    for axb in bx[:-1]:
        plt.setp(axb.get_xticklabels(), visible=False)
    # one row per experiment: legend columns fill top-down, so interleave and pad
    h, lab = bx[0].get_legend_handles_labels()
    blank = matplotlib.lines.Line2D([], [], color="none")
    order = [(h[0], lab[0]), (h[3], lab[3]), (h[1], lab[1]), (h[4], lab[4]),
             (h[2], lab[2]), (h[5], lab[5]), (blank, ""), (h[6], lab[6])]
    bx[0].legend(*zip(*order), frameon=False, fontsize=6.5, ncol=4, loc="lower left",
                 bbox_to_anchor=(-0.02, 1.02), handlelength=2.6, columnspacing=1.0)
    bx[-1].set_xlabel("Hours from first deposited read")
    for letter, x in (("a", 0.005), ("b", 0.47)):
        fig.text(x, 0.975, letter, fontsize=10, fontweight="bold", ha="left", va="top")
    fig.savefig(figdir / f"{FIG_ID}.pdf")
    fig.savefig(figdir / f"{FIG_ID}.png", dpi=600)
    plt.close(fig)

    write_csv(figdir / f"{FIG_ID}.csv", pts)
    head = [r for r in rows if r["include_in_headline"] == "1"]
    rng = lambda key, sel=head: [min(r[key] for r in sel), max(r[key] for r in sel)]
    s1 = [r for r in head if r["sample_id"].startswith("lowinput_s1")]
    s2 = [r for r in head if r["sample_id"].startswith("lowinput_s2")]
    ratio = lambda key, sel: [min(r[key] / r["E_bases_obs"] for r in sel),
                              max(r[key] / r["E_bases_obs"] for r in sel)]
    payload = {
        "id": FIG_ID,
        "title": "Carrier-depletion model against the seven runs",
        "caption": (
            "(a) Predicted versus observed base-weighted enrichment, for 465 bp per carrier "
            "capture (squares) and for intact lambda with the measured mean bases per "
            "lambda-derived read (circles). The dashed line is 1:1. D6321 r0 is drawn hollow "
            "and is excluded from headline statistics. (b) Change over each run: community "
            "reads per lambda-derived read (top), share of lambda-derived reads of 2 kb or "
            "longer (middle), and median length of the lambda-derived reads the instrument "
            "unblocked, from its adaptive-sampling decision records (bottom), in 15-min bins "
            "of time from the first deposited read. Lambda-derived reads are those assigned "
            "to lambda plus those tied between lambda and another organism."),
        "source_files": ["decomposition.tsv", "time_resolved.tsv", "cos_ends.tsv",
                         "<sample>.reads.npz", "<sample>.instrument.npz",
                         "assets/measurements.tsv", "summary/per_sample.tsv"],
        "software": software(),
        "metrics": {
            "n_plotted_points": len(pts),
            "headline_excludes": [r["library"] for r in rows if r["include_in_headline"] != "1"],
            "D6311_simple_over_observed": ratio("E_bases_simple", s1),
            "D6311_intact_over_observed": ratio("E_bases_intact", s1),
            "D6321_intact_over_observed": ratio("E_bases_intact", s2),
            "D6321_simple_over_observed": ratio("E_bases_simple", s2),
            "f_cos": {"D6311": rng("f_cos", s1), "D6321": rng("f_cos", s2)},
            "L_eff_implied_bp": rng("L_eff_implied"),
            "phi": rng("phi_cos_adaptable"), "beta": rng("beta_breaks_per_genome"),
            "escape_frac_reads": rng("escape_frac_reads"),
            "escape_frac_lambda_bases": rng("escape_frac_lambda_bases"),
            "lambda_mean_len_bp": rng("lambda_mean_len"),
            "eject_mode_bp": rng("eject_mode_bp"),
            "unblocked_median_bp": rng("unblocked_median_bp"),
            "first_6h_vs_last_12h": {k: [rng(f"{k}_first6h"), rng(f"{k}_last12h")] for k in (
                "escape_frac", "kept_frac", "unblocked_median_bp", "sample_per_lambda_derived")},
            "D6311_per_organism": slopes,
            "per_library": [{k: r[k] for k in (
                "library", "include_in_headline", "E_bases_obs", "E_bases_simple", "E_bases_intact",
                "lambda_mean_len", "escape_frac_reads", "escape_frac_lambda_bases", "f_cos",
                "L_eff_implied", "phi_cos_adaptable", "beta_breaks_per_genome", "beta_upper_bound",
                "unblocked_median_bp")} for r in rows],
        },
    }
    (figdir / f"{FIG_ID}.json").write_text(json.dumps(payload, indent=2, default=float))


def table(rows, figdir):
    """Supplementary Table S15, unrounded; the manuscript rounds."""
    cols = [("library", "library"), ("include_in_headline", "include_in_headline"),
            ("E_bases_obs", "E_bases_observed"), ("E_bases_simple", "E_bases_465bp"),
            ("E_bases_intact", "E_bases_intact"), ("lambda_mean_len", "mean_lambda_derived_read_bp"),
            ("escape_frac_reads", "escaped_reads_frac"),
            ("escape_frac_lambda_bases", "escaped_lambda_bases_frac"),
            ("f_cos", "f_cos"), ("L_eff_implied", "L_eff_bp"), ("phi_cos_adaptable", "phi"),
            ("beta_breaks_per_genome", "beta"), ("beta_upper_bound", "beta_bound")]
    trows = [{new: r[old] for old, new in cols} for r in rows]
    write_csv(figdir / f"{TABLE_ID}.csv", trows)
    payload = {
        "id": TABLE_ID, "display_type": "calculation",
        "title": "Carrier-depletion model, per library",
        "caption": (
            "E is base-weighted enrichment: observed; predicted with 465 bp per carrier "
            "capture (the pooled median carrier read); and predicted for intact lambda with "
            "the measured mean bases per lambda-derived read. Escapes are lambda-derived "
            "reads of 2 kb or longer, as a share of lambda-derived reads and of their bases. "
            "f_cos, L_eff, phi and beta are defined in Supplementary Section S10; the bound on "
            "beta is -ln F. D6321 r0 is shown for completeness and excluded from all headline "
            "statistics. Every column counts lambda-derived reads: those assigned to lambda "
            "plus those tied between lambda and another organism."),
        "source_files": ["decomposition.tsv"],
        "software": software(),
        "metrics": {"n_rows": len(trows)},
    }
    (figdir / f"{TABLE_ID}.json").write_text(json.dumps(payload, indent=2, default=float))


def write_csv(path, rows):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})


def software():
    return [f"python {sys.version.split()[0]}", f"numpy {np.__version__}",
            f"matplotlib {matplotlib.__version__}", "minimap2 (tools image; see minimap2.version)"]


if __name__ == "__main__":
    main(*sys.argv[1:])
