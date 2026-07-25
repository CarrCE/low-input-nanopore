#!/usr/bin/env python3
"""
Estimate what a no-adaptive-sampling control would have yielded, from Mojarro
et al. 2019.

No matched control library was sequenced for this study with adaptive sampling
disabled. Mojarro et al. 2019 is the closest available substitute, and it is a
close one: the same genomic carrier (1 ug of lambda gDNA), the same picogram
sample regime, a standard 48 h MinION run with NO real-time selective
sequencing, and a published whole-run read and base total. It is therefore a
carrier-sequencing experiment differing from this study principally in the
absence of adaptive sampling.

The comparison uses exactly the enrichment definition applied to this study's
own data (see bin/compute_metrics.py):

    enrichment = (sample share of output) / (sample share of input mass)

What it shows is a sign change, not merely a smaller number. Without adaptive
sampling the sample is *depleted* on a base basis, because carrier molecules
sequence to full length while the scarce sample molecules do not dominate the
base budget. With depletion-mode adaptive sampling the carrier is truncated on
rejection and the relationship inverts. The read-length statistics on both sides
are what make that mechanism legible rather than assumed.

This is an estimate across experiments, not a control. Every difference is
recorded in the JSON sidecar.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

# --- Mojarro et al. 2019, Astrobiology 19(9):1139-1152, doi:10.1089/ast.2018.1929
# Table 1 "Low-Input Carrier Sequencing Metrics" + Methods.
MOJARRO = {
    "sample_pg": 2.0,               # abstract: "2 pg of purified B. subtilis spore DNA"
    "carrier_pg": 1_000_000.0,      # Methods: 1 ug lambda gDNA; protocol replicated unmodified
    "total_reads": 1_303_007,       # Table 1, "All reads" = whole 48 h run
    "total_bases": 8_698_026_598,
    "target_reads": 5,              # Table 1, "B. subtilis reads"
    "target_bases": 5_270,
}


def poisson_ci(k: int, alpha: float = 0.05):
    """Garwood exact Poisson interval; k is small here so this matters."""
    from scipy.stats import chi2
    lo = chi2.ppf(alpha / 2, 2 * k) / 2 if k > 0 else 0.0
    hi = chi2.ppf(1 - alpha / 2, 2 * k + 2) / 2
    return float(lo), float(hi)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--this-study", nargs="*", default=[],
                   help="name=enrichment_bases,enrichment_reads (e.g. lowinput_s2=43.01,18.39)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--basename", default="estimated_control")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    m = MOJARRO
    f_in = m["sample_pg"] / (m["sample_pg"] + m["carrier_pg"])
    f_out_bases = m["target_bases"] / m["total_bases"]
    f_out_reads = m["target_reads"] / m["total_reads"]
    enr_bases = f_out_bases / f_in
    enr_reads = f_out_reads / f_in

    lo, hi = poisson_ci(m["target_reads"])
    enr_bases_lo = enr_bases * lo / m["target_reads"]
    enr_bases_hi = enr_bases * hi / m["target_reads"]

    mean_run_read = m["total_bases"] / m["total_reads"]
    mean_target_read = m["target_bases"] / m["target_reads"]

    print("Estimated no-adaptive-sampling control (Mojarro et al. 2019)")
    print(f"  input sample fraction        {f_in:.4e}")
    print(f"  output fraction, bases       {f_out_bases:.4e}  -> enrichment {enr_bases:.3f}x")
    print(f"  output fraction, reads       {f_out_reads:.4e}  -> enrichment {enr_reads:.3f}x")
    print(f"  95% Poisson CI on 5 reads    [{lo:.2f}, {hi:.2f}]"
          f"  -> bases enrichment [{enr_bases_lo:.3f}, {enr_bases_hi:.3f}]x")
    print(f"  whole-run mean read length   {mean_run_read:,.0f} bp (carrier-dominated)")
    print(f"  target mean read length      {mean_target_read:,.0f} bp")
    print("  => without adaptive sampling the sample is DEPLETED on a base basis")

    comparisons = []
    for spec in args.this_study:
        name, vals = spec.split("=")
        eb, er = (float(x) for x in vals.split(","))
        comparisons.append({
            "experiment": name,
            "enrichment_bases": eb,
            "enrichment_reads": er,
            "fold_over_estimated_control_bases": eb / enr_bases,
            "fold_over_estimated_control_reads": er / enr_reads,
            "fold_over_estimated_control_bases_ci": [eb / enr_bases_hi, eb / enr_bases_lo],
        })
        print(f"  {name}: {eb:.2f}x bases -> {eb/enr_bases:.0f}x the estimated control "
              f"(95% CI {eb/enr_bases_hi:.0f}-{eb/enr_bases_lo:.0f}x); "
              f"{er:.2f}x reads -> {er/enr_reads:.1f}x")

    payload = {
        "id": args.basename,
        "title": "Estimated no-adaptive-sampling control from Mojarro et al. 2019",
        "caption": (
            "No control library was sequenced for this study with adaptive sampling "
            "disabled. Mojarro et al. 2019 provides the closest published substitute: "
            "the same 1 ug lambda genomic carrier, a 2 pg DNA sample, and a standard "
            "48 h MinION run with no real-time selective sequencing. Applying this "
            "study's enrichment definition to their Table 1 gives 0.30x on a base "
            "basis - that is, without adaptive sampling the sample's share of output "
            "bases is about three-fold LOWER than its share of input mass, because "
            "carrier molecules sequence to full length (whole-run mean 6,675 bp) while "
            "target molecules do not (mean 1,054 bp). Depletion-mode adaptive sampling "
            "inverts that relationship, truncating carrier reads on rejection. The "
            "comparison is an estimate across experiments, not a control."),
        "source": ("Mojarro A, Hachey J, Bailey R, Brown M, Doebler R, Ruvkun G, "
                   "Zuber MT, Carr CE. Astrobiology 2019;19(9):1139-1152. "
                   "doi:10.1089/ast.2018.1929. Table 1 and Methods."),
        "inputs": m,
        "estimated_control": {
            "input_sample_fraction": f_in,
            "output_fraction_bases": f_out_bases,
            "output_fraction_reads": f_out_reads,
            "enrichment_bases": enr_bases,
            "enrichment_reads": enr_reads,
            "enrichment_bases_95ci": [enr_bases_lo, enr_bases_hi],
            "whole_run_mean_read_bp": mean_run_read,
            "target_mean_read_bp": mean_target_read,
        },
        "comparisons": comparisons,
        "caveats": [
            "Different flow cell chemistry: R9.4 (Mojarro) vs R10.4.1 (this study).",
            "Different basecaller: Albacore 1.10 vs dorado sup v5.0.0.",
            "Different sample: B. subtilis spore DNA vs ZymoBIOMICS mock communities.",
            "Different target-detection method: CarrierSeq bioinformatic filtering vs "
            "competitive alignment against known reference genomes.",
            "The estimate rests on 5 target reads, so its Poisson interval is wide; the "
            "interval is reported and never collapsed to a point.",
            "Carrier mass matches (1 ug lambda) and both are picogram-to-nanogram "
            "carrier-sequencing experiments, which is what makes the comparison "
            "meaningful at all.",
            "This is not a substitute for a matched control run on the same platform.",
        ],
    }
    (outdir / f"{args.basename}.json").write_text(json.dumps(payload, indent=2))

    with open(outdir / f"{args.basename}.csv", "w") as fh:
        fh.write("quantity,value\n")
        fh.write(f"input_sample_fraction,{f_in!r}\n")
        fh.write(f"output_fraction_bases,{f_out_bases!r}\n")
        fh.write(f"output_fraction_reads,{f_out_reads!r}\n")
        fh.write(f"enrichment_bases,{enr_bases!r}\n")
        fh.write(f"enrichment_reads,{enr_reads!r}\n")
        fh.write(f"enrichment_bases_ci_lo,{enr_bases_lo!r}\n")
        fh.write(f"enrichment_bases_ci_hi,{enr_bases_hi!r}\n")
        for c in comparisons:
            fh.write(f"{c['experiment']}_fold_over_control_bases,"
                     f"{c['fold_over_estimated_control_bases']!r}\n")
    print(f"\n[control] wrote {args.basename}.json/.csv in {outdir}")


if __name__ == "__main__":
    main()
