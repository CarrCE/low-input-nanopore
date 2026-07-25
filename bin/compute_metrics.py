#!/usr/bin/env python3
"""
Turn per-organism read/base counts into the enrichment metrics reported in the
manuscript.

Metric definitions (stated explicitly because the comparison to prior work
hinges on them)
---------------------------------------------------------------------------
input_sample_fraction
    library_dna_ng / (library_dna_ng + carrier_dna_ng). The fraction of the
    mass entering library prep that is sample rather than carrier.

output_sample_fraction
    sample bases (or reads) / all basecalled bases (or reads), where the
    denominator is every read in the FASTQ including carrier, contaminant,
    ambiguous and unmapped. Counts come from the assignment step, which
    accounts for every read exactly once, so this denominator is exact.

enrichment
    output_sample_fraction / input_sample_fraction. This is the quantity
    reported as ">100x enrichment": how much depletion-mode adaptive sampling
    raised the sample's share of the output above its share of the input.

reads_per_fg / bases_per_fg
    reads (or bases) assigned to an organism, divided by the femtograms of that
    organism's DNA that entered library prep, i.e.
        library_dna_ng * 1e6 * theoretical_dna_fraction
    For the "All organisms" row the denominator is library_dna_ng * 1e6.
    This matches the convention used for the prior-study comparison, where the
    denominator is DNA into library prep ("post-extraction").

bases
    Full read length of every read assigned to the organism, NOT the aligned
    span. This is deliberate: the prior studies being compared against count
    Kraken2-classified read lengths, so counting aligned bases here would
    understate this study relative to them. Aligned bases are reported
    alongside as `aligned_bases` for internal use (coverage, identity).
"""
from __future__ import annotations

import argparse
import json
import math
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True, help="<prefix>.counts.tsv from assign_reads.py")
    p.add_argument("--reference-tsv", required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--experiment", required=True)
    p.add_argument("--replicate", required=True)
    p.add_argument("--library-dna-ng", default="",
                   help="sample DNA into library prep; blank if unquantified")
    p.add_argument("--carrier-dna-ng", default="0")
    p.add_argument("--mode", default="competitive")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-json", required=True)
    return p.parse_args()


def read_tsv(path):
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def fnum(x, default=None):
    """Parse a possibly-blank numeric field."""
    if x is None or str(x).strip() == "":
        return default
    try:
        return float(x)
    except ValueError:
        return default


def main():
    args = parse_args()

    counts = {r["organism"]: r for r in read_tsv(args.counts)}
    refs = read_tsv(args.reference_tsv)

    theoretical = {r["organism"]: fnum(r.get("theoretical_dna_fraction"), 0.0) for r in refs}
    roles = {r["organism"]: r["role"] for r in refs}

    library_ng = fnum(args.library_dna_ng)
    carrier_ng = fnum(args.carrier_dna_ng, 0.0) or 0.0

    # ---- totals over every read in the FASTQ -------------------------------
    total_reads = sum(int(r["reads"]) for r in counts.values())
    total_bases = sum(int(r["read_bases"]) for r in counts.values())

    def bucket(role):
        rs = bs = 0
        for org, r in counts.items():
            org_role = roles.get(org, r.get("role", "unknown"))
            if org.startswith("ambiguous:"):
                org_role = "ambiguous"
            elif org == "unassigned":
                org_role = "unassigned"
            if org_role == role:
                rs += int(r["reads"])
                bs += int(r["read_bases"])
        return rs, bs

    sample_reads, sample_bases = bucket("sample")
    carrier_reads, carrier_bases = bucket("carrier")
    contam_reads, contam_bases = bucket("contaminant")
    ambig_reads, ambig_bases = bucket("ambiguous")
    unassigned_reads, unassigned_bases = bucket("unassigned")

    input_sample_fraction = (library_ng / (library_ng + carrier_ng)
                             if library_ng is not None and (library_ng + carrier_ng) > 0
                             else None)
    out_frac_bases = sample_bases / total_bases if total_bases else None
    out_frac_reads = sample_reads / total_reads if total_reads else None

    enrichment_bases = (out_frac_bases / input_sample_fraction
                        if input_sample_fraction else None)
    enrichment_reads = (out_frac_reads / input_sample_fraction
                        if input_sample_fraction else None)

    # ---- per-organism rows -------------------------------------------------
    out_rows = []

    def emit(organism, role, reads, read_bases, aligned_bases, theo_frac):
        dna_fg = (library_ng * 1e6 * theo_frac
                  if library_ng is not None and theo_frac is not None else None)
        measured = read_bases / sample_bases if sample_bases and role == "sample" else None
        out_rows.append({
            "sample_id": args.sample_id,
            "experiment": args.experiment,
            "replicate": args.replicate,
            "mode": args.mode,
            "organism": organism,
            "role": role,
            "reads": reads,
            "bases": read_bases,
            "aligned_bases": aligned_bases,
            "theoretical_dna_fraction": theo_frac,
            "measured_base_fraction": measured,
            "dna_fg": dna_fg,
            "reads_per_fg": (reads / dna_fg) if dna_fg else None,
            "bases_per_fg": (read_bases / dna_fg) if dna_fg else None,
        })

    for ref in refs:
        org = ref["organism"]
        c = counts.get(org)
        reads = int(c["reads"]) if c else 0
        rb = int(c["read_bases"]) if c else 0
        ab = int(c["aligned_bases"]) if c else 0
        emit(org, ref["role"], reads, rb, ab, theoretical.get(org))

    # Ambiguous and unassigned classes are carried through explicitly so the
    # per-organism table always reconciles to the FASTQ totals.
    for org, r in sorted(counts.items()):
        if org.startswith("ambiguous:") or org == "unassigned":
            emit(org, "ambiguous" if org.startswith("ambiguous:") else "unassigned",
                 int(r["reads"]), int(r["read_bases"]), int(r["aligned_bases"]), None)

    # Aggregate "All organisms" row: the headline per-fg number.
    all_dna_fg = library_ng * 1e6 if library_ng is not None else None
    out_rows.append({
        "sample_id": args.sample_id, "experiment": args.experiment,
        "replicate": args.replicate, "mode": args.mode,
        "organism": "All organisms", "role": "sample_total",
        "reads": sample_reads, "bases": sample_bases,
        "aligned_bases": sum(int(counts[o]["aligned_bases"]) for o in counts
                             if roles.get(o) == "sample"),
        "theoretical_dna_fraction": 1.0,
        "measured_base_fraction": 1.0 if sample_bases else None,
        "dna_fg": all_dna_fg,
        "reads_per_fg": (sample_reads / all_dna_fg) if all_dna_fg else None,
        "bases_per_fg": (sample_bases / all_dna_fg) if all_dna_fg else None,
    })

    cols = ["sample_id", "experiment", "replicate", "mode", "organism", "role",
            "reads", "bases", "aligned_bases", "theoretical_dna_fraction",
            "measured_base_fraction", "dna_fg", "reads_per_fg", "bases_per_fg"]
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(
                "" if r.get(c) is None else
                (f"{r[c]:.10g}" if isinstance(r[c], float) else str(r[c]))
                for c in cols) + "\n")

    # ---- reconciliation check ---------------------------------------------
    # Every read in the FASTQ must land in exactly one bucket. If this fails the
    # per-fg numbers are not trustworthy, so fail loudly rather than publish.
    bucket_reads = (sample_reads + carrier_reads + contam_reads
                    + ambig_reads + unassigned_reads)
    if bucket_reads != total_reads:
        sys.exit(f"error: read accounting does not reconcile: buckets sum to "
                 f"{bucket_reads:,} but total is {total_reads:,}")

    summary = {
        "sample_id": args.sample_id,
        "experiment": args.experiment,
        "replicate": args.replicate,
        "mode": args.mode,
        "input": {
            "library_dna_ng": library_ng,
            "carrier_dna_ng": carrier_ng,
            "input_sample_fraction": input_sample_fraction,
        },
        "totals": {
            "reads": total_reads,
            "bases": total_bases,
        },
        "by_role": {
            "sample":      {"reads": sample_reads,     "bases": sample_bases},
            "carrier":     {"reads": carrier_reads,    "bases": carrier_bases},
            "contaminant": {"reads": contam_reads,     "bases": contam_bases},
            "ambiguous":   {"reads": ambig_reads,      "bases": ambig_bases},
            "unassigned":  {"reads": unassigned_reads, "bases": unassigned_bases},
        },
        "output_sample_fraction": {
            "reads": out_frac_reads,
            "bases": out_frac_bases,
        },
        "enrichment": {
            "reads": enrichment_reads,
            "bases": enrichment_bases,
        },
        "per_fg": {
            "reads_per_fg": (sample_reads / all_dna_fg) if all_dna_fg else None,
            "bases_per_fg": (sample_bases / all_dna_fg) if all_dna_fg else None,
        },
    }
    with open(args.out_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    def pct(x):
        return "n/a" if x is None else f"{100 * x:.4f}%"

    print(f"[metrics] {args.sample_id} ({args.mode})")
    print(f"  total            {total_reads:>12,} reads  {total_bases:>15,} bases")
    print(f"  sample           {sample_reads:>12,} reads  {sample_bases:>15,} bases")
    print(f"  carrier          {carrier_reads:>12,} reads  {carrier_bases:>15,} bases")
    print(f"  contaminant      {contam_reads:>12,} reads  {contam_bases:>15,} bases")
    print(f"  ambiguous        {ambig_reads:>12,} reads  {ambig_bases:>15,} bases")
    print(f"  unassigned       {unassigned_reads:>12,} reads  {unassigned_bases:>15,} bases")
    print(f"  input sample fraction  {pct(input_sample_fraction)}")
    print(f"  output sample fraction {pct(out_frac_bases)} (bases)")
    if enrichment_bases:
        print(f"  ENRICHMENT             {enrichment_bases:.1f}x (bases), "
              f"{enrichment_reads:.1f}x (reads)")


if __name__ == "__main__":
    main()
