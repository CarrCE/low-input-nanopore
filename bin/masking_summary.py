#!/usr/bin/env python3
"""Collect the per-sample human-masking statistics into one table.

`--mask_human` writes a `<sample>.human_stats.json` per replicate. This gathers
them into `results/summary/human_masking.tsv`, which is the Supplementary table
of what masking did, plus a JSON sidecar carrying the rule, the provenance and
a publication-ready caption.

Deliberately reports BOTH read and base percentages. Masking replaces bases
with N in place, so no read is ever removed: "0.33% of reads masked" and
"0.02% of bases masked" are both true and mean different things, and quoting
only the first would overstate what left the dataset.

Usage:
  masking_summary.py --stats results/*/human/*.human_stats.json \
                     --out results/summary/human_masking.tsv
"""

import argparse
import glob
import json
import os
import sys

# Order matters: this is the column order of the published table.
COLUMNS = [
    "sample_id",
    "reads",
    "bases",
    "hrrt_flagged",
    "hrrt_flagged_pct",
    "rescued_reads",
    "masked_reads",
    "masked_reads_pct",
    "masked_bases",
    "masked_bases_pct",
    "fully_masked_reads",
    "partially_masked_reads",
    "chimeric_masks",
    "rescued_median_bp",
    "masked_median_bp",
]


def row_from(stats):
    lens = stats.get("flagged_read_length", {})
    return {
        "sample_id": stats["sample_id"],
        "reads": stats["reads"],
        "bases": stats["bases"],
        "hrrt_flagged": stats["hrrt_flagged"],
        "hrrt_flagged_pct": stats["hrrt_flagged_pct"],
        "rescued_reads": stats["rescued_reads"],
        "masked_reads": stats["masked_reads"],
        "masked_reads_pct": stats["masked_reads_pct"],
        "masked_bases": stats["masked_bases"],
        "masked_bases_pct": stats["masked_bases_pct"],
        "fully_masked_reads": stats["fully_masked_reads"],
        "partially_masked_reads": stats["partially_masked_reads"],
        "chimeric_masks": stats["chimeric_masks"],
        "rescued_median_bp": lens.get("rescued_median", 0),
        "masked_median_bp": lens.get("masked_median", 0),
    }


def totals(rows):
    """Pooled row. Percentages are recomputed from the pooled counts, never
    averaged: the replicates differ in size by more than a factor of two, so a
    mean of per-sample percentages is not the percentage of the pool."""
    t = {c: 0 for c in COLUMNS}
    t["sample_id"] = "all"
    for c in ("reads", "bases", "hrrt_flagged", "rescued_reads", "masked_reads",
              "masked_bases", "fully_masked_reads", "partially_masked_reads",
              "chimeric_masks"):
        t[c] = sum(r[c] for r in rows)
    t["hrrt_flagged_pct"] = round(100.0 * t["hrrt_flagged"] / t["reads"], 4) if t["reads"] else 0.0
    t["masked_reads_pct"] = round(100.0 * t["masked_reads"] / t["reads"], 4) if t["reads"] else 0.0
    t["masked_bases_pct"] = round(100.0 * t["masked_bases"] / t["bases"], 4) if t["bases"] else 0.0
    # A pooled median cannot be derived from per-sample medians; leave blank
    # rather than print a mean of medians and let it be read as one.
    t["rescued_median_bp"] = ""
    t["masked_median_bp"] = ""
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", nargs="+", required=True,
                    help="human_stats.json files, or globs")
    ap.add_argument("--out", required=True, help="output TSV")
    args = ap.parse_args()

    paths = []
    for pat in args.stats:
        hits = sorted(glob.glob(pat)) if any(c in pat for c in "*?[") else [pat]
        if not hits:
            sys.exit(f"error: no files matched {pat}")
        paths.extend(hits)

    loaded = []
    for p in paths:
        with open(p) as fh:
            loaded.append(json.load(fh))

    rows = sorted((row_from(s) for s in loaded), key=lambda r: r["sample_id"])
    if not rows:
        sys.exit("error: no statistics files to summarise")
    rows.append(totals(rows))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")

    # The rule and provenance are identical across samples by construction;
    # assert that rather than quietly reporting the first one.
    rules = {json.dumps(s["rule"], sort_keys=True) for s in loaded}
    provs = {json.dumps(s["provenance"], sort_keys=True) for s in loaded}
    if len(rules) != 1:
        sys.exit("error: samples were masked under different rules; refusing "
                 "to summarise them in one table")
    prov = json.loads(sorted(provs)[0])
    if len(provs) != 1:
        # Only the assignments filename varies legitimately (it is per sample).
        varying = {k for p in provs for k, v in json.loads(p).items()
                   if v != prov[k]}
        if varying - {"assignments"}:
            sys.exit(f"error: samples disagree on provenance fields: {sorted(varying)}")

    total = rows[-1]
    meta = {
        "id": "human_masking",
        "title": "Human read masking, per replicate",
        "caption": (
            "Human sequence masking applied before public deposition. Reads "
            "were screened with NCBI's Human Read Removal Tool (HRRT) and a "
            "flagged read was released intact only where competitive "
            "assignment attributed it to a community organism; everything "
            "else was masked with N over the interval no organism accounts "
            f"for. Masking preserves read count and read length exactly, so "
            f"the {total['masked_reads']:,} masked reads "
            f"({total['masked_reads_pct']}% of {total['reads']:,}) remain in "
            f"the deposited files with {total['masked_bases']:,} bases "
            f"({total['masked_bases_pct']}%) replaced. HRRT flagged "
            f"{total['hrrt_flagged']:,} reads ({total['hrrt_flagged_pct']}%); "
            f"{total['rescued_reads']:,} of those were rescued by positive "
            "attribution. Reads confidently assignable to human were masked; "
            "this does not establish that no human sequence remains."
        ),
        "rule": json.loads(sorted(rules)[0]),
        "provenance": prov,
        "source": [os.path.relpath(p) for p in paths],
        "columns": COLUMNS,
    }
    side = os.path.splitext(args.out)[0] + ".json"
    with open(side, "w") as fh:
        json.dump(meta, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(f"[masking-summary] {len(rows) - 1} samples -> {args.out}")
    print(f"[masking-summary] pooled: {total['hrrt_flagged']:,} flagged "
          f"({total['hrrt_flagged_pct']}%), {total['rescued_reads']:,} rescued, "
          f"{total['masked_reads']:,} masked ({total['masked_reads_pct']}%), "
          f"{total['masked_bases']:,} bases masked ({total['masked_bases_pct']}%)")


if __name__ == "__main__":
    main()
