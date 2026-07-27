#!/usr/bin/env python3
"""
Check assets/measurements.tsv against the samplesheets, and report what is
still PENDING.

The measurements file carries the experimental quantities every per-femtogram
and enrichment number divides by. This asserts the properties that make it
trustworthy:

  * every sample in a samplesheet has a measurements row, and vice versa
  * no headline sample is missing a mass or a basis
  * every basis is one of the recognised values, so a typo cannot silently
    become an unrecognised-but-accepted provenance
  * a sample whose mass basis is not a real measurement is not marked
    include_in_headline
  * every exclusion carries a stated reason

That last pair is the point of the file: a nominal or extrapolated denominator
feeding a headline statistic should be impossible to introduce quietly.

Exit codes: 0 clean, 1 a check failed, 2 unresolved PENDINGs remain (which is
expected while the file is still being filled in, and is reported separately
from a genuine inconsistency).

Usage:
    python3 bin/check_measurements.py
    python3 bin/check_measurements.py --strict   # PENDING is a failure too
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

MEASURED = {"qubit_hs", "qubit_hs_measured", "qubit_hs_inferred",
            "qubit_hs_shared_tube", "raw_fluorescence_extrapolated"}
NOT_MEASURED = {"nominal_unmeasured", "PENDING"}
SAMPLE_BASES = MEASURED | NOT_MEASURED
CARRIER_BASES = {"direct_addition", "dilution_series", "PENDING"}

# A basis that is an extrapolation below an assay's validated range is a
# measurement, but a weaker one. Feeding it to a headline statistic is a
# decision, so it is surfaced rather than silently allowed.
WEAK = {"raw_fluorescence_extrapolated", "qubit_hs_shared_tube"}


def read_tsv(path):
    with open(path) as fh:
        rows = [l for l in fh if not l.startswith("#") and l.strip()]
    return list(csv.DictReader(rows, delimiter="\t"))


def read_samplesheets(d):
    ids = {}
    for p in sorted(Path(d).glob("*.csv")):
        if p.name in ("test.csv", "raghavendra_2023.csv", "all.csv"):
            continue
        with open(p) as fh:
            rows = [l for l in fh if not l.startswith("#") and l.strip()]
        for r in csv.DictReader(rows):
            if r.get("sample_id"):
                ids[r["sample_id"]] = p.name
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--measurements", default=None)
    ap.add_argument("--samplesheets", default=None)
    ap.add_argument("--strict", action="store_true",
                    help="treat unresolved PENDING values as failures")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    meas_path = Path(args.measurements) if args.measurements else root / "assets/measurements.tsv"
    sheet_dir = Path(args.samplesheets) if args.samplesheets else root / "assets/samplesheets"

    rows = read_tsv(meas_path)
    sheets = read_samplesheets(sheet_dir)
    errors, pendings, warnings = [], [], []

    ids = {r["sample_id"] for r in rows}
    for sid, src in sheets.items():
        if sid not in ids:
            errors.append(f"{sid} is in {src} but has no measurements row")
    for sid in ids - set(sheets):
        errors.append(f"{sid} has a measurements row but is in no samplesheet")

    for r in rows:
        sid = r["sample_id"]
        head = r["include_in_headline"].strip() == "1"
        sb, cb = r["sample_dna_basis"].strip(), r["carrier_dna_basis"].strip()

        if sb not in SAMPLE_BASES:
            errors.append(f"{sid}: unrecognised sample_dna_basis {sb!r}")
        if cb not in CARRIER_BASES:
            errors.append(f"{sid}: unrecognised carrier_dna_basis {cb!r}")

        for field in ("sample_dna_ng", "carrier_dna_ng"):
            v = r[field].strip()
            if v == "PENDING" or v == "":
                (pendings if not head else errors if v == "" else pendings).append(
                    f"{sid}: {field} is {'blank' if not v else 'PENDING'}"
                    + (" and this sample is in the headline" if head else ""))
            else:
                try:
                    if float(v) <= 0:
                        errors.append(f"{sid}: {field} is not positive ({v})")
                except ValueError:
                    errors.append(f"{sid}: {field} is not a number ({v!r})")

        for field, val in (("sample_dna_basis", sb), ("carrier_dna_basis", cb)):
            if val == "PENDING":
                pendings.append(f"{sid}: {field} is PENDING"
                                + (" and this sample is in the headline" if head else ""))

        # The check this file exists for.
        if head and sb == "nominal_unmeasured":
            errors.append(f"{sid}: contributes to the headline but its sample mass "
                          f"is nominal, not measured")
        if head and sb in WEAK:
            warnings.append(f"{sid}: contributes to the headline on a mass that is "
                            f"a weaker measurement than a direct in-range reading ({sb})")
        if not head and not r.get("include_reason", "").strip():
            errors.append(f"{sid}: excluded from the headline with no stated reason")

    for label, items in (("ERROR", errors), ("WARNING", warnings), ("PENDING", pendings)):
        for m in items:
            print(f"{label}: {m}")

    n_head = sum(1 for r in rows if r["include_in_headline"].strip() == "1")
    print(f"\n{len(rows)} rows, {n_head} in the headline; "
          f"{len(errors)} errors, {len(warnings)} warnings, {len(pendings)} pending")

    if errors:
        return 1
    if pendings:
        return 1 if args.strict else 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
