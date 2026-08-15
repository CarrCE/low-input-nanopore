#!/usr/bin/env python3
"""
Enforce the display-item contract (docs/display-items.md).

For every figure produced by this repository, check that:

  1. a vector PDF exists;
  2. a JSON sidecar exists and carries id / title / caption / source_files /
     software;
  3. a CSV exists, and its row count equals the number of points the figure
     actually draws, as declared by `metrics.n_plotted_points` in the JSON.

Check 3 is the one that matters. A CSV of summary statistics next to a figure
looks correct and is not: it does not let a reader redraw the figure. That
mistake was made twice in this project -- read-length quantiles standing in for
a histogram, and a per-organism statistics table standing in for depth profiles
-- so the invariant is now machine-checked rather than trusted.

Additional statistics tables are welcome and are ignored here; they live at
`<id>_summary.csv` so the canonical `<id>.csv` keeps a single meaning.

Exit status is non-zero if any item fails, so this can gate a release.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REQUIRED_JSON_KEYS = ("id", "title", "caption", "source_files", "software")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # "comparison/figures" was the pre-move path. It has not existed since the
    # comparison module was split into bin/ and assets/, so this scan printed
    # "[skip] comparison/figures does not exist" and still reported conformance
    # -- the per-femtogram comparison figure was never checked by the gate
    # written to check it.
    p.add_argument("--dirs", nargs="+", default=["results/summary", "results/comparison"],
                   help="directories to scan for display items")
    # A directory named explicitly is one the caller expects to exist; skipping
    # it silently is how the above went unnoticed.
    p.add_argument("--allow-missing-dirs", action="store_true",
                   help="treat a named directory that does not exist as empty "
                        "rather than as an error")
    p.add_argument("--no-require-pdf", dest="require_pdf", action="store_false",
                   help="do not require a vector PDF for plotted items")
    p.set_defaults(require_pdf=True)
    return p.parse_args()


def count_csv_rows(path: Path) -> int:
    with open(path, newline="") as fh:
        return max(0, sum(1 for _ in csv.reader(fh)) - 1)


def main():
    args = parse_args()
    failures, checked = [], 0

    for d in args.dirs:
        base = Path(d)
        if not base.is_dir():
            if args.allow_missing_dirs:
                print(f"[skip] {d} does not exist")
                continue
            failures.append((d, ["directory does not exist; nothing in it was "
                                 "checked. Generate it, or pass "
                                 "--allow-missing-dirs."]))
            continue
        for js in sorted(base.glob("*.json")):
            item = js.stem
            if item.endswith("_summary"):
                continue
            checked += 1
            problems = []
            meta = json.loads(js.read_text())

            for k in REQUIRED_JSON_KEYS:
                if not meta.get(k):
                    problems.append(f"JSON missing/empty '{k}'")

            csv_path = base / f"{item}.csv"
            pdf_path = base / f"{item}.pdf"

            if not csv_path.exists():
                problems.append("no canonical CSV")
            if args.require_pdf and not pdf_path.exists():
                # estimated_control is a calculation, not a figure; it declares
                # no plotted points and is exempt.
                if (meta.get("metrics") or {}).get("n_plotted_points") is not None:
                    problems.append("no vector PDF")

            if meta.get("display_type") == "calculation":
                # `problems` accumulated above (missing title/caption/software)
                # used to be discarded here, so a calculation item passed with
                # any metadata defect as long as its CSV existed.
                if problems:
                    failures.append((item, problems))
                    print(f"[FAIL] {item:<24} {'':>7}  " + "; ".join(problems))
                else:
                    print(f"[ok]   {item:<24} {'':>7}  calculation (no figure); JSON+CSV present")
                continue

            declared = (meta.get("metrics") or {}).get("n_plotted_points")
            if csv_path.exists():
                rows = count_csv_rows(csv_path)
                if declared is None:
                    problems.append(
                        "JSON does not declare metrics.n_plotted_points, so the CSV "
                        "cannot be checked against the figure")
                elif rows != declared:
                    problems.append(
                        f"CSV has {rows} rows but the figure draws {declared} points "
                        "-- the CSV is not the plotted data")
                else:
                    print(f"[ok]   {item:<24} {rows:>7,} plotted points, PDF+CSV+JSON present")

            if problems:
                failures.append((item, problems))

    print()
    for item, problems in failures:
        print(f"[FAIL] {item}")
        for p in problems:
            print(f"         - {p}")

    print(f"\n{checked - len(failures)}/{checked} display items conform "
          f"to docs/display-items.md")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
