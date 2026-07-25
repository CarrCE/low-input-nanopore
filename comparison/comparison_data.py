#!/usr/bin/env python3
"""
Loader for the low-input comparison module.

Reads the two committed TSVs (`prior_studies.tsv`, `this_study.tsv`) and, when
asked, supersedes the seeded this-study rows with live pipeline output from
`results/<sample_id>/<mode>/<sample_id>.metrics.tsv`.

Nothing here opens the legacy .xlsx. The workbook is extracted exactly once by
`extract_workbook.py`; the TSVs are the versioned source of truth.

Conventions
-----------
    reads_per_fg = reads / (dna_pg * 1000)
    bases_per_fg = bases / (dna_pg * 1000)
    dna_pg is DNA into LIBRARY PREP (post-extraction) for every study.

Classifier selection
--------------------
Prior-study rows are duplicated per read/base-assignment method. Exactly one
`classifier` per study is plotted; see `DEFAULT_CLASSIFIERS`. This exists
because the legacy workbook mixed classifiers across the two axes for Zorzano
et al. 2025 (see README.md, defect (a)).

Run this file directly for a data-integrity check:
    python3 comparison_data.py --help
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

PRIOR_TSV = HERE / "prior_studies.tsv"
THIS_STUDY_TSV = HERE / "this_study.tsv"
DEFAULT_RESULTS_DIR = REPO / "results"

NUMERIC_COLUMNS = ["replicate_idx", "reads", "bases", "dna_pg",
                   "reads_per_fg", "bases_per_fg"]

STUDIES = ["Mojarro et al. 2019", "B. Raghavendra et al. 2023", "Zorzano et al. 2025"]

#: Which classifier variant is plotted for each prior study by default.
DEFAULT_CLASSIFIERS = {
    "Mojarro et al. 2019": "published_table1",
    "Zorzano et al. 2025": "kraken2_q1",          # internally consistent on both axes
    "B. Raghavendra et al. 2023": "kraken2_q10",  # only variant with base counts
}

ZORZANO_CHOICES = ["kraken2_q1", "published_squeezemeta", "legacy_hybrid_workbook"]
RAGHAVENDRA_CHOICES = ["kraken2_q10", "published"]

#: experiment id -> figure series label.
ROUND_BY_EXPERIMENT = {"lowinput_s1": "Round 1", "lowinput_s2": "Round 2"}

CONDITION_BY_EXPERIMENT = {
    "lowinput_s1": "D6311 community DNA + lambda carrier, depletion-mode adaptive sampling",
    "lowinput_s2": "D6321 whole cells, post-extraction, depletion-mode adaptive sampling",
}

HEADLINE_ORGANISM = "All organisms"


# ---------------------------------------------------------------------------
def to_float(s) -> float:
    """Correctly-rounded string -> float64.

    Python's float() is correctly rounded; pandas' `to_numeric` fast parser is
    not, and silently drops the last significant digit of values such as
    0.0004218604439799141. These tables are the audit trail for a figure, so
    every value must round-trip exactly.
    """
    if s is None:
        return float("nan")
    if isinstance(s, (int, float)) and not isinstance(s, bool):
        return float(s)
    t = str(s).strip()
    if t == "" or t.lower() in {"nan", "na", "n/a", "n.d."}:
        return float("nan")
    try:
        return float(t)
    except ValueError:
        return float("nan")


def _coerce(df: pd.DataFrame) -> pd.DataFrame:
    for c in NUMERIC_COLUMNS:
        if c in df.columns:
            df[c] = df[c].map(to_float).astype("float64")
    if "verified" in df.columns:
        df["verified"] = (df["verified"].astype(str).str.upper()
                          .map({"TRUE": True, "FALSE": False}).fillna(True))
    return df


def read_table(path: Path) -> pd.DataFrame:
    if not Path(path).exists():
        sys.exit(f"[error] missing input table: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    return _coerce(df)


def check_per_fg(df: pd.DataFrame, label: str, rtol: float = 1e-9) -> list[str]:
    """Verify the stored per-fg columns match reads/(dna_pg*1000)."""
    problems = []
    for _, r in df.iterrows():
        for count_col, rate_col in (("reads", "reads_per_fg"), ("bases", "bases_per_fg")):
            n, pg, rate = r.get(count_col), r.get("dna_pg"), r.get(rate_col)
            if pd.isna(n) or pd.isna(pg) or pd.isna(rate):
                continue
            expect = n / (pg * 1000.0)
            if abs(expect - rate) > rtol * max(1.0, abs(expect)):
                problems.append(
                    f"{label}: {r['study']} / {r['condition']} / {r['classifier']}: "
                    f"{rate_col}={rate!r} but {count_col}/(dna_pg*1000)={expect!r}")
    return problems


# ---------------------------------------------------------------------------
def load_prior(path: Path = PRIOR_TSV,
               zorzano_classifier: str = "kraken2_q1",
               raghavendra_classifier: str = "kraken2_q10"):
    """Return (plottable_df, notes). Notes are human-readable strings to print."""
    df = read_table(path)
    notes = []

    chosen = dict(DEFAULT_CLASSIFIERS)
    chosen["Zorzano et al. 2025"] = zorzano_classifier
    chosen["B. Raghavendra et al. 2023"] = raghavendra_classifier

    samples = df[df["row_type"] == "sample"].copy()
    n_aggregate = int((df["row_type"] == "aggregate").sum())
    if n_aggregate:
        notes.append(f"{n_aggregate} aggregate row(s) held in prior_studies.tsv are "
                     "never plotted (whole-study summaries, not samples).")

    keep = pd.Series(False, index=samples.index)
    for study, cls in chosen.items():
        available = set(samples.loc[samples["study"] == study, "classifier"])
        if cls not in available:
            sys.exit(f"[error] classifier {cls!r} not present for {study!r}. "
                     f"Available: {sorted(available)}")
        keep |= (samples["study"] == study) & (samples["classifier"] == cls)
    sel = samples[keep].copy()

    plottable = sel.dropna(subset=["reads_per_fg", "bases_per_fg"]).copy()
    dropped = sel[sel[["reads_per_fg", "bases_per_fg"]].isna().any(axis=1)]
    for study, grp in dropped.groupby("study"):
        notes.append(
            f"EXCLUDED {len(grp)} {study} row(s) with classifier "
            f"{chosen[study]!r}: no base counts are published, so bases/fg is "
            "undefined and the point cannot be placed on the y axis.")

    for study, cls in chosen.items():
        n = int((plottable["study"] == study).sum())
        notes.append(f"{study}: classifier={cls}, {n} point(s) plotted.")

    return plottable, notes


# ---------------------------------------------------------------------------
def _live_rows(results_dir: Path) -> pd.DataFrame:
    """Scan results/<sample_id>/<mode>/<sample_id>.metrics.tsv for headline rows."""
    rows = []
    for path in sorted(Path(results_dir).glob("*/*/*.metrics.tsv")):
        m = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
        if "organism" not in m.columns:
            continue
        head = m[m["organism"] == HEADLINE_ORGANISM]
        for _, r in head.iterrows():
            experiment = r.get("experiment", "")
            if experiment not in ROUND_BY_EXPERIMENT:
                continue          # smoke tests and other experiments are not figure data
            reads = to_float(r.get("reads"))
            bases = to_float(r.get("bases"))
            dna_fg = to_float(r.get("dna_fg"))
            rpf = to_float(r.get("reads_per_fg"))
            bpf = to_float(r.get("bases_per_fg"))
            if pd.isna(rpf) or pd.isna(bpf):
                continue          # unquantified input mass (e.g. lowinput_s2_r0)
            rows.append({
                "study": "This study",
                "study_short": "This study",
                "condition": CONDITION_BY_EXPERIMENT.get(experiment, experiment),
                "organism": HEADLINE_ORGANISM,
                "replicate_idx": to_float(r.get("replicate")),
                "replicate_label": f"replicate {r.get('replicate')}",
                "row_type": "sample",
                "classifier": f"minimap2_{r.get('mode', 'competitive')}",
                "reads": reads,
                "bases": bases,
                "dna_pg": dna_fg / 1000.0,
                "dna_basis": "qubit_library_input",
                "reads_per_fg": rpf,
                "bases_per_fg": bpf,
                "verified": True,
                "source": "pipeline_run",
                "source_detail": str(path),
                "provenance_note": (
                    f"Live pipeline output; '{HEADLINE_ORGANISM}' row of "
                    f"{path.name}. Supersedes the legacy_spreadsheet seed value "
                    f"for experiment {experiment}."),
                "round": ROUND_BY_EXPERIMENT[experiment],
                "experiment": experiment,
                "sample_id": r.get("sample_id", ""),
                "mode": r.get("mode", ""),
            })
    return pd.DataFrame(rows)


def load_this_study(path: Path = THIS_STUDY_TSV, results_dir: Path | None = None):
    """Seeded rows, with live pipeline rows superseding them per experiment."""
    df = read_table(path)
    notes = []
    if results_dir is None:
        notes.append(f"This study: {len(df)} seeded row(s) from {Path(path).name} "
                     "(source=legacy_spreadsheet). Pass --results-dir to use live "
                     "pipeline output instead.")
        return df, notes

    results_dir = Path(results_dir)
    if not results_dir.exists():
        notes.append(f"WARNING: --results-dir {results_dir} does not exist; "
                     "falling back to the seeded values.")
        return df, notes

    live = _live_rows(results_dir)
    if live.empty:
        notes.append(f"WARNING: no headline rows found under {results_dir}; "
                     "falling back to the seeded values.")
        return df, notes

    superseded = sorted(set(live["experiment"]))
    kept = df[~df["experiment"].isin(superseded)]
    merged = pd.concat([kept, live], ignore_index=True)
    notes.append(f"This study: {len(live)} live row(s) from {results_dir} superseded "
                 f"the seeded values for experiment(s) {', '.join(superseded)}; "
                 f"{len(kept)} seeded row(s) retained.")
    return _coerce(merged), notes


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Integrity check for the comparison module's data tables.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--prior-tsv", type=Path, default=PRIOR_TSV)
    ap.add_argument("--this-study-tsv", type=Path, default=THIS_STUDY_TSV)
    ap.add_argument("--results-dir", type=Path, default=None,
                    help="if given, load live pipeline metrics from this directory")
    args = ap.parse_args()

    prior = read_table(args.prior_tsv)
    mine = read_table(args.this_study_tsv)
    problems = check_per_fg(prior, "prior_studies.tsv") + \
        check_per_fg(mine, "this_study.tsv")

    print(f"prior_studies.tsv : {len(prior)} rows, "
          f"{prior['study'].nunique()} studies, "
          f"{prior['classifier'].nunique()} classifier variants")
    print(prior.groupby(["study", "classifier"]).size().to_string())
    print()
    _, notes = load_this_study(args.this_study_tsv, args.results_dir)
    for n in notes:
        print(f"  {n}")
    print()
    unverified = prior[~prior["verified"]]
    print(f"unverified rows: {len(unverified)}")
    for _, r in unverified.iterrows():
        print(f"  ! {r['study']} / {r['condition']} / {r['classifier']}")
    print()
    if problems:
        print("PER-FG CONSISTENCY PROBLEMS:")
        for p in problems:
            print("  " + p)
        sys.exit(1)
    print("per-fg consistency: OK (all rows satisfy count/(dna_pg*1000))")


if __name__ == "__main__":
    main()
