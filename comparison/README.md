# `comparison/` — low-input nanopore performance vs prior work

This module replaces a hand-maintained Excel workbook
(`2026-04-29 Low-Input Comparison Calculations.xlsx`) with versioned, auditable
text tables plus a script that regenerates the comparison figure.

The figure places this study's per-femtogram sequencing performance against
Mojarro et al. 2019, Basapathi Raghavendra et al. 2023 and Zorzano et al. 2025
on log-log axes: **reads / fg DNA into library prep** (x) against
**bases / fg DNA into library prep** (y).

## Contents

| File | Role |
|---|---|
| `prior_studies.tsv` | One row per prior-study sample **per classifier variant**. Every row carries explicit provenance. Committed; never regenerated at runtime. |
| `this_study.tsv` | Round 1 / Round 2 headline values, schema-compatible with `prior_studies.tsv`. Seeded from the workbook (`source=legacy_spreadsheet`) so the module runs standalone. |
| `comparison_data.py` | Loader: reads the TSVs, selects one classifier per study, and optionally supersedes the seeded this-study rows with live pipeline output. Run directly for a data-integrity check. |
| `plot_comparison.py` | Regenerates the figure and its display-item outputs. |
| `extract_workbook.py` | **Build-time only.** The one-time xlsx -> TSV extraction, committed so the derivation is auditable. Nothing in the runtime path opens the workbook. |
| `figures/` | Generated outputs (PDF / PNG / CSV / JSON). |

## Formulas

Every per-fg value in both TSVs, for every study including this one, uses:

```
reads_per_fg = reads / (dna_pg * 1000)
bases_per_fg = bases / (dna_pg * 1000)
```

`dna_pg` is **DNA into library prep** (i.e. post-extraction), never DNA into
extraction and never cell-equivalent mass. `bases` is the full read length of
every classified read, not the aligned span — this matches how the prior
studies count Kraken2-classified read lengths (see the metric definitions at
the top of `bin/compute_metrics.py`).

Where a prior study's DNA mass was below the fluorometer detection limit
(Zorzano et al. 2025, 7 of 9 conditions), `dna_pg` is back-calculated from that
paper's own main-text yield of 1.046 Mbp per ng extracted DNA:
`dna_ng = total_bases / 1.046e6`. The `dna_basis` column records which of the
two applies to each row.

The dashed contours on the figure are iso-improvement circles around the
Round 1 mean: geometric-mean fold improvement expressed as Euclidean distance
in log10-log10 space, at 10x / 100x / 1000x.

## Known defects, and how each is handled

These were all present in the workbook. None of them is silently repaired; each
is represented explicitly in the data and surfaced by the plot script.

### (a) Zorzano et al. 2025 — the two axes used different classifiers

In the workbook the x value came from `'Zorzano et al. 2025'!Q`, which is the
paper's own SqueezeMeta **"Reported Hits"** (column G) divided by DNA mass, while
the y value came from `'Zorzano et al. 2025'!R`, which is **our Kraken2 minQ1
hit-bases** (column I) divided by the same mass. Those are two different
classifiers, and they disagree substantially:

| Condition | Reported hits (SqueezeMeta) | Kraken2 minQ1 hits |
|---|---:|---:|
| Blank | 28 | 213 |
| Microbialite (2800 a) | 7,806 | 1,186 |
| Microbialites (2800 a), 0.86 MGy | 1,150 | 497 |
| Microbialites (2800 a), 10.45 MGy | 938 | 691 |
| Stromatolites (541 Ma) | 858 | 683 |
| Stromatolites (541 Ma), 10.45 MGy | 697 | 653 |
| Oxide iron formation (2930 Ma) | 254 | 663 |
| Oxide iron formation (2930 Ma), 10.45 MGy | 59 | 124 |
| Carbonate formation (2930 Ma) | 406 | 645 |

The disagreement is not even one-directional: SqueezeMeta reports 6.6x more hits
than Kraken2 for the Microbialite, and 2.6x fewer for the Oxide iron formation.
Mixing them means the plotted Zorzano points do not correspond to any single
analysis.

`prior_studies.tsv` therefore carries **three** variants of each Zorzano
condition, distinguished by the `classifier` column, and the plot script selects
one with `--zorzano-classifier`:

| `classifier` | x from | y from | Notes |
|---|---|---|---|
| `kraken2_q1` **(default)** | Kraken2 minQ1 hits | Kraken2 minQ1 hit bases | Internally consistent. Both axes from the same reanalysis of the published raw reads (wf-metagenomics v2.14.1, kraken2, PlusPF-8, minQ 1). |
| `published_squeezemeta` | Reported Hits | hits x mean read length | Consistent in the sense that both axes derive from the paper's own hit calls, but the paper never publishes hit **bases**, so y is a derived estimate assuming hits have the same length distribution as non-hits. That assumption is demonstrably wrong — our Kraken2 minQ1 hits average 2.7x to 21.5x the overall mean read length, depending on condition — so treat this y as a rough sketch only. |
| `legacy_hybrid_workbook` | Reported Hits | Kraken2 minQ1 hit bases | The defective combination, kept **only** so the previous figure can be reproduced for audit. Verified to reproduce the old script's 24 plotted points bit-for-bit. Do not publish. |

Switching from `legacy_hybrid_workbook` to the default `kraken2_q1` moves the
Zorzano points horizontally only — the y values were already Kraken2 — but it
changes the conclusion. Under the workbook's hybrid the largest Zorzano
improvement was the **Blank** (7,659x reads); under `kraken2_q1` it is the
**Microbialite (2800 a)** (4,090x reads, up from the workbook's 614x for that
same condition). The auto-selected figure callout follows.

### (b) Mojarro et al. 2019 — the counts have no recorded citation

The workbook's `'Mojarro et al. 2019'` sheet contains, as hand-typed literals
with the source/URL columns left blank:

```
Reads obtained    5
Bases obtained    5270
```

Only the DNA mass (2 pg) carries an attribution. The sheet's own
first-principles estimate — 10,000 *B. spizizenii* cells x 4.4276 fg DNA/cell x
5% simulated extraction efficiency = 2.214 pg — corroborates the **mass**, but
nothing anywhere in the file supports the read or base counts.

Handling:

* The row is in `prior_studies.tsv` with `verified=FALSE`,
  `source=unsourced_literal`, `classifier=published_unverified`, and a
  `provenance_note` flagging it as UNVERIFIED / citation needed.
* `plot_comparison.py` prints a loud banner to stderr whenever an unverified
  row is plotted, naming the row and echoing the note.
* `--drop-unverified` removes it from the figure entirely.

**This must be resolved before submission**: either trace 5 reads / 5,270 bases
to a specific table or figure in Mojarro et al. 2019 and update
`source_detail`/`verified`, or drop the point.

### (c) Basapathi Raghavendra et al. 2023 — reanalysis disagrees with the paper

Our Kraken2 minQ10 reanalysis of the authors' own deposited reads
(<https://zenodo.org/records/8208597>) does not reproduce the pass-read counts
in the paper's Tables 2 and 3:

| Condition | Rep / species | Published pass reads | Kraken2 minQ10 |
|---|---|---:|---:|
| 10 pg *E. coli* lenticule discs | 1 / *E. coli* | 4 | 4 |
| 10 pg *E. coli* lenticule discs | 2 / *E. coli* | 3 | 2 |
| 10 pg *E. coli* lenticule discs | 3 / *E. coli* | 1 | 13 |
| 10 pg YSC-2 yeast | 1 / *S. cerevisiae* | 4 | 11 |
| 10 pg YSC-2 yeast | 2 / *S. cerevisiae* | 4 | 2 |
| 10 pg YSC-2 yeast | 3 / *S. cerevisiae* | 2 | 2 |
| 10 pg *E. coli* + 2 pg *S. cerevisiae* | 1 / *E. coli* | 2 | 2 |
| 10 pg *E. coli* + 2 pg *S. cerevisiae* | 1 / *S. cerevisiae* | 2 | 2 |
| 10 pg *E. coli* + 2 pg *S. cerevisiae* | 2 / *E. coli* | 120 | 32 |
| 10 pg *E. coli* + 2 pg *S. cerevisiae* | 2 / *S. cerevisiae* | 17 | 157 |
| 10 pg *S. cerevisiae* + 2 pg *E. coli* | 1 / *E. coli* | 88 | 29 |
| 10 pg *S. cerevisiae* + 2 pg *E. coli* | 1 / *S. cerevisiae* | 243 | 136 |
| 10 pg *S. cerevisiae* + 2 pg *E. coli* | 2 / *E. coli* | 78 | 35 |
| 10 pg *S. cerevisiae* + 2 pg *E. coli* | 2 / *S. cerevisiae* | 263 | 10 |

Both are in `prior_studies.tsv`, as `classifier=kraken2_q10` and
`classifier=published`, each with its own `source` and `provenance_note`.
`--raghavendra-classifier` selects between them.

The default is `kraken2_q10` for a hard reason, not a preference: **the paper
publishes pass reads only, no base counts.** The `published` rows therefore have
`bases` and `bases_per_fg` empty, cannot be placed on the y axis, and the loader
excludes them with an explicit message — running
`--raghavendra-classifier published` produces a figure with no Raghavendra
points. That flag exists to make the limitation visible, not to produce an
alternative figure.

The last row is worth singling out: the paper's 263 *S. cerevisiae* pass reads
for `APN068` equals the **total** number of reads in that file, which strongly
suggests a transcription error rather than a classification disagreement. This
is recorded in that row's `provenance_note`.

### (d) Minor: a stale label in the workbook

`'Per Base Comparison'!B7` is labelled "All conditions (hits, assuming same
length as non-hits)", but its formula is `SUM('Zorzano et al. 2025'!I4:I12)` —
a direct sum of measured Kraken2 minQ1 hit bases, with no read-length
extrapolation involved. The label is wrong; the number (12,592,526 bases) is
right. Carried in `prior_studies.tsv` as a `row_type=aggregate` row with the
correction in its `provenance_note`. Aggregate rows are never plotted.

### (e) Open item: Round 2 input mass does not match the samplesheet

The seeded Round 2 value comes from the workbook's `S2` sheet: 3 D6321 preps
pooled into one extraction, 10 uL elution, 2 uL to Qubit, leaving
8 uL x 0.0279 ng/uL = **0.2232 ng into library prep**, giving 7,257 reads and
18,899,147 bases. `assets/samplesheets/lowinput_s2.csv` instead lists three
replicates at 0.26 / 0.32 / 0.40 ng. The legacy value does not correspond 1:1 to
any of them. This is flagged in the row's `provenance_note` and is resolved the
moment live pipeline output is used (`--results-dir`), which supersedes the
seeded row with per-replicate values.

## Schema

Both TSVs share these columns:

| Column | Meaning |
|---|---|
| `study` | Full study name, used for colour/marker lookup. |
| `study_short` | Short form used in figure callouts. |
| `condition` | Sample/condition label as published. |
| `organism` | Target organism, or `metagenome` / `All organisms`. |
| `replicate_idx` | 0-based replicate index **within** `(study, condition)`. |
| `replicate_label` | Human-readable replicate description. |
| `row_type` | `sample` (plottable) or `aggregate` (whole-study summary, never plotted). |
| `classifier` | Which read/base assignment produced the counts. Exactly one variant per study is plotted. |
| `reads`, `bases` | Raw counts. Empty when the source does not publish them. |
| `dna_pg` | DNA into library prep, picograms. |
| `dna_basis` | How `dna_pg` was obtained (measured vs back-calculated). |
| `reads_per_fg`, `bases_per_fg` | Derived per the formulas above. |
| `verified` | `FALSE` when the counts have no traceable citation. Triggers the plot warning. |
| `source` | Coarse provenance class: `published_table`, `kraken2_reanalysis`, `unsourced_literal`, `legacy_spreadsheet`, `pipeline_run`. |
| `source_detail` | Exact provenance: workbook sheet and cell references, accessions, tool versions, URLs. |
| `provenance_note` | Prose explanation, including any defect flags. |

`this_study.tsv` adds `round` (`Round 1` / `Round 2`, the figure series),
`experiment`, `sample_id` and `mode`.

## Regenerating the figure

```bash
cd comparison
python3 plot_comparison.py
```

Writes to `comparison/figures/`:

* `low_input_comparison.pdf` — vector, Type 42 (editable) fonts
* `low_input_comparison.png` — 600 dpi
* `low_input_comparison.csv` — every plotted point with full provenance, plus
  fold-improvement columns against the Round 1 mean
* `low_input_comparison.json` — display-item sidecar with
  `id` / `title` / `caption` / `source_files` / `software` / `metrics`.
  `source_files` records a SHA-256 for each input, so a figure can always be
  tied back to the exact table that produced it.

Useful options (`--help` on every script lists them all):

```bash
# use live pipeline output instead of the seeded this-study values
python3 plot_comparison.py --results-dir ../results

# reproduce the old spreadsheet figure exactly, for audit
python3 plot_comparison.py --zorzano-classifier legacy_hybrid_workbook \
                           --basename low_input_comparison_legacy

# plot Zorzano on the paper's own hit calls for both axes
python3 plot_comparison.py --zorzano-classifier published_squeezemeta

# drop the uncited Mojarro point
python3 plot_comparison.py --drop-unverified

# write somewhere else
python3 plot_comparison.py --outdir /path/to/manuscript/figures
```

### Using live pipeline output

With `--results-dir`, the loader scans
`results/<sample_id>/<mode>/<sample_id>.metrics.tsv` and takes the row where
`organism == "All organisms"` as each replicate's headline value. Rows are
mapped to figure series by `experiment` (`lowinput_s1` -> Round 1,
`lowinput_s2` -> Round 2; other experiments, e.g. the smoke test, are ignored).
Replicates with no quantified input mass — `lowinput_s2_r0`, whose Qubit reading
was off-scale low and whose per-fg fields are therefore blank — are skipped.

Superseding is per experiment: if any live rows are found for `lowinput_s1`, all
seeded `lowinput_s1` rows are dropped. Seeded rows for experiments with no live
output are retained, and the `source` column in the output CSV shows which of
`legacy_spreadsheet` / `pipeline_run` each point came from.

### Data integrity check

```bash
python3 comparison_data.py
```

Prints the row counts per study and classifier, lists unverified rows, and
re-derives every `reads_per_fg` / `bases_per_fg` from `reads`, `bases` and
`dna_pg` to confirm the stored values are consistent. Exits non-zero on any
inconsistency.

### Re-extracting from the workbook

Only if the workbook itself is corrected:

```bash
python3 extract_workbook.py --workbook "/path/to/2026-04-29 Low-Input Comparison Calculations.xlsx"
```

This rewrites both TSVs and asserts that every emitted per-fg value reproduces
the workbook's own computed cell. Commit the regenerated TSVs. The workbook
lives under the gitignored `ignore/` tree, so this step is not reproducible from
a clean clone — which is the point of committing the TSVs.

## Requirements

Python 3.12+, and only `numpy`, `pandas`, `matplotlib` (plus `openpyxl` for the
build-time extractor) — all present in the repo's analysis container. All
scripts derive their paths from `Path(__file__).resolve().parent` and accept
`--outdir`, so the module can be run from any working directory or relocated
wholesale.

## Style notes

The figure keeps the visual language of the original manuscript script:
colourblind-safe Okabe-Ito palette (Mojarro `#0072B2` circles, Basapathi
Raghavendra `#D55E00` squares, Zorzano `#009E73` triangles), black stars for
Round 1 and a black diamond for Round 2, dashed iso-improvement contours,
per-study callouts naming the largest improvement, marginal rugs on the top and
right axes, a dotted constant-bases-per-read reference line, and italicised
species names via mathtext. Unlike the original, Round 1 and Round 2 marker
positions and label text are computed from the data rather than hardcoded, so
the figure follows the TSVs (or live results) automatically.
