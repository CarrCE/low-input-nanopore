# Low-input nanopore performance vs prior work

Code in `bin/comparison/`, inputs in `assets/comparison/`, outputs written to
`results/comparison/`.

Every value is held in versioned text tables that cite the published table or
the reanalysis each row rests on, and the figure regenerates from them. An
earlier revision of this comparison lived in a hand-maintained spreadsheet; it
is not in the repository, and nothing here depends on it.

The figure places this study's per-femtogram sequencing performance against
Mojarro et al. 2019, Basapathi Raghavendra et al. 2023 and Zorzano et al. 2025
on log-log axes: **reads / fg DNA into library prep** (x) against
**bases / fg DNA into library prep** (y).

## Contents

| File | Role |
|---|---|
| `prior_studies.tsv` | One row per prior-study sample **per classifier variant**. Every row carries explicit provenance. Committed; never regenerated at runtime. |
| `this_study.tsv` | Round 1 / Round 2 headline values, schema-compatible with `prior_studies.tsv`. A committed snapshot of pipeline output (`source=pipeline_run_snapshot`, written by `seed_this_study.py`) so the module runs standalone and plots the published values. |
| `comparison_data.py` | Loader: reads the TSVs, selects one classifier per study, and optionally supersedes the seeded this-study rows with live pipeline output. Run directly for a data-integrity check. |
| `plot_comparison.py` | Regenerates the figure and its display-item outputs. |
| `seed_this_study.py` | Regenerates `this_study.tsv` as a snapshot of pipeline output, so the figure can be redrawn without the reads and the seeded and live paths cannot disagree. |
| `kraken2_db.manifest.tsv` | The pinned classification databases behind every `kraken2_*` row: URL, size, checksum and where the checksum came from. Authoritative; both scripts below read it. |
| `fetch_kraken2_db.sh` | Downloads and verifies those databases. Idempotent; refuses to proceed on any size or checksum mismatch. |
| `run_kraken2_reanalysis.sh` | Headless `wf-metagenomics` run pinned to v2.14.1, against the verified local databases. The command-line equivalent of the EPI2ME desktop run that produced the committed numbers. |
| `fetch_raghavendra.sh` | Stages the 10 deposited Basapathi Raghavendra FASTQs from Zenodo. |
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

The defective combination itself — x from `published_squeezemeta`, y from
`kraken2_q1` — is no longer carried as rows of its own. It was, briefly, so the
earlier figure could be reproduced for audit; but a mixture of two variants that
are both present adds no information, and 9 rows that must not be plotted are a
hazard out of proportion to their use. Pair the two variants on `condition` if
you ever need it back.

The mixing changed the conclusion, which is why it is documented at all. It
moved the Zorzano points horizontally only — the y values were already Kraken2.
Under the hybrid the largest Zorzano improvement was the **Blank** (7,659x
reads); under `kraken2_q1` it is the **Microbialite (2800 a)** (4,090x reads,
against 614x for that same condition under the hybrid). The figure's
auto-selected callout follows.

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
| `source` | Coarse provenance class: `published_table`, `kraken2_reanalysis`, `minimap2_reanalysis_this_pipeline`, `legacy_spreadsheet`, `pipeline_run`, `pipeline_run_snapshot`. |
| `source_detail` | Exact provenance: workbook sheet and cell references, accessions, tool versions, URLs. |
| `provenance_note` | Prose explanation, including any defect flags. |

`this_study.tsv` adds `round` (`Round 1` / `Round 2`, the figure series),
`experiment`, `sample_id` and `mode`.

## Regenerating the figure

```bash
cd comparison
python3 plot_comparison.py
```

Writes to `results/comparison/`:

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

### Refreshing the this-study snapshot

`this_study.tsv` is a committed snapshot of pipeline output, so the figure can be
redrawn without the reads. Refresh it whenever the pipeline output changes:

```bash
python3 bin/comparison/seed_this_study.py --results-dir results
```

It asks the loader for exactly the rows `--results-dir` would produce, so the
seeded and live paths cannot disagree, and commits nothing that a reader cannot
trace. Earlier versions of this file were extracted from a hand-maintained
spreadsheet; that spreadsheet is not in the repository, and the extractor that
read it has been removed along with the last reference to it.

## Classification provenance: what the `kraken2_*` rows ran against

Every row in `prior_studies.tsv` with `classifier=kraken2_q1` or
`kraken2_q10` is a reanalysis of published raw reads, and a reanalysis is only
as reproducible as the database it classified against. Kraken2 assigns a read to
the lowest common ancestor of its k-mer hits, so the database build determines
both what can be found and where in the taxonomy it lands. Two builds of "PlusPF"
are not interchangeable and swapping them changes counts.

The pin lives in **`kraken2_db.manifest.tsv`**:

| | |
|---|---|
| workflow | `epi2me-labs/wf-metagenomics` v2.14.1, commit `a57ff73c22b77c2754b7910cd8d24ab6056ed8cc` |
| classifier | kraken2 |
| database | Kraken2 PlusPF-8, **2024-12-28** build (`k2_pluspf_08gb_20241228.tar.gz`, 5,925,280,339 B, MD5 `01b8b1eb…`) |
| taxonomy | NCBI new_taxdump **2025-01-01** (139,761,991 B, MD5 `171470a1…`, SHA-256 `7ff98c65…`) |

Those two URLs are not a reconstruction. wf-metagenomics v2.14.1 hard-codes them
in its `database_sets` map under the key `PlusPF-8`, so `--database_set PlusPF-8`
on that revision resolves to exactly these files; the manifest records what the
reanalysis used rather than a plausible substitute.

One trap worth naming: NCBI serves the taxonomy at two paths. The one under
`pub/taxonomy/new_taxdump/new_taxdump.zip` is **rolling** — rewritten daily, and
today it returns a file with today's date. Only the dated copy under
`taxdump_archive/` is stable, and that is what is pinned. Substituting the
rolling path would silently reclassify against a different taxonomy.

```bash
bash bin/comparison/fetch_kraken2_db.sh                    # ~5.5 GiB, verified
bash bin/comparison/run_kraken2_reanalysis.sh \
     --fastq data/raghavendra_2023 --out results_kraken2 --min-qual 10
```

`fetch_kraken2_db.sh` stops on any size or checksum mismatch rather than
redownloading, because a mismatch means the URL has started serving something
other than the pinned build — the exact event this pin exists to catch. The
runner additionally re-resolves the workflow tag and stops if it no longer points
at the recorded commit; upstream tags can move.

**What this does and does not claim.** The committed `kraken2_q1` and
`kraken2_q10` values were produced through the EPI2ME desktop application, which
records its parameters only inside the run directory it creates. The scripts here
are the command-line equivalent with the same pinned inputs, added so the
provenance lives in version control; the committed numbers have **not** been
regenerated through them. Anyone re-running should expect agreement and should
report it if not.

Note also that `--min_read_qual` in wf-metagenomics thresholds on the mean
**Phred** value of a read, while this repository's own filter averages in
error-probability space per the ONT convention. They are different thresholds
with the same name; the value is passed through unchanged so that a re-run
reproduces the committed rows rather than a corrected version of them.

## Requirements

Python 3.12+, and only `numpy`, `pandas`, `matplotlib` (plus `openpyxl` for the
build-time extractor) — all present in the repo's analysis container. All
scripts derive their paths from `Path(__file__).resolve().parent` and accept
`--outdir`, so the module can be run from any working directory or relocated
wholesale.

The two shell scripts need more: `fetch_kraken2_db.sh` wants curl, tar, unzip and
one of md5sum/md5, plus ~13 GiB of free space while unpacking;
`run_kraken2_reanalysis.sh` additionally wants nextflow, docker and git, and at
least 8 GiB of RAM available to the workflow (the PlusPF-8 index is memory-mapped
unless `--kraken2_memory_mapping` is set). Neither is needed to regenerate the
figure — the committed TSVs cover that.

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

## Mapping-based re-analysis of Basapathi Raghavendra et al. 2023

### Why this exists

Defect (c) above compares two Raghavendra variants — the paper's own WIMP /
Centrifuge calls (`published`) and our Kraken2 minQ10 reanalysis (`kraken2_q10`)
— against *this* study's numbers, which come from competitive **minimap2
alignment against a reference set containing the exact organisms present**.
That is not a like-for-like comparison, and the asymmetry runs in our favour:

* Kraken2 against PlusPF-8 must decide, for every read, which of ~thousands of
  genomes (or none) a read belongs to, using exact k-mer matches. On noisy
  ONT reads at low depth this is conservative — reads that came from *E. coli*
  routinely fail to accumulate enough discriminating k-mers and land in the
  unclassified bin.
* Direct alignment to a two-genome index asks a much easier question: does this
  read align to *E. coli* or to *S. cerevisiae*? minimap2 tolerates the ~5-10%
  ONT error rate that k-mer classification does not, and there is no
  possibility of the read being assigned to some third organism.

So comparing our minimap2 numbers against their Kraken2 numbers inflates our
apparent advantage by an unknown factor that has nothing to do with library
chemistry or adaptive sampling — the thing we actually claim. Running their
deposited reads through this pipeline, with the same assignment rule, the same
`bases` convention and the same `reads / (dna_pg * 1000)` denominator, removes
that confound. Whatever advantage survives is real.

Note that a mapping re-analysis is only possible for Raghavendra et al. because
the organisms are known a priori (two defined species at stated masses). The
same treatment cannot be applied to Zorzano et al. 2025, whose samples are
environmental metagenomes with no known reference set — the Zorzano comparison
remains classifier-to-classifier and that limitation stands.

### Carrier: there is none

**Raghavendra et al. used no genomic carrier.** Methods, "Library preparation":

> The gDNA for the nanopore library was prepared using nuclease-free water and
> the low concentration of the DNA to be tested was directly pipetted from the
> previously diluted samples, into DNA LoBind sterile eppendorf tubes. The whole
> amount of DNA (either from *E. coli*, *S. cerevisiae*, the mock community DNA,
> or obtained from each MMS-2 extraction), was diluted in 25 µl for library
> preparation.

SQK-LSK114 at halved reagent volumes on an R10.4.1 Flongle, applied directly to
the diluted sample. No lambda and no other filler DNA appears anywhere in the
paper. The reference sets therefore contain **no `role=carrier` row**, and every
samplesheet row has `carrier_dna_ng=0`.

Consequence for the metrics: `input_sample_fraction = library / (library +
carrier) = 1.0`, so `enrichment` is identically equal to
`output_sample_fraction` and carries no information for this study. **Only
`reads_per_fg` and `bases_per_fg` are comparable** between Raghavendra et al.
and this study; the enrichment axis is ours alone, because only we add a
carrier to deplete against.

### Strains and accessions

| Organism | Accession | Basis |
|---|---|---|
| *Escherichia coli* | `GCA_900706755.1` | **Exact strain.** Methods, "Test samples": "*Escherichia coli* (NCTC 9001 Lenticule disc, Sigma Aldrich, UK) with an average genomic size of ~ 5 Mb". NCTC 9001 was sequenced by the Sanger Institute's NCTC 3000 project (BioProject PRJEB6403, BioSample SAMEA2517362, assembly 27731_A01): 5,040,580 bp in 5 contigs, contig N50 1,476,437 bp. The 5.04 Mb total corroborates the paper's stated "~5 Mb". |
| *Saccharomyces cerevisiae* | `GCF_000146045.2` | **Assumption, documented as such.** Methods: "yeast from *Saccharomyces cerevisiae*, Type II (YSC-2, 51,475, Sigma Aldrich) with an average genomic size of ~ 12 Mb". Sigma YSC-2 is a dried commercial yeast preparation sold by enzymatic activity, not a strain-defined culture-collection deposit — neither the supplier nor the paper gives a strain, and no assembly is attributable to it. We use the standard reference genome (strain S288C, assembly R64, complete, 12,071,326 bp), which matches the stated "~12 Mb". |

Two points worth knowing before editing these files:

* The *E. coli* accession is a **GenBank (`GCA_`) accession on purpose**. The
  paired RefSeq assembly `GCF_900706755.1` is **suppressed**, so
  `datasets download genome accession GCF_900706755.1` fails. `main.nf` and
  `bin/build_reference_set.py` both accept `GC[AF]_`, and
  `assets/references/lowinput_s1.tsv` already sets the precedent of using a
  contig-level `GCA_` assembly when it is the exact strain.
* The yeast strain assumption is low-risk here: the only competitor in the index
  is a bacterium, so strain-level divergence within *S. cerevisiae* cannot
  reassign reads between the two organisms. It can only change how many reads
  align at all.

*Homo sapiens* is deliberately **not** included as a contaminant reference, even
though the paper reports human reads in every sample and in the negative
controls (Table 2: 0-2 reads; Table 3: 1, 7, 19 and 53 reads). Two reasons.
(1) Cost: references are fetched and indexed once per reference set, and this
study needs four, so GRCh38 would be downloaded and indexed four times.
(2) Direction of the error: adding a competitor to a *competitive* assignment
can only move reads away from *E. coli* / *S. cerevisiae*, never toward them, so
omitting human makes these counts an **upper bound** on the prior study — the
conservative direction with respect to any advantage we claim. Unmapped human
reads fall into the `unassigned` bin, which `bin/compute_metrics.py` reconciles
explicitly, so they stay visible. Each reference-set preamble carries the exact
row to add if the stricter accounting is wanted; the Mix2 samples are where it
would matter most.

### Files

| File | Contents |
|---|---|
| `bin/comparison/fetch_raghavendra.sh` | Downloads `MinION low detectability.zip` from <https://zenodo.org/records/8208597>, verifies size and MD5, extracts the 10 `fastq_pass` files and stages them as `data/raghavendra_2023/<alias>.fastq.gz`. |
| `assets/samplesheets/raghavendra_2023.csv` | 10 rows, `include_in_headline=0` throughout. |
| `assets/references/raghavendra_2023_ecoli.tsv` | *E. coli* only, fraction 1.0 (Ec_R1-R3). |
| `assets/references/raghavendra_2023_scerevisiae.tsv` | *S. cerevisiae* only, fraction 1.0 (YSC_R1-R3). |
| `assets/references/raghavendra_2023_mix1.tsv` | 10 pg *E. coli* (0.833333) + 2 pg *S. cerevisiae* (0.166667). |
| `assets/references/raghavendra_2023_mix2.tsv` | 10 pg *S. cerevisiae* (0.833333) + 2 pg *E. coli* (0.166667). |

Four reference sets are needed rather than two because
`theoretical_dna_fraction` is what turns a read count into a per-fg number, and
it differs by sample design. Mix1 and Mix2 are mirror images — same two
organisms, fractions swapped — and are kept as separate files precisely so the
fractions cannot be attached to the wrong organism.

`library_dna_ng` is the **total** sample DNA into library prep: 0.010 ng for the
single-organism libraries and 0.012 ng (10 pg + 2 pg) for the mixes. Combined
with the fractions above this reproduces the per-organism masses already in
`prior_studies.tsv` — 10,000 fg and 2,000 fg — to within 4 parts per million,
so the new rows and the existing `kraken2_q10` / `published` rows share a
denominator and are directly comparable. These masses are the paper's own
stated dilution targets, not measurements: "The Qubit 4.0. fluorometer
sensitivity is limited to 10 pg/µl", so nothing at this level was quantified.

### How to run

```bash
bash bin/comparison/fetch_raghavendra.sh
./run.sh -profile docker --samplesheet assets/samplesheets/raghavendra_2023.csv
```

The fetch script is idempotent — a second run with all 10 FASTQs already
present and gzip-valid is a no-op and does not touch the network. It caches the
archive under `data/raghavendra_2023/.archive/` while working and deletes it
afterwards; set `KEEP_ARCHIVE=1` to keep it, or `ZIP=/path/to/archive.zip` to
extract from a copy you already have instead of downloading 771 MiB.

The run itself is small: the 10 deposited `fastq_pass` files hold between 2 and
520 reads each.

### Feeding the results back into `prior_studies.tsv`

The resulting rows enter `prior_studies.tsv` with
**`classifier=minimap2_competitive`**, `source=pipeline_run`, alongside — not
instead of — the existing `kraken2_q10` and `published` rows for the same
samples. **All three variants are retained.** They answer three different
questions and none supersedes the others:

| `classifier` | What it is | Has `bases`? |
|---|---|---|
| `published` | The paper's own WIMP / Centrifuge pass-read counts, Tables 2 and 3. | No — the paper publishes reads only, so these rows cannot be plotted on the y axis. |
| `kraken2_q10` | Our Kraken2 / PlusPF-8 minQ10 reanalysis of the deposited reads. Currently the plot default. | Yes. |
| `minimap2_competitive` | This pipeline, same assignment rule and same `bases` convention as our own data. The only variant that is methodologically like-for-like with `this_study.tsv`. | Yes. |

Keeping all three makes the size of the classifier effect measurable: the gap
between `kraken2_q10` and `minimap2_competitive` on identical reads *is* the
mapping-versus-classification advantage, quantified rather than assumed, and it
is the correct amount to discount our headline comparison by.

### Two provenance findings from staging the data

Recorded here because both bear on rows already in `prior_studies.tsv`:

1. **The "263" typo is confirmed.** `prior_studies.tsv` suspected that Table 3's
   263 *S. cerevisiae* pass reads for `APN068` was a transcription of that
   file's total read count. The deposited `APN068_pass_...fastq.gz` contains
   **exactly 263 reads**. The suspicion is now a verified fact.
2. **The "Sample total reads in the deposited file" notes are mislabelled.**
   Those `provenance_note` values (224, 73, 313, 337, 37, 190, 411, 1180, 2110,
   4700) are the paper's MinKNOW **pass + fail** totals from Tables 2 and 3, not
   the read counts of the deposited `fastq_pass` files, which are 8, 4, 27, 23,
   8, 2, 7, 455, 520 and 263 respectively. The counts themselves are correct and
   correctly attributed to the paper; only the phrase "in the deposited file" is
   wrong. Both numbers are now carried in the `notes` column of
   `assets/samplesheets/raghavendra_2023.csv`.

### Author Correction: does not affect these numbers

There is an Author Correction to this paper — Sci Rep **15**:14107 (2025),
<https://doi.org/10.1038/s41598-025-98000-4>, published 23 April 2025. It reads
in full:

> The original version of this Article contained errors in Figure 4 a and b,
> where the data points were interchanged. The negative control data points were
> therefore incorrectly represented as these were stated as non-zero, which
> contradicted the quantifications of the other conditions.

**Assessment: no effect on Table 2 or Table 3, and therefore none on any
Raghavendra row in `prior_studies.tsv`.** Figure 4a/4b are the DNA-yield curves
from the MMS-2 Mars-simulant regolith incubation experiment — a different
experiment, a different sample type, and a fluorometric mass measurement rather
than a read count. Tables 2 and 3, which are the sole source of the `published`
rows and the design of the 10 samples re-analysed here, are untouched. The
correction changes no read count, no input mass and no organism assignment. The
version consulted is the corrected one (its final page reads
"© The Author(s) 2023, corrected publication 2025"), and the fetched Zenodo
reads are unchanged by it.
