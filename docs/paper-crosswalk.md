# Paper crosswalk

Which output of this repository backs which display item in the paper.

Numbering is deliberately absent from the rest of the documentation: a
"Supplementary Table S11" that is renumbered in review turns every reference to
it into a silent error, and the repository has no way to notice. Outputs are
named by what they are and where they land. This file is the one place the two
are tied together, so there is exactly one thing to update when numbering
changes.

**The `paper` column is deliberately blank.** The paper of record is the first
author's own manuscript and its numbering is not final. Fill the column in when
it is; until then the `draft` column records where each output appears in the
AI-drafted manuscript that is being used for number-checking, which is *not* the
document being submitted.

## Figures

| output | what it shows | draft | paper |
|---|---|---|---|
| `results/summary/abundance.pdf` + `.csv` + `.json` | theoretical vs measured community abundance | main text, abundance figure | |
| `results/summary/readlengths.pdf` + `.csv` + `.json` | read-length distribution by assignment class | main text, read-length figure | |
| `results/comparison/low_input_comparison.pdf` + `.csv` + `.json` | performance per femtogram against prior studies | main text, per-femtogram figure | |
| `results/summary/coverage.pdf` + `.csv` + `.json` | coverage uniformity, depth in 1 kb bins | SI, coverage-uniformity figure | |
| `results/summary/pooled_coverage.pdf` + `.csv` + `.json` | coverage uniformity pooled across replicates | SI, pooled-coverage figure | |
| `results/summary/mode_delta.pdf` + `.csv` + `.json` | reads retained, competitive vs sequential | SI, mode-comparison figure | |

Each figure ships a `.csv` of **the points it draws** and a `.json` carrying its
caption, source files, software versions and the derived quantities quoted in
the text. A number in a figure can be audited from the CSV without rerunning
anything. See [`display-items.md`](display-items.md).

## Tables

| output | what it shows | draft | paper |
|---|---|---|---|
| `results/summary/experiment_summary.tsv` | experiment-level summary, both decision rules | SI, experiment-summary table | |
| `results/summary/per_sample.tsv` | per-replicate results | main text and SI, per-replicate tables | |
| `results/summary/per_organism.tsv` | per-organism reads, bases and per-femtogram yields | source for the abundance figure and spike-in tables | |
| `results/summary/coverage_attribution.tsv` | coverage summary for every organism-replicate pair reaching 1× | SI, coverage table | |
| `results/summary/mode_delta_summary.csv` | community reads destroyed by sequential subtraction | SI, mode-delta table | |
| `results/q10/summary/experiment_summary.tsv` | quality-matched rerun at Q10 | SI, Q10 table | |
| `results/contaminant_divergence/summary.tsv` | divergence of the carrier-derived *E. coli* from K-12 MG1655 | SI, contaminant-reference table | |
| `assets/references/*.tsv` | reference sets: every organism, its role and accession | SI, reference-set table | |
| `assets/measurements.tsv` | sequencing runs and the DNA masses behind every per-femtogram and enrichment figure | SI, runs table | |
| `assets/comparison/prior_studies.tsv` | prior-study values used in the comparison figure | source for the per-femtogram figure | |

## Four tables that need an extra target

`results-of-record/` was produced by `make all`, `q10`, `raghavendra`,
`divergence`, the display items, `verify` and `check`. These four SI tables come
from targets outside that set, so their outputs are **not** in that tree and
must be regenerated deliberately:

| output | target | note |
|---|---|---|
| `results/summary/human_masking.tsv` | `make masksummary` | reads `results/*/human/*.human_stats.json`, which only a `--mask_human` run writes. The published table came from the 13 Aug 2026 masking run. |
| `results/summary/sequencing_summary.tsv` | `make seqsummary` | needs the FASTQs as well as the assignments. Note that reads fetched from the SRA archive have no `qs:f:` header tag, so `median_qscore` comes back empty from a fetched copy. |
| `results/summary/attribution_threshold.*` | `make threshold` | the distribution behind `--min_aln_frac`; regenerate this rather than quoting the threshold from memory if it is ever questioned |
| software versions table | `make versions` | what the built images report when asked, not what was pinned |

That gap is stated rather than papered over: a reader who runs the documented
sequence and then looks for `human_masking.tsv` should know why it is not there.
