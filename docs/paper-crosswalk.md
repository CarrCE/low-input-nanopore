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
document being submitted. Draft numbering as of 2 Nov 2026 (SI restructured
into sections S1-S9).

## Figures

| output | what it shows | draft | paper |
|---|---|---|---|
| `results/summary/mode_delta.pdf` + `.csv` + `.json` | reads retained, competitive vs sequential | Fig. 2 | |
| `results/summary/abundance.pdf` + `.csv` + `.json` | theoretical vs measured community abundance | Fig. 3 | |
| `results/summary/readlengths.pdf` + `.csv` + `.json` | (a) read-length distribution by assignment class; (b) carrier read lengths per experiment, by read and by base | Fig. 4 | |
| `results/comparison/low_input_comparison.pdf` + `.csv` + `.json` | performance per femtogram against prior studies | Fig. 5 | |
| `results/summary/carrier-model.pdf` + `.csv` + `.json` | carrier-depletion model against the runs (`make carriermodel`) | Fig. S1 | |
| `results/summary/coverage.pdf` + `.csv` + `.json` | coverage uniformity, depth in 1 kb bins, per library | Fig. S2 | |
| `results/summary/pooled_coverage.pdf` + `.csv` + `.json` | coverage uniformity pooled across replicates | Fig. S3 | |

Fig. 1 (study design) is drawn in the manuscript and has no output here.

Each figure ships a `.csv` of **the points it draws** and a `.json` carrying its
caption, source files, software versions and the derived quantities quoted in
the text. A number in a figure can be audited from the CSV without rerunning
anything. See [`display-items.md`](display-items.md).

## Tables

| output | what it shows | draft | paper |
|---|---|---|---|
| `assets/measurements.tsv` | sequencing runs and the DNA masses behind every per-femtogram and enrichment figure | Table S1 | |
| `assets/references/*.tsv` | reference sets: every organism, its role and accession | Table S2 | |
| `results/summary/attribution_threshold.*` | the distribution behind `--min_aln_frac` (`make threshold`, from the floor-off run) | Table S3 | |
| `results/summary/floor_effect.json`, `floor_effect_by_organism.tsv` | every number in the attribution-floor paragraph: what the floor excludes, from whom, Q10's overlap, and the effect on each experiment (`make floor`) | Section S2 text | |
| `results/summary/floor_complexity.json`, `floor_complexity_references.tsv` | what the excluded reads are by sequence, and which reference genomes hold simple repeats (`make floorcomplexity`) | Section S2 text | |
| `results/summary/experiment_summary.tsv` + `results/q10/summary/experiment_summary.tsv` | experiment-level summary, both decision rules, all reads and Q10 | Table S4 | |
| `results/summary/human_masking.tsv` | human read masking (`make masksummary`) | Table S5 | |
| `results/summary/mode_delta_summary.csv` | community reads destroyed by sequential subtraction | Table S6 | |
| `results/contaminant_divergence/summary.tsv` | divergence of the carrier-derived *E. coli* from K-12 MG1655 | Table S7 | |
| `results/summary/sequencing_summary.tsv` | yield, read length and quality per replicate (`make seqsummary`) | Table S8 | |
| `results/summary/per_sample.tsv` | per-replicate results | Table 1 (D6321); Table S9 (D6311) | |
| `results/summary/readlengths_summary.csv` | read-length statistics by assignment class | Table S10 | |
| `results/summary/carrier-model-table.csv` + `.json` | carrier-depletion model, per library (`make carriermodel`) | Table S11 | |
| `results/summary/coverage_attribution.tsv` | coverage summary for every organism-replicate pair reaching 1× | Table S12 | |
| `results/<sample>/coverage/` (`bin/coverage_dropouts.py`) | *Listeria* coverage dropouts | Table S13 | |
| `assets/comparison/prior_studies.tsv` | prior-study values; the Basapathi Raghavendra et al. mixtures | Table S14; source for Fig. 5 | |
| software versions table | what the built images report (`make versions`) | Table S15 | |
| `results/summary/per_organism.tsv` | per-organism reads, bases and per-femtogram yields | source for Fig. 3 and the spike-in values | |
| `results/comparison/estimated_control.csv` + `.json` | estimated no-adaptive-sampling control | Section S6 text | |
| `results/carrier_model/*.tsv` | carrier-model analysis tables ([`carrier-model.md`](carrier-model.md)) | Section S6 text | |

## Tables that need an extra target

`results-of-record/` was produced by `make all`, `q10`, `raghavendra`,
`divergence`, the display items, `verify` and `check`. These SI items come
from targets outside that set, so their outputs are **not** in that tree and
must be regenerated deliberately:

| output | target | note |
|---|---|---|
| `results/summary/human_masking.tsv` | `make masksummary` | reads `results/*/human/*.human_stats.json`, which only a `--mask_human` run writes. The published table came from the 13 Aug 2026 masking run. |
| `results/summary/sequencing_summary.tsv` | `make seqsummary` | needs the FASTQs as well as the assignments. Note that reads fetched from the SRA archive have no `qs:f:` header tag, so `median_qscore` comes back empty from a fetched copy. |
| `results/summary/attribution_threshold.*` | `make threshold` | the distribution behind `--min_aln_frac`, read from the floor-off run `results/nofloor` (`make nofloor`): in the published run the reads below the floor are already unassigned. Regenerate this rather than quoting the threshold from memory. |
| `results/summary/floor_effect.*` | `make floor` | needs `make all`, `make nofloor` and `make q10`. Counts the reads the assigner actually excluded, which is slightly stricter than `attribution_threshold.py`'s view (winning organism only). The paper's floor paragraph was first computed on the reads before human masking; regenerate it here, on the deposited reads. |
| `results/fragmentation/summary/fragmentation.json` | `make fragmentation` | needs `make all` and `make nofloor`. Backs what the paper says about the cause of the floor's eukaryote losses: genome content, not assembly quality; see [`reference-fragmentation.md`](reference-fragmentation.md). About 1 h 12 m. |
| `results/summary/floor_complexity.json`, `floor_complexity_references.tsv` | `make floorcomplexity` | needs `make all` and `make nofloor`, and the deposited (masked) reads, never the originals in `data/`. The SI S2 statements on what the excluded reads are (masked share, compressibility, low-entropy windows, 1-3 bp tandem repeats) and on which reference genomes hold simple repeats. |
| software versions table | `make versions` | what the built images report when asked, not what was pinned |
| `results/summary/carrier-model*`, `results/carrier_model/` | `make carriermodel` | needs the masked reads and the MinKNOW instrument-data record as well as a finished run; see [`carrier-model.md`](carrier-model.md) |

That gap is stated rather than papered over: a reader who runs the documented
sequence and then looks for `human_masking.tsv` should know why it is not there.
