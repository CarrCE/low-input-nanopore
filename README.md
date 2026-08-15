# low-input-nanopore

Depletion-mode adaptive sampling lets an Oxford Nanopore flow cell eject reads
that match a known carrier, so a tiny amount of sample DNA can be carried into
library prep on a large mass of cheap genomic carrier (here, lambda) and still
dominate the sequencing output. This repository quantifies *how much*: it takes
basecalled reads from libraries built as `sample DNA + ~1000 ng lambda carrier`,
assigns every read to exactly one organism — community member, carrier,
carrier-derived *E. coli* contaminant, ambiguous, or unassigned — and reports the
enrichment achieved (the sample's share of output bases divided by its share of
input mass) together with yield normalised to input DNA mass (reads and bases per
femtogram). Two datasets are analysed: `lowinput_s1` (ZymoBIOMICS Microbial
Community DNA Standard II, log distribution, D6311; 3 replicates) and
`lowinput_s2` (ZymoBIOMICS Spike-in Control II, low microbial load, D6321, whole
cells; 4 replicates, r0–r3).

### "Sub-picogram" means per organism, not per library

Worth stating up front, because the two scales are three orders of magnitude
apart and easy to conflate. The **libraries** are sub-nanogram: about 1.21 ng of
sample DNA for `lowinput_s1`. The **detections** are sub-picogram, and that is
the claim: in a log-distributed community, the rarest members enter library prep
at femtogram masses and are still recovered.

Both figures below are computed from `results/summary/per_organism.tsv`, not
asserted:

| organism | design fraction | mass into library prep | reads recovered |
|---|---|---|---|
| *Staphylococcus aureus* | 1×10⁻⁶ | **1.21 fg** | 3, 3, 4 — all three `lowinput_s1` replicates |
| *Cryptococcus neoformans* | 9×10⁻⁶ | **10.9 fg** | 4, 0, 9 — two of three |

*C. neoformans* is absent from one replicate, which is what detection at ten
femtograms looks like: a handful of reads, and Poisson sampling doing the rest.
Reporting it as detected in all three would overstate what these data show. The
per-organism yields these rest on are in
[Metric definitions](#metric-definitions).

---

## Quickstart

Requirements on the host: **Nextflow (>=23.10.0)** and **Docker**. Nothing else —
every tool the pipeline calls lives in a pinned container.

Hardware, for the shipped defaults: **48 GB of memory available to Docker and 14
cores**, because `MAP_COMPETITIVE` asks for them (`conf/base.config`). On a
smaller machine pass `--max_memory 16.GB --max_cpus 8`; the mapping is slower
but nothing else changes. Budget **~110 GB** for the seven FASTQs and a further
**~65 GB** for Nextflow's `work/` and `results/`.

Measured on an M3 Max with 14 cores to Docker: the analysis itself
(`make all`) takes **1 h 32 m**, and the plotting and check targets add about
ten minutes. Two further targets cost roughly an hour and a half each — `q10`,
a second full pass over every replicate, and `divergence`, which runs breseq
under amd64 emulation. **A complete reproduction is 4 h 44 m**; the per-target
breakdown is [below](#reproducing-the-published-analysis).

### Getting the reads

The sequencing reads are public: BioProject `PRJNA1513130`, runs
`SRR40180147`–`SRR40180153`, released 15 August 2026.

**ENA has not mirrored them yet, and `--fetch_from_sra` resolves through ENA**,
so for now download the seven FASTQs from NCBI and place them in `data/` —
hard-link or copy, never symlink to a directory outside the repository, for the
reason in [`docs/operating-notes.md`](docs/operating-notes.md#getting-reads-into-data-link-or-copy-never-symlink-outward).
Steps 2 and 3 then run from a fresh clone. Once ENA catches up,
`--fetch_from_sra` does this for you; see [Data
availability](#data-availability).

One result needs no reads at all — the per-femtogram comparison figure is built
from committed tables: `make images && make comparison`.

### The four commands

```bash
# 1. Build the two local images (minimap2/samtools/seqkit/datasets, and python)
make images

# 2. Smoke test: 40,000 reads through the whole pipeline, a couple of minutes
#    (needs data/test/test_s2.fastq -- see `make demo-data`)
./run.sh -profile docker,test

# 3. The published analysis: both experiments, both assignment modes
./run.sh -profile docker --samplesheet assets/samplesheets/all.csv --mode both

# 4. Everything downstream of the pipeline, then the checks (see below)
make attribution coverage modedelta comparison estcontrol
make assigneddepth poolcov
make seqsummary runmeta versions divergence
make measurements verify
```

`make test` and `make all` are shorthands for steps 2 and 3. Before a long run,
read [`docs/operating-notes.md`](docs/operating-notes.md) — a short list of
things that each cost somebody an hour.

### Reproducing the published analysis

Step 3 is **the** analysis: one Nextflow run, both experiments, both assignment
modes. It produces the abundance and read-length figures and the study-level tables
directly. Everything in step 4 is either a plotter that reads
`results/summary/`, or a check — no second pipeline is involved.

Two are ordered — `coverage` requires `attribution`, and `poolcov` requires
`assigneddepth` — and the Makefile declares both, so `make coverage` on a fresh
tree builds what it needs rather than failing on a missing input.

All times below are measured on a complete cold run from a fresh clone (M3 Max,
14 cores to Docker), not estimated.

| step | target | produces | cost |
|---|---|---|---|
| 3 | `all` | abundance and read-length figures, `results/summary/` tables | **1 h 32 m** |
| 4 | `attribution` → `coverage` | coverage-uniformity figure, coverage table | 4 s |
| 4 | `modedelta` | mode-comparison figure | 3 s |
| 4 | `comparison` | per-femtogram comparison figure | 2 s |
| 4 | `estcontrol` | estimated no-AS control | 2 s |
| 4 | `assigneddepth` → `poolcov` | pooled-coverage figure | 5 m |
| 4 | `seqsummary`, `runmeta`, `versions` | sequencing-summary, runs and software tables | 5 m (`runmeta` reads all ~110 GB) |
| 4 | `q10` | the **quality-matched rerun** table, at Q10 | **1 h 27 m** |
| 4 | `raghavendra` | the committed prior-study rows behind the comparison figure | 38 s |
| 4 | `divergence` | the contaminant-divergence table | **1 h 33 m** (breseq under amd64 emulation) |
| 4 | `measurements`, `verify` | the two checks | seconds |
| | | **total** | **4 h 44 m** |

Two of these are easy to skip and both are published results. **`make q10`** is
a second complete pass over all seven replicates in both modes, costing about as
much as `make all`. **`make divergence`** is over a third of the total on its
own; it is the only step outside the Nextflow DAG, running breseq over the
contaminant reads of a finished run in its own container, and it produces
`results/contaminant_divergence/summary.tsv` — the evidence that the stock
reference is close enough to the carrier-derived contaminant to be used as one.
A reproduction that skips either is incomplete.

`make verify` goes last: it checks display items that step 4 builds. See
[`docs/operating-notes.md`](docs/operating-notes.md#make-verify-goes-last).

### Verified

On 15 August 2026 the sequence above was run end to end from a fresh `git clone`
with containers rebuilt from `docker/`, on a machine holding no prior state from
this project, **reading the seven deposited masked FASTQs** — each md5-verified
against `assets/deposited_files.tsv` before the run. Every reported value
reproduced:

- `results/summary/experiment_summary.tsv` byte-identical to the published
  tables; `per_organism.tsv` and `coverage_attribution.tsv` identical
- all twelve values of the Q10 analysis to printed precision, and its quoted
  retention ranges — 89.76–98.38% of reads and 77.20–91.16% of bases against the
  stated 89.8–98.4% and 77.2–91.2%
- all seven rows of the contaminant-divergence table, including the per-replicate
  SNP counts, every variant call identical
- all fourteen prior-study rows in `assets/comparison/prior_studies.tsv`
  re-derived from the deposited Basapathi Raghavendra reads
- `make verify` 7/7 and `make check` clean

The two pipeline invocations took 1.54 h and 1.45 h of wall time on the
reference platform (54 and 52 tasks; see `results/pipeline_info/`), the
downstream targets minutes more.

This is stated because "the code is available" and "the code reproduces the
paper" are different claims, and only the second is worth much.

**What a reader can and cannot yet check.** The run above read the deposited
files from local disk. The deposit is public — runs `SRR40180147`–`SRR40180153`,
released 15 Aug 2026 — but ENA has not yet mirrored it, so `--fetch_from_sra`
cannot resolve them today. Downloading the seven FASTQs by hand from NCBI and
supplying them through the samplesheet reproduces everything above; the one-command
path becomes available when ENA catches up. When it does, expect the archive's
regenerated copy rather than the submitted bytes — see
[Data availability](#data-availability) for what that does and does not change.

Each target writes to its own subdirectory of `results/`, so a clone has exactly
one output tree and no run can overwrite another's tables. Which target produces
what, and why `make all` rather than `make s1` plus `make s2`, are in
[`docs/operating-notes.md`](docs/operating-notes.md#where-each-target-writes).

`./run.sh` is a thin wrapper around `nextflow run` that passes every argument
through; it exists because Nextflow cannot run from a project path containing
spaces. If your clone path has none, `nextflow run .` works identically. See
[`docs/operating-notes.md`](docs/operating-notes.md#why-runsh-instead-of-nextflow-run).

---

## Design principles

These mirror the conventions used across the lab's analysis repositories.

1. **Docker-first.** The host runs Nextflow and Docker, nothing else. There is no
   conda environment, no `pip install`, no host-installed `minimap2`. Every
   process carries a `container` directive pinned to an immutable tag
   (`conf/base.config`), so a fresh clone reproduces the published numbers from
   the same inputs. Images are built from pinned source (see `docker/`) rather
   than pulled from a floating upstream tag.
2. **Every artifact carries provenance.** `BUILD_REFERENCE` emits
   `reference_provenance.json` recording each genome's accession, size, contig
   count and a SHA-256 of its sequence. Nextflow writes a timestamped timeline,
   report, trace and DAG to `results/pipeline_info/` on every run. Reference
   accessions are read out of the same reference-set TSV that defines the
   analysis, so the download and the analysis cannot drift apart.
3. **The repository is CODE ONLY.** `data/`, `refs/`, `results/` and `work/` are
   gitignored. Reads are supplied by the user, reference genomes are fetched by
   the pipeline, and results are regenerated. The only things tracked are the
   workflow, the scripts, the container definitions, the samplesheets, the
   reference-set definitions, and this documentation.
4. **Self-contained, relative paths.** Samplesheets reference `data/...` and
   `assets/references/...` relative to the repository root; relative paths are
   resolved against `projectDir` at parse time. `params.outdir` and
   `params.refdir` default to `${projectDir}/results` and `${projectDir}/refs`.
   Nothing points outside the clone.

---

## Repository layout

| Path | Contents |
|---|---|
| `main.nf` | The workflow: fetch genomes, build the combined reference, map, assign, compute metrics, profile coverage |
| `nextflow.config` | Manifest, all `params` defaults, `docker`/`singularity`/`test` profiles, provenance reporting |
| `conf/base.config` | Container pinning, per-label CPU/memory/time, `check_max` clamping |
| `run.sh` | Launcher that works around Nextflow's inability to handle spaces in the project path |
| `Makefile` | Runs: `all` (the published analysis), `s1`, `s2`, `raghavendra`, `q10`, `test`. Display items: `comparison`, `coverage`, `modedelta`, `poolcov`, `estcontrol`. Checks: `verify`, `check`, `measurements`. Records: `versions`, `runmeta`, `seqsummary`, `attribution`, `assigneddepth`. Plus `images`, `demo-data`, `clean`, `help` |
| `bin/build_reference_set.py` | Concatenates fetched genomes into one FASTA + contig→organism→role map + genome sizes + provenance |
| `bin/assign_reads.py` | Competitive per-read assignment from a qname-grouped BAM; emits counts, per-read calls, read lengths |
| `bin/compute_metrics.py` | Enrichment and per-femtogram metrics; enforces the read-accounting reconciliation |
| `bin/sequencing_summary.py` | `make seqsummary` — per-replicate reads, bases, median read length and median ONT qscore, before and after depleting carrier and contaminant. Quality comes from the FASTQ `qs` tag, joined positionally to the assignments with a per-read identity assertion |
| `bin/coverage_attribution.py` | `make attribution` — per-replicate alignment depth beside the depth attributable to awarded reads. The 1x interpretability threshold is applied to the latter: a pair whose awarded reads cover 0.83x has not been sequenced deeply enough to characterise, whatever its 65x alignment depth reads. Also excludes the smoke test, which the coverage globs otherwise sweep in |
| `bin/assigned_depth.sh` | `make assigneddepth` — per-base depth of the reads assignment *awarded* to each organism, recovered from a finished run without re-mapping the full FASTQ. `COVERAGE_PROFILE` measures every primary alignment instead, which for a member sharing sequence with an abundant relative is largely that relative's reads |
| `bin/pool_coverage.py` | `make poolcov` — sums per-base depth across the replicates of an experiment, so every community member can be assessed rather than only those deep enough in a single library. Reports alignment depth **and** the depth attributable to reads assignment awarded, because for a member sharing sequence with an abundant relative these differ by orders of magnitude |
| `bin/plot_pooled_coverage.py` | The pooled-coverage figure — uniformity across every member, pooled; marks the members whose depth is not attributable |
| `bin/coverage_dropouts.py` | Locates low-coverage regions in a depth profile and annotates them against a GFF3; answers *where* and *what*, where the coverage summary answers only *how uneven* |
| `assets/references/lowinput_s1.tsv` | Reference set A: D6311 community (10 organisms) + lambda carrier + K-12 contaminant |
| `assets/references/lowinput_s2.tsv` | Reference set B: D6321 spike-in (3 organisms) + lambda carrier + K-12 contaminant |
| `assets/samplesheets/*.csv` | Per-dataset sample definitions: keys and FASTQ paths |
| `assets/measurements.tsv` | The experimental quantities the headline numbers divide by — DNA masses, each with the basis it was obtained on, and `include_in_headline`. Read directly by the workflow; authoritative for every sample it names. Kept separate from the samplesheets so a measurement is never confused with a local file path |
| `bin/check_measurements.py` | `make measurements` — asserts the two files agree and that nothing nominal or unmeasured feeds a headline statistic |
| `bin/software_versions.sh` | `make versions` — records what the built images actually contain, by asking them rather than by transcribing the Dockerfiles |
| `bin/run_metadata.py` | `make runmeta` — per-run acquisition id, basecalling model and acquisition window, derived from the reads. The window is a minimum/maximum over *every* read: a dorado FASTQ is not sorted by time, so reading the first header gives the wrong date (it did, by two days, in a manuscript table). Slow by necessity — a full pass over every FASTQ |
| `docker/tools/` | minimap2 2.28, htslib/samtools 1.21, seqkit 2.8.2, NCBI `datasets`, built natively for amd64 and arm64 |
| `docker/analysis/` | Python 3.12 + pinned `requirements.txt` (pysam, pandas, numpy, matplotlib, scipy, openpyxl) |
| `docs/benchmarks.md` | Native arm64 vs emulated amd64 timing, projected run cost, `breseq` notes |
| `docs/ecoli-partitioning.md` | Why the community's *E. coli* must not be subtracted with the contaminant, and how competitive assignment handles it |
| `docs/TODO.md` | Honest list of known gaps |
| `docs/paper-crosswalk.md` | Which output backs which figure and table in the paper |
| `docs/history-rewrite.md` | The 15 Aug 2026 history rewrite: what was removed and why |
| `bin/comparison/` | The prior-work comparison: figure generator, the loader, the Kraken2/Raghavendra fetch and reanalysis scripts, and `seed_this_study.py` |
| `assets/comparison/` | Its inputs: `prior_studies.tsv` (one row per prior-study sample per classifier variant, each citing its published source), `this_study.tsv` (a committed snapshot of pipeline output), `kraken2_db.manifest.tsv` (the pinned database) |
| `docs/comparison.md` | How that comparison is built, and the four data defects it repairs |
| `data/`, `refs/`, `results/`, `work/` | **Gitignored.** Input reads, fetched genomes, outputs, Nextflow scratch |

---

## The method

### Competitive assignment vs sequential subtraction

The classic way to handle a carrier is **sequential subtraction**: map everything
to the carrier and delete what sticks; map the survivors to the contaminant and
delete what sticks; map whatever is left to the community. This is cheap and
intuitive, and it is wrong whenever a community member shares sequence with the
thing being subtracted. In `lowinput_s1` it does exactly that: the community
contains *E. coli* B-1109 and the contaminant is *E. coli* K-12 MG1655, so the
subtraction step deletes the community's own *E. coli* along with the carrier
background. See `docs/ecoli-partitioning.md`.

**Competitive assignment** maps each read once, with `minimap2 -ax map-ont -N 10
--secondary=yes`, against a *single combined index* holding the carrier, the
contaminant, and every community member. For each read the best alignment score
(`AS`) achieved against each organism is compared, and the read is awarded to the
winner only if it beats the runner-up organism by a margin of
`max(10, 0.01 * winning AS)`. Reads whose top two organisms are effectively tied
are reported as an explicit `ambiguous:A|B` class naming both, rather than being
silently awarded to one of them. A read that only one organism aligned at all is
assigned to that organism: ambiguity is a statement about two organisms being
indistinguishable, not about a single alignment being weak.

### What the modes mean

| `--mode` | Meaning | Status |
|---|---|---|
| `competitive` (default) | One minimap2 pass against the combined index; best organism wins by AS margin, ties reported as ambiguous | **Implemented** |
| `sequential` | The classic subtraction chain: remove carrier, then remove contaminant, then map survivors to the community | **Implemented** |
| `both` | Run both and emit a per-organism delta table quantifying what subtraction costs (`results/summary/mode_delta.*`) | **Implemented** |

`bin/contaminant_divergence.sh` measures how far the carrier-derived
contaminant actually is from the stock reference, using a finished run's
assignments rather than re-mapping. Across all seven replicates 86-92% of
contaminant reads map to stock MG1655, with no structural variants and at
worst ~100 SNPs in 4.64 Mb, so the stock reference is adequate and the
consensus below is not needed for these data. See `docs/TODO.md` §6.

`--breseq_consensus` (sequential and both only) subtracts against a
reference-guided consensus of the contaminant actually present in the carrier
prep, built with breseq, instead of against the stock reference. That is what
the original `lowinput_s1` analysis did, so the flag exists to reproduce it
faithfully. It needs real depth: below `--breseq_min_depth` (default 10×) breseq
predicts missing coverage across the whole reference and returns a deleted
genome rather than a consensus, so the bundled test profile — ~0.3× contaminant
depth — cannot exercise it. The accounting is covered by `make check`. All seven
replicates have since been run through it, and the consensus route did not
change the conclusion; see `docs/TODO.md` §6.

### The read-accounting guarantee

Every read in the input FASTQ lands in **exactly one** of five buckets:

| Bucket | Meaning |
|---|---|
| `sample` | Assigned to a community member (`role=sample` in the reference set) |
| `carrier` | Assigned to lambda |
| `contaminant` | Assigned to *E. coli* K-12 (carried over from lambda production) |
| `ambiguous` | Top two organisms within the margin; reported as `ambiguous:A\|B`, naming both |
| `unassigned` | No alignment survived the MAPQ floor, or the read was unmapped |

`assign_reads.py` streams a qname-grouped BAM and emits one call per read, so
the per-organism counts sum to the FASTQ read count by construction.
`compute_metrics.py` then re-checks it:

```
if bucket_reads != total_reads:
    sys.exit("error: read accounting does not reconcile: ...")
```

If the buckets do not sum to the total, **the run fails**. This is deliberate: the
denominator of `output_sample_fraction` is the total read/base count, so a leak
anywhere in the assignment silently inflates the headline enrichment. Failing
loudly is preferred to publishing an unverifiable number.

---

## Metric definitions

These are the definitions implemented in `bin/compute_metrics.py`; they are
stated explicitly because the comparison to prior work hinges on them.

**`input_sample_fraction`**
`library_dna_ng / (library_dna_ng + carrier_dna_ng)` — the fraction of the mass
entering library prep that is sample rather than carrier.

**`output_sample_fraction`**
sample bases (or reads) / all basecalled bases (or reads), where the denominator
is *every* read in the FASTQ including carrier, contaminant, ambiguous and
unmapped. Counts come from the assignment step, which accounts for every read
exactly once, so this denominator is exact.

**`enrichment`**
`output_sample_fraction / input_sample_fraction`: how much the sample's share of
the output exceeds its share of the input. Reported both base-weighted and
read-weighted. An earlier preliminary analysis of `lowinput_s1` reported this as
">100x"; correcting the carrier mass and removing carrier-derived *E. coli* from
the community lowers it to **53.6x** (manuscript, Supplementary Section S2).
That figure moved once more, from 70.7x, when the attribution floor
(`--min_aln_frac`) stopped an organism claiming a whole read on a sliver of
alignment; see [Metric definitions](#metric-definitions) above and
`bin/attribution_threshold.py` for the distribution behind the threshold.

**`reads_per_fg` / `bases_per_fg`**
reads (or bases) assigned to an organism, divided by the femtograms of that
organism's DNA that entered library prep, i.e.
`library_dna_ng * 1e6 * theoretical_dna_fraction`. For the `All organisms` row
the denominator is `library_dna_ng * 1e6`. This matches the convention used for
the prior-study comparison, where the denominator is DNA into library prep
("post-extraction").

**`bases`** — a deliberate choice.
`bases` is the **full read length** of every read assigned to the organism, *not*
the aligned span. The prior studies being compared against count Kraken2-classified
read lengths, so counting aligned bases here would understate this study relative
to them. Aligned bases are reported alongside as `aligned_bases` for internal use
(coverage, identity). Read length is taken from `infer_read_length()`, which
counts hard clips and is therefore the only value consistent across primary,
secondary and supplementary records of the same read.

Masses come from `assets/measurements.tsv`, not from the samplesheets. Each is
paired with a `basis` column recording how it was obtained, and
`include_in_headline` is an explicit field there. A row with no usable mass
emits blank per-fg fields rather than a fabricated number.

The separation is deliberate. A samplesheet mixes local file paths — meaningless
in someone else's clone — with experimental facts that belong to the paper.
Keeping a second copy of a mass beside a path is how a stale carrier value once
survived next to a Methods section that contradicted it. The pipeline now aborts
if a samplesheet carries a mass that disagrees with `measurements.tsv`, and
`make measurements` fails if an unmeasured value feeds a headline statistic or
an exclusion has no stated reason.

---

## Outputs

Everything lands under `params.outdir` (default `results/`).

| Path | Produced by | Contents |
|---|---|---|
| `results/references/<set>/combined.fasta` | `BUILD_REFERENCE` | The single index that competitive assignment maps against |
| `results/references/<set>/contig_map.tsv` | `BUILD_REFERENCE` | `contig` → `organism` → `role` → `accession`; contig names are left untouched so alignments stay traceable to the NCBI record |
| `results/references/<set>/genome_sizes.tsv` | `BUILD_REFERENCE` | Per-organism genome size, contig count, theoretical DNA fraction |
| `results/references/<set>/reference_provenance.json` | `BUILD_REFERENCE` | Accession, size, contig count and SHA-256 of the sequence of every genome used |
| `results/<sample_id>/alignments/<sample_id>.qname.bam` | `MAP_COMPETITIVE` | qname-grouped BAM (only when `--keep_bams`, the default) |
| `results/<sample_id>/competitive/<sample_id>.counts.tsv` | `ASSIGN_READS` | Per-class `reads`, `read_bases`, `aligned_bases` |
| `results/<sample_id>/competitive/<sample_id>.assignments.tsv.gz` | `ASSIGN_READS` | Per-read: organism, call, role, best AS, runner-up AS, margin, read length, aligned bases |
| `results/<sample_id>/competitive/<sample_id>.readlengths.tsv.gz` | `ASSIGN_READS` | Per-read length by assignment class, for the adaptive-sampling ejection-signature analysis |
| `results/<sample_id>/competitive/<sample_id>.metrics.tsv` | `COMPUTE_METRICS` | Per-organism table: reads, bases, aligned bases, theoretical vs measured fraction, `dna_fg`, `reads_per_fg`, `bases_per_fg`, plus an `All organisms` row |
| `results/<sample_id>/competitive/<sample_id>.summary.json` | `COMPUTE_METRICS` | Totals by role, `input_sample_fraction`, `output_sample_fraction`, `enrichment`, per-fg headline numbers |
| `results/<sample_id>/coverage/<sample_id>.depth.tsv.gz` | `COVERAGE_PROFILE` | Per-base depth over community (`role=sample`) contigs only; carrier depth is uninformative and would dominate the file |
| `results/<sample_id>/human/<sample_id>.masked.fastq.gz` | `MASK_HUMAN` | The deposited reads: same records, same lengths, human-attributable intervals replaced with `N` (only with `--mask_human`) |
| `results/<sample_id>/human/<sample_id>.human_mask.tsv.gz` | `MASK_HUMAN` | Per-read audit trail for every HRRT-flagged read: call, organism, masked intervals, reason |
| `results/<sample_id>/human/<sample_id>.human_stats.json` | `MASK_HUMAN` | Counts, percentages, the rule applied and its provenance |
| `results/pipeline_info/{timeline,report,trace,dag}_<timestamp>.*` | Nextflow | Provenance for every run |

Fetched genomes are cached in `params.refdir` (default `refs/`) via `storeDir`, so
re-runs do not re-download from NCBI.

---

## Human read masking

Human sequence has to be removed before reads from a human-handled sample can
be deposited publicly. Run it with:

```bash
make all NF_ARGS="--mask_human true"
make masksummary          # per-replicate statistics table
```

It is **off by default**, so any run without it is bit-identical to one from
before the feature existed.

**The screen is not the decision.** NCBI's Human Read Removal Tool (HRRT) flags
candidates, and HRRT alone is far too blunt for a microbial community. Its
false-positive rate tracks evolutionary distance from human: across the seven
replicates it flags **18.6%** of reads assigned to *S. cerevisiae* (663 of
3,567) and **14 of the 20** assigned to *C. neoformans*, against **0.20%** of
those assigned to any bacterial community member (1,441 of 706,347) — a
hundred-fold difference, because a conserved region is a conserved region
whichever genome it sits in. So a flagged read is released intact only when this
pipeline's own competitive assignment **positively attributes** it to a
community organism. The inverse rule — "mask unless it failed to look human" —
reads almost identically and is backwards for a privacy filter; do not
refactor toward it.

Three further properties are worth knowing before reading the outputs:

- **Nothing is deleted.** Masking replaces bases with `N` in place, so read
  count, read order and every read length are preserved exactly. The
  per-record length invariant is asserted at run time and a violation is fatal.
- **A rescued read can still be partly masked.** The kept interval is the part
  some organism accounts for; anything human-exclusive goes, and a read with
  at least `--chimera_min_bp` (default 150) of human-exclusive sequence is
  masked regardless of its assignment. Whole-read best-hit assignment cannot
  catch those on its own — the read wins its comparison on the microbial half.
- **The flag list is the only thing taken from HRRT.** Its own masked output is
  discarded; the intervals come from a separate `map-ont` alignment of the
  flagged reads against T2T-CHM13v2.0 (`--human_accession`, cached in
  `refs/human/`).

`make check` covers the rule with 18 logic assertions and 12 more against a
committed 1,520-read fixture (`assets/testdata/`) that includes GIAB HG002
reads, conserved-region reads a real HRRT run flagged, and synthetic chimeras
spliced at recorded offsets. Neither suite needs the network or a mapping pass.

A third check looks the other way. The fixture is built from reads **before**
masking, so it can hold sequence the deposit withholds — and nothing that tests
the masker can see that, because the masker is asserted against the fixture's
own records. `tests/fixture_deposit_agreement.py` compares every study-derived
fixture read against its released form: base for base when the deposited reads
are present (`make check DEPOSIT=/path/to/lowinput_s1_r1.fastq.gz`), and against
a committed list of the ids the deposit masks when they are not.

**What may and may not be claimed.** Reads confidently assignable to human were
masked. That is not the same as "no human sequence remains" — recall is roughly
91.5% on alignable reads — and identity to CHM13 is a measure of similarity,
not of read quality.

The `--mask_human` run needs the amd64 `scrubber` image (HRRT has no arm64
build) and about 9 GB of memory for `MAP_HUMAN`, whose footprint is set by the
7.55 GB human index rather than by read volume.

---

## Containers

Pinned in `conf/base.config`. Nothing runs on host-installed software.

The `:0.1.0` on the image names is the **image** version and moves only when
an image's contents change; it is deliberately independent of the pipeline
version in `nextflow.config` and `CITATION.cff`. A v1.0.0 pipeline running
`tools:0.1.0` means the tools image has not needed a rebuild, not that
something is out of step.

| Label | Image | Source | Used by |
|---|---|---|---|
| `tools` | `low-input-nanopore/tools:0.1.0` | `docker/tools/Dockerfile` — built locally by `make images`; minimap2 2.28, htslib 1.21, samtools 1.21, seqkit 2.8.2, NCBI `datasets` v2, on `debian:bookworm-20241111-slim` | `FETCH_GENOMES`, `FETCH_READS`, `MAP_COMPETITIVE`, `COVERAGE_PROFILE` |
| `analysis` | `low-input-nanopore/analysis:0.1.0` | `docker/analysis/Dockerfile` — built locally by `make images`; `python:3.12-slim-bookworm` + pinned `requirements.txt` (numpy 2.1.3, pandas 2.2.3, matplotlib 3.9.2, openpyxl 3.1.5, pysam 0.22.1, scipy 1.14.1) | `BUILD_REFERENCE`, `ASSIGN_READS`, `COMPUTE_METRICS` |
| `breseq` | `quay.io/biocontainers/breseq:0.40.1--h3be2455_0` | Upstream biocontainer, forced to `--platform linux/amd64`, run as root with `-e HOME=/tmp` | The optional sequential/`--breseq_consensus` path |
| `scrubber` | `low-input-nanopore/scrubber:0.1.0` | `docker/scrubber/Dockerfile` — built locally by `make images`; NCBI `sra-human-scrubber` pinned by digest, plus `procps` | `HRRT_SCREEN`, on the optional `--mask_human` path |

The two local images are built from pinned source rather than pulled from
bioconda because bioconda publishes **linux/amd64 only**, while this study's
reference platform is Apple Silicon (arm64). `breseq` deliberately stays on the
emulated amd64 image: the arm64 conda build of its `bowtie2` dependency has a
known SIMD/SIGILL risk, and that step runs on a small read subset where
emulation cost is immaterial. It must run as root (`breseq` fails under a mapped
non-root user), which is why that label overrides the global
`docker.runOptions` user mapping.

For the measured native-arm64-vs-emulated-amd64 comparison (minimap2 was 1.41x
faster native — a ~40% saving, not the order of magnitude emulation is often
assumed to cost) and the projected full-run cost, see
**[`docs/benchmarks.md`](docs/benchmarks.md)**.

A `singularity` profile is also defined in `nextflow.config`.

---

## Data availability

The seven replicates are deposited in the NCBI Sequence Read Archive under
**BioProject [PRJNA1513130](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1513130)**,
one BioSample per replicate (`SAMN62407365`–`SAMN62407371`, listed per row in
the samplesheets). The deposited FASTQs are the **human-masked** files described
above, and every published value was computed from them — see
[Human read masking](#human-read-masking).

**This repository's git history was rewritten once, on 15 Aug 2026 and before it
was ever public**, to remove six reads whose deposited form is masked — the test
fixture had been holding them unmasked. What was removed, how it happened, and
how the result was verified are recorded in
[`docs/history-rewrite.md`](docs/history-rewrite.md). No published result is
affected; the rewrite touched only `assets/testdata/`, which no part of the
pipeline reads.

The deposit was **released on 15 August 2026**. The seven runs are
`SRR40180147`–`SRR40180153`, carried per row in the samplesheets'
`sra_accession` column; each run's read and base counts were checked against
[`assets/deposited_files.tsv`](assets/deposited_files.tsv) and agree exactly,
totalling 60,416,747 reads and 54,642,391,422 bases.

> **ENA has not yet mirrored the project, so `--fetch_from_sra` cannot resolve
> these accessions yet.** The records are public at NCBI; ENA, which this
> repository fetches from, receives INSDC submissions on a delay of a few days.
> Until it does, run without the flag and supply the FASTQs through the
> samplesheet's `fastq` column. `bin/fetch_ena_reads.sh` says which of the two
> situations it is in rather than reporting both as "not released".

### Fetching the reads by accession

```bash
nextflow run . -profile docker \
  --samplesheet assets/samplesheets/all.csv \
  --fetch_from_sra
```

The samplesheet's local `fastq` column is ignored; every sample is fetched
instead, into the gitignored `sra/` cache (`--sradir`). A repository clone plus
a network is then the whole input.

Resolution goes through ENA's portal API rather than the SRA toolkit, for two
reasons: the toolkit publishes no `linux/arm64` build, and this study's
reference platform is Apple Silicon; and ENA *can* serve the **submitted**
file — the exact bytes uploaded — alongside the archive's regenerated copy.
Only the submitted copy can be checked against
[`assets/deposited_files.tsv`](assets/deposited_files.tsv), which records the
md5, size, read count and base count of each deposited file. `bin/fetch_ena_reads.sh`
prefers it, verifies both the archive's own checksum and ours, and says which
copy it used. If only the regenerated copy exists it falls back with a warning
and skips our checksum, because SRA rebuilds `fastq_ftp` from its own archive:
the same sequence, different bytes, and possibly different read names.

**Expect that fallback here.** This submission was made to NCBI, and NCBI
serves only the `.sra` archive format for these runs — no original submitted
file — so the regenerated copy is likely the only one either archive will ever
offer. A public download therefore reproduces the *sequence* but not the
*bytes*, and the md5s in `deposited_files.tsv` record what was uploaded rather
than what can be downloaded. The reads themselves are unaffected: read and base
counts were verified against NCBI's own records per run. One consequence worth
anticipating is that SRA does not preserve ONT header tags, so `qs:f:` is
absent from a fetched FASTQ and `bin/sequencing_summary.py` will report an
empty `median_qscore`.

The samplesheets carry both a BioSample accession and, since the 15 Aug 2026
release, the run accession itself. `sra_accession` wins when present, so
resolution is now unambiguous rather than relying on a BioSample mapping to
exactly one run.

`tests/sra_fetch.sh` asserts every branch of that logic against saved portal
responses, with no network: which copy is chosen, that a checksum mismatch
deletes the file rather than caching a corrupt one, and that the two accession
tables agree with each other.

Reference genomes are *not* redistributed here either: they are fetched from NCBI
by accession at run time, driven by the same reference-set TSV that defines the
analysis, and cached in the gitignored `refs/`.

---

## Citation

If you use this pipeline, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff) (GitHub renders a formatted citation from it via
"Cite this repository").

Authors: Jordan McKaig (School of Earth and Atmospheric Sciences; Daniel
Guggenheim School of Aerospace Engineering, Georgia Institute of Technology),
Jenna Sailor (Wallace H. Coulter Department of Biomedical Engineering, Georgia
Institute of Technology), Christopher E. Carr (Daniel Guggenheim School of
Aerospace Engineering; School of Earth and Atmospheric Sciences, Georgia
Institute of Technology).

Corresponding author: Christopher E. Carr.

## License

MIT — see [`LICENSE`](LICENSE). Copyright (c) 2026 Christopher E. Carr / Georgia
Institute of Technology.

Nothing third-party is vendored here: every tool is pulled at build or run time
by pinned digest, and no reference genome is redistributed.
[`THIRD_PARTY.md`](THIRD_PARTY.md) lists what is used and under what terms,
explains why breseq being GPL-2.0 does not reach this MIT code, and carries the
required notices for the two datasets that **are** redistributed — the GIAB
HG002 reads in `assets/testdata/` (CC0, with a statement of how they were
modified) and the prior-study values in `assets/comparison/`.
