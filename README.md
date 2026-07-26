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

---

## Quickstart

Requirements on the host: **Nextflow (>=23.10.0)** and **Docker**. Nothing else —
every tool the pipeline calls lives in a pinned container.

```bash
# 1. Build the two local images (minimap2/samtools/seqkit/datasets, and python)
make images

# 2. Smoke test: 40,000 reads through the whole pipeline, a couple of minutes
./run.sh -profile docker,test

# 3. The real runs
./run.sh -profile docker --samplesheet assets/samplesheets/lowinput_s2.csv
./run.sh -profile docker --samplesheet assets/samplesheets/lowinput_s1.csv
```

`make test`, `make s1` and `make s2` are shorthands for the same three commands.

### Why `./run.sh` instead of `nextflow run`

`./run.sh` is a thin wrapper around `nextflow run` and passes every argument
straight through. It exists for one reason: **Nextflow cannot run from a project
path that contains spaces.** Nextflow writes an `export PATH=...` line and an
inner `bash <path>` into each task wrapper without quoting them, so a project
directory such as `.../My Research Data/...` breaks the wrapper
before any command executes. That is a Nextflow limitation, not something the
pipeline can fix internally.

When the repository sits at a path containing spaces, `run.sh` creates a stable,
space-free symlink to the repo under `$TMPDIR` (named from a hash of the real
path, so repeated runs reuse it and `-resume` keeps working) and launches
Nextflow through that link. Nextflow then generates space-free `PATH` and
work-directory references, while still bind-mounting the real (spaced) location
for staged inputs — a path it does escape correctly. If `TMPDIR` itself contains
spaces, `run.sh` stops with an instruction to set a space-free `TMPDIR`.

If your clone path has **no** spaces, `run.sh` simply execs `nextflow run` in
place, and a direct

```bash
nextflow run . -profile docker --samplesheet assets/samplesheets/lowinput_s2.csv
```

works identically. Use `run.sh` anyway if you want the same command to be
portable between both situations.

Related: the workflow stages `bin/*.py` as explicit process inputs rather than
relying on Nextflow's `bin/` PATH injection, which is the other thing that breaks
under a spaced project path.

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
| `Makefile` | `images`, `test`, `s1`, `s2`, `demo-data`, `clean`, `help` |
| `bin/build_reference_set.py` | Concatenates fetched genomes into one FASTA + contig→organism→role map + genome sizes + provenance |
| `bin/assign_reads.py` | Competitive per-read assignment from a qname-grouped BAM; emits counts, per-read calls, read lengths |
| `bin/compute_metrics.py` | Enrichment and per-femtogram metrics; enforces the read-accounting reconciliation |
| `bin/sequencing_summary.py` | `make seqsummary` — per-replicate reads, bases, median read length and median ONT qscore, before and after depleting carrier and contaminant. Quality comes from the FASTQ `qs` tag, joined positionally to the assignments with a per-read identity assertion |
| `bin/coverage_dropouts.py` | Locates low-coverage regions in a depth profile and annotates them against a GFF3; answers *where* and *what*, where the coverage summary answers only *how uneven* |
| `assets/references/lowinput_s1.tsv` | Reference set A: D6311 community (10 organisms) + lambda carrier + K-12 contaminant |
| `assets/references/lowinput_s2.tsv` | Reference set B: D6321 spike-in (3 organisms) + lambda carrier + K-12 contaminant |
| `assets/samplesheets/*.csv` | Per-dataset sample definitions: keys and FASTQ paths |
| `assets/measurements.tsv` | The experimental quantities the headline numbers divide by — DNA masses, each with the basis it was obtained on. Kept separate from the samplesheets so a measurement is never confused with a local file path. **Draft; see the PENDING entries.** |
| `bin/check_measurements.py` | `make measurements` — asserts the two agree and that nothing nominal or unmeasured feeds a headline statistic |
| `docker/tools/` | minimap2 2.28, htslib/samtools 1.21, seqkit 2.8.2, NCBI `datasets`, built natively for amd64 and arm64 |
| `docker/analysis/` | Python 3.12 + pinned `requirements.txt` (pysam, pandas, numpy, matplotlib, scipy, openpyxl) |
| `docs/benchmarks.md` | Native arm64 vs emulated amd64 timing, projected run cost, `breseq` notes |
| `docs/ecoli-partitioning.md` | Why the community's *E. coli* must not be subtracted with the contaminant, and how competitive assignment handles it |
| `docs/TODO.md` | Honest list of known gaps |
| `modules/local/` | Reserved for extracted process modules |
| `comparison/` | Reserved for the prior-study reanalysis (see `docs/TODO.md`) |
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
consensus below is not needed for these data. See `docs/TODO.md` item 2.

`--breseq_consensus` (sequential and both only) subtracts against a
reference-guided consensus of the contaminant actually present in the carrier
prep, built with breseq, instead of against the stock reference. That is what
the original `lowinput_s1` analysis did, so the flag exists to reproduce it
faithfully. It needs real depth: below `--breseq_min_depth` (default 10×) breseq
predicts missing coverage across the whole reference and returns a deleted
genome rather than a consensus, so the bundled test profile — ~0.3× contaminant
depth — cannot exercise it. The accounting is covered by `make check`; **no full
replicate has been run through it yet**, see `docs/TODO.md` item 2.

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
`output_sample_fraction / input_sample_fraction`. This is the quantity reported
as ">100x enrichment": how much depletion-mode adaptive sampling raised the
sample's share of the output above its share of the input. Reported both
base-weighted and read-weighted.

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

Rows where the input mass is unquantified (blank `library_dna_ng`, e.g.
`lowinput_s2_r0`) emit blank per-fg fields rather than a fabricated number, and
carry `include_in_headline=0` in the samplesheet.

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
| `results/pipeline_info/{timeline,report,trace,dag}_<timestamp>.*` | Nextflow | Provenance for every run |

Fetched genomes are cached in `params.refdir` (default `refs/`) via `storeDir`, so
re-runs do not re-download from NCBI.

---

## Containers

Pinned in `conf/base.config`. Nothing runs on host-installed software.

| Label | Image | Source | Used by |
|---|---|---|---|
| `tools` | `low-input-nanopore/tools:0.1.0` | `docker/tools/Dockerfile` — built locally by `make images`; minimap2 2.28, htslib 1.21, samtools 1.21, seqkit 2.8.2, NCBI `datasets` v2, on `debian:bookworm-20241111-slim` | `FETCH_GENOMES`, `MAP_COMPETITIVE`, `COVERAGE_PROFILE` |
| `analysis` | `low-input-nanopore/analysis:0.1.0` | `docker/analysis/Dockerfile` — built locally by `make images`; `python:3.12-slim-bookworm` + pinned `requirements.txt` (numpy 2.1.3, pandas 2.2.3, matplotlib 3.9.2, openpyxl 3.1.5, pysam 0.22.1, scipy 1.14.1) | `BUILD_REFERENCE`, `ASSIGN_READS`, `COMPUTE_METRICS` |
| `breseq` | `quay.io/biocontainers/breseq:0.40.1--h3be2455_0` | Upstream biocontainer, forced to `--platform linux/amd64`, run as root with `-e HOME=/tmp` | The optional sequential/`--breseq_consensus` path (not yet wired) |

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

**The reads are not yet deposited.** Until they are, every samplesheet row must
carry a local `fastq` path; the pipeline fails with an explicit message if that
column is blank.

The samplesheets reserve an `sra_accession` column, and `nextflow.config` reserves
`params.fetch_from_sra`, for the point at which the FASTQs are deposited to
SRA/ENA. Setting `--fetch_from_sra` today exits with an error — the retrieval
path is not implemented. Deposition and the wiring of that parameter are tracked
as item 1 in **[`docs/TODO.md`](docs/TODO.md)**.

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
