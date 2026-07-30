# Benchmarks and platform notes

## Container architecture: native arm64 vs emulated amd64

All the usual bioinformatics containers (`quay.io/biocontainers/*`,
`ensemblorg/datasets-cli`) publish **linux/amd64 only**. This study's reference
platform is an Apple M3 Max (arm64), so those images run under Docker Desktop's
Rosetta-backed emulation.

To decide whether to build native images, both paths were timed on the same
input: a 400,000,001-byte slice of `lowinput_s2_r1.fastq` (260,835 reads),
mapped with `minimap2 -ax map-ont -t 8 --secondary=no` against a combined
lambda + *E. coli* K-12 reference (4.75 Mb).

| Image | Arch | minimap2 wall-clock | CPU time | Peak RSS |
|---|---|---|---|---|
| `low-input-nanopore/tools:0.1.0` | native arm64 | **6.64 s** | 47.4 s | 2.53 GB |
| `quay.io/biocontainers/minimap2:2.28--he4a0461_2` | amd64, emulated | 9.35 s | 68.5 s | 2.54 GB |

Native is **1.41x faster** — a ~40% saving, not the order of magnitude that
emulation is often assumed to cost. Rosetta handles minimap2's SIMD workload
well.

The native image is used anyway because:

1. 40% is material across the full ~45 Gbp dataset.
2. Building from pinned source produces a reproducible image on **both**
   architectures, so the pipeline does not depend on Rosetta being present
   (it is unavailable on arm64 Linux, e.g. Graviton or Ampere CI runners).
3. It removes a dependency on upstream biocontainer tags remaining published.

Both are identical in output: same minimap2 `2.28-r1209`, same CMD, same
260,835 sequences mapped.

## Measured full-run cost

An earlier version of this document projected "30-40 minutes of minimap2
wall-clock" for the whole dataset by scaling the micro-benchmark above. That
projection was wrong, and it is worth recording why: the benchmark used
`--secondary=no` against a 4.75 Mb reference, whereas the pipeline runs
`-N 10 --secondary=yes` against a 13.8 Mb reference, because the assignment
margin needs the runner-up organism's alignment score.

Measured on the real `lowinput_s2` run (4 replicates, ~73 GB of FASTQ,
M3 Max, 14 CPUs available to Docker; from `results/pipeline_info/trace_*.txt`):

| Process | Per replicate (range) | Peak RSS | Notes |
|---|---|---|---|
| `MAP_COMPETITIVE` | 7 m 21 s – 55 m 15 s | 8.5–9.3 GB | ~800–940% CPU |
| `ASSIGN_READS` | 2 m 14 s – 19 m | ~32 MB | streams the BAM; memory-flat |
| `COVERAGE_PROFILE` | 43 m 12 s – 58 m 5 s | 6.2 GB | **was the bottleneck — see below** |
| `COVERAGE_SUMMARY` | ~3 s | 0.7–0.9 GB | after vectorisation |
| `COMPUTE_METRICS` | < 100 ms | ~12 MB | |
| `AGGREGATE` | 35.7 s | 4.5 GB | pools ~35 M read-length rows |

Total for `lowinput_s2`: **3 h 28 m wall-clock, 44.7 CPU-hours** — confirmed from
`trace_2026-07-25_01-54-14.txt` (17 fresh tasks, 3 h 27 m 47 s).

> **That total is pre-fix and is no longer the cost of running this pipeline.**
> It is retained because it is what the two fixes below were measured against.
> The table above it is the same run, so it is pre-fix too: `COVERAGE_PROFILE` at
> 43–58 min per replicate is the defect, not the current behaviour. Note also
> its scope — 4 replicates in **one** mode, not the 7 replicates in both modes
> that `make all` runs.

### Post-fix cost

From `trace_2026-07-25_09-17-25.txt`, the 7-replicate both-modes run (14 cores to
Docker, M3 Max):

| Process | n | Per task | Notes |
|---|---:|---|---|
| `MAP_COMPETITIVE` | 7 | 6 m 04 s – 8 m 30 s | 800–1080% CPU, 8.4–10 GB RSS |
| `ASSIGN_READS` | 14 | 2 m 09 s – 3 m 13 s | single-threaded, ~30 MB |
| `COVERAGE_PROFILE` | 7 | 1 m 32 s – 2 m 08 s | was 43–58 min |
| `COVERAGE_SUMMARY` | 7 | 2 – 13 s | was ~18 min |
| `AGGREGATE` | 3 | ~1 m 20 s | one per mode plus the pooled table |

`MAP_COMPETITIVE` is `process_high`, so `cpus = max_cpus` and **one mapping task
runs at a time**. That makes it the serialised spine and a hard floor on wall
clock: 48 m 43 s of mapping for seven replicates, measured. Everything
downstream is either single-threaded or seconds long and overlaps with the next
replicate's mapping.

**Caveat on the ~1–1.5 h figure quoted in the README.** No trace in
`results/pipeline_info/` is a complete from-scratch `make all`: the 7-replicate
traces are 51-of-52 cached, and the one with 28 freshly-executed tasks had
`AGGREGATE` fail. So the estimate is a reconstruction from the measured per-task
times plus the serialisation structure above, not an end-to-end measurement.
The closest real datapoint is 46 m 45 s of wall clock for 4 fresh mappings with
the rest cached, which is consistent with it. Treat 1–1.5 h as an estimate with
the mapping spine as its floor, and if you do run one from cold, replace this
paragraph with the trace.

Add for a full reproduction of every display item: `assigneddepth` (~30–60 min,
one FASTQ pass per replicate), `runmeta` (a full pass over ~110 GB), and
`divergence` (breseq under amd64 emulation over all seven replicates — the
slowest single step, and unmeasured here). Those bring the end-to-end total to
roughly **4–6 h**, which is the figure the README quotes.

### Two bottlenecks found and fixed

Both were pipeline defects rather than intrinsic cost, and both are worth
knowing about if you adapt this code:

1. **`COVERAGE_PROFILE` sorted the entire BAM.** Restricting `samtools view` by
   region requires a coordinate-sorted, indexed BAM; this one is qname-grouped
   straight from minimap2, so the region form failed and the fallback sorted all
   ~10 M reads — of which >95% are carrier and irrelevant to community coverage.
   Filtering on RNAME in the SAM stream before sorting means the sort only sees
   community reads. Per-organism statistics are unchanged; the emitted depth file
   also shrinks (9.1 M rows instead of 13.8 M) because non-community `@SQ` lines
   are dropped from the header.

2. **`COVERAGE_SUMMARY` parsed depth line-by-line in Python**, which took over
   18 minutes per replicate. Chunked parsing with pandas' C engine plus
   `numpy.bincount` produces byte-identical output in **~3 seconds**.

The lesson generalises: in a carrier-based protocol the carrier is the
overwhelming majority of every intermediate file, so any step that does not
filter it out early pays for it repeatedly.

Live timings for any run are always recorded in `results/pipeline_info/`.

## Note on `breseq`

`breseq` is only used in the optional `--mode sequential --breseq_consensus`
path. It stays on an amd64 image: the arm64 conda build of its `bowtie2`
dependency has a known SIMD/SIGILL risk, and the step runs on a small subset of
reads where the emulation cost is irrelevant. That container must run as root
(`breseq` fails as a non-root user), which is why the `withLabel: breseq`
process overrides the global `docker.runOptions` user mapping.
