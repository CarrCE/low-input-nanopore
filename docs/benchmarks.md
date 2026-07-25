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

## Projected full-run cost

Scaling the native figure to the full dataset (7 replicates, ~119 GB of FASTQ)
gives roughly **30-40 minutes of minimap2 wall-clock at 8 threads** against a
small reference, and proportionally more against the 10-organism `lowinput_s1`
reference set. Mapping is not the bottleneck for this study; per-read
assignment and coverage summarisation dominate. Actual timings are recorded in
`results/pipeline_info/` by the Nextflow trace.

## Note on `breseq`

`breseq` is only used in the optional `--mode sequential --breseq_consensus`
path. It stays on an amd64 image: the arm64 conda build of its `bowtie2`
dependency has a known SIMD/SIGILL risk, and the step runs on a small subset of
reads where the emulation cost is irrelevant. That container must run as root
(`breseq` fails as a non-root user), which is why the `withLabel: breseq`
process overrides the global `docker.runOptions` user mapping.
