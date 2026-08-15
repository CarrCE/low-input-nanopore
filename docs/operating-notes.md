# Operating notes

Things worth knowing before a long run. None of it is needed to understand what
the pipeline does — that is the README — but each of these has cost somebody an
hour, and a four-hour reproduction is an expensive place to discover them.

## Getting reads into `data/`: link or copy, never symlink outward

`data/readme.md` is tracked, so `data/` already exists in a fresh clone.
`ln -s <dir> data` therefore creates `data/data`, not `data`. Link the contents
instead:

```bash
SRC=/path/to/your/fastqs
ln "$SRC"/*.fastq data/                       # hard links: no copy, no extra disk
mkdir -p data/test && ln "$SRC/test_s2.fastq" data/test/
```

**Hard-link or copy — do not symlink to a directory outside the repository.**
The Nextflow processes resolve symlinks and stage the real file, so `make all`
works either way. But `assigneddepth`, `poolcov`, `seqsummary`, `runmeta` and
`divergence` run `docker run -v "$REPO":/repo`, and inside that container a
symlink pointing outside the repo is a dangling link. Those five steps fail with
`no such file or directory` *after* the pipeline has already succeeded — the
worst time to find out.

Hard links require the reads to be on the same filesystem as the clone. If they
are not, copy.

## Run the experiments together, not one after the other

Use `make all`. `AGGREGATE` sees only the samples of the run that invokes it,
and `make s1` and `make s2` write to the same `results/summary/`, so running one
after the other leaves the study-level tables and the abundance, read-length and
mode-comparison figures describing whichever ran last. `make s1` and `make s2`
are for working on one dataset at a time, not for assembling a full result.

`--mode both` is required for the mode-comparison figure, which contrasts the
two assignment rules.

The smoke test is safe from this: the test profile publishes to `results/smoke/`
rather than `results/`. That matters because `make check` requires `make test`,
which would otherwise put a smoke test one command away from replacing every
published table and figure with 40,000 reads' worth of output.

## `make verify` goes last

It asserts that every display item satisfies `docs/display-items.md`, and the
mode-comparison, coverage and per-femtogram figures are built after the
pipeline. Run straight after the pipeline it reports the figures it has not yet
been given a chance to build.

## Where each target writes

Runs that must not overwrite the study tables publish into their own
subdirectory of `results/`, so a clone has exactly one output tree.

| target | writes to | is it a published result? |
|---|---|---|
| `make all` | `results/summary/` | Yes — the study-level tables and figures |
| `make q10` | `results/q10/` | **Yes.** The quality-matched rerun at Q10. Reproducing the paper means running this. |
| `make raghavendra` | `results/raghavendra/` | Indirectly — it is the *provenance* of the prior-study rows committed to `assets/comparison/prior_studies.tsv`. The comparison figure is drawn from that committed table, so it reproduces without this; running it re-derives those numbers from the reads. |
| `make test` | `results/smoke/` | No — 40,000 reads |

`make q10` is a second full pass over all seven replicates in both modes and
costs about as much as `make all`. `make raghavendra` is cheap: the deposited
chunks under `data/raghavendra_2023/` total ~1.3 MB.

## Why `./run.sh` instead of `nextflow run`

`./run.sh` is a thin wrapper around `nextflow run` that passes every argument
straight through. It exists for one reason: **Nextflow cannot run from a project
path that contains spaces.** It writes an `export PATH=...` line and an inner
`bash <path>` into each task wrapper without quoting them, so a project
directory such as `.../My Research Data/...` breaks the wrapper before any
command executes. That is a Nextflow limitation, not something this pipeline can
fix internally.

When the repository sits at a path containing spaces, `run.sh` creates a stable,
space-free symlink to the repo under `$TMPDIR` — named from a hash of the real
path, so repeated runs reuse it and `-resume` keeps working — and launches
Nextflow through that link. Nextflow then generates space-free `PATH` and
work-directory references while still bind-mounting the real location for staged
inputs, a path it does escape correctly. If `TMPDIR` itself contains spaces,
`run.sh` stops and asks for a space-free one.

If your clone path has no spaces, `run.sh` simply execs `nextflow run` in place,
and

```bash
nextflow run . -profile docker --samplesheet assets/samplesheets/lowinput_s2.csv
```

works identically. Use `run.sh` anyway if you want one command that is portable
between both situations.

Related: the workflow stages `bin/*.py` as explicit process inputs rather than
relying on Nextflow's `bin/` PATH injection, which is the other thing a spaced
project path breaks.

## The work directory is deliberately outside the repository

`run.sh` puts Nextflow's `work/` under `$TMPDIR`. It reaches tens of GB of
intermediate BAMs — over 300 GB for a full run including `q10` — and a clone in
a cloud-synced folder would push all of it through the sync client. Override
with `NXF_WORK` to place it on a specific disk. Clear it between `make all` and
`make q10` if space is tight.
