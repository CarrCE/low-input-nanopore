# ---------------------------------------------------------------------------
# low-input-nanopore
#
# Convenience wrapper around `./run.sh` (which itself wraps `nextflow run`; see
# README.md for why that indirection exists) and the two local container builds.
#
# Everything is quoted because this repository is routinely cloned into a path
# containing spaces (".../My Research Data/..."). Target names
# never contain paths, so make itself stays happy; the recipes quote $(ROOT).
#
# Kept compatible with the GNU Make 3.81 that ships with macOS: no .ONESHELL,
# so every recipe line stands alone as its own shell command.
# ---------------------------------------------------------------------------

SHELL       := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c   # honoured by make >= 3.82; ignored by 3.81

ROOT := $(CURDIR)

TOOLS_IMAGE    ?= low-input-nanopore/tools:0.1.0
ANALYSIS_IMAGE ?= low-input-nanopore/analysis:0.1.0
SCRUBBER_IMAGE ?= low-input-nanopore/scrubber:0.1.0

# Attribution floor; keep in step with params.min_aln_frac in nextflow.config.
MIN_ALN_FRAC ?= 0.10

# Container engine profile passed to Nextflow; override with `make s2 PROFILE=singularity`.
PROFILE ?= docker
# Extra arguments appended to every pipeline run, e.g. `make test NF_ARGS=-resume`.
NF_ARGS ?=

# Smoke-test subsample (see assets/samplesheets/test.csv).
DEMO_SOURCE ?= data/lowinput_s2_r1.fastq
DEMO_OUT    ?= data/test/test_s2.fastq
DEMO_READS  ?= 40000

.DEFAULT_GOAL := help

.PHONY: help images test check verify measurements masksummary threshold seqsummary versions runmeta poolcov attribution comparison coverage modedelta estcontrol assigneddepth divergence s1 s2 all raghavendra q10 demo-data clean

help: ## Show this help
	@printf 'low-input-nanopore -- make targets\n\n'
	@grep -hE '^[a-zA-Z0-9_-]+:.*##' "$(firstword $(MAKEFILE_LIST))" \
	  | awk 'BEGIN {FS = ":.*##[ ]*"} {printf "  %-11s %s\n", $$1, $$2}'
	@printf '\nVariables: PROFILE=%s  NF_ARGS=%s\n' '$(PROFILE)' '$(NF_ARGS)'
	@printf 'Images:    %s  %s\n' '$(TOOLS_IMAGE)' '$(ANALYSIS_IMAGE)'

images: ## Build the local Docker images (tools + analysis + scrubber)
	@echo "==> building $(TOOLS_IMAGE)"
	@docker build -t "$(TOOLS_IMAGE)" "$(ROOT)/docker/tools"
	@echo "==> building $(ANALYSIS_IMAGE)"
	@docker build -t "$(ANALYSIS_IMAGE)" "$(ROOT)/docker/analysis"
	@echo "==> building $(SCRUBBER_IMAGE)"
	@# amd64 explicitly: HRRT has no arm64 build, so on Apple Silicon this
	@# image is emulated. Only --mask_human runs need it.
	@docker build --platform linux/amd64 -t "$(SCRUBBER_IMAGE)" "$(ROOT)/docker/scrubber"
	@echo "==> done: $(TOOLS_IMAGE) $(ANALYSIS_IMAGE) $(SCRUBBER_IMAGE)"

test: ## Run the smoke test (40k reads, test profile)
	@"$(ROOT)/run.sh" -profile $(PROFILE),test $(NF_ARGS)

s1: ## Run lowinput_s1 (D6311 log-distributed DNA, 3 replicates)
	@"$(ROOT)/run.sh" -profile $(PROFILE) \
	    --samplesheet assets/samplesheets/lowinput_s1.csv $(NF_ARGS)

s2: ## Run lowinput_s2 (D6321 low-microbial-load cells, r0-r3)
	@"$(ROOT)/run.sh" -profile $(PROFILE) \
	    --samplesheet assets/samplesheets/lowinput_s2.csv $(NF_ARGS)

# THIS is the invocation behind the published display items, and `make s1` /
# `make s2` are not: AGGREGATE only ever sees the samples of the run that
# invokes it, and every target here writes to the same default results/summary.
# Running s2 then s1 therefore leaves results/summary describing lowinput_s1
# alone, and `make test` leaves it describing 40,000 smoke-test reads. --mode
# both is required for Figure S1, which compares the two assignment rules.
all: ## Run both experiments in both modes -- the published analysis
	@"$(ROOT)/run.sh" -profile $(PROFILE) \
	    --samplesheet assets/samplesheets/all.csv --mode both $(NF_ARGS)

# Secondary analyses publish UNDER results/, not alongside it: results/raghavendra
# and results/q10 rather than results_raghavendra/ and results_q10/. One output
# tree per clone is easier to reason about, to archive and to delete, and
# `results/` alone in .gitignore then covers every one of them.
#
# Nesting is safe against every glob in this repository, which was checked
# rather than assumed. Each consumer keys on a <sid>/<sid>.* shape --
# coverage_attribution.py wants <sid>/coverage/<sid>.coverage_summary.tsv,
# assigned_depth.sh and contaminant_divergence.sh want
# <sid>/competitive/<sid>.assignments.tsv.gz, comparison_data.py globs
# */<mode>/*.metrics.tsv, and `make coverage` globs results/*/coverage/*. A
# directory named q10 satisfies none of them, so it is skipped rather than
# mistaken for a replicate.
raghavendra: ## Reanalyse the Basapathi Raghavendra 2023 reads (-> results/raghavendra)
	@"$(ROOT)/run.sh" -profile $(PROFILE) \
	    --samplesheet assets/samplesheets/raghavendra_2023.csv \
	    --outdir results/raghavendra $(NF_ARGS)

q10: ## Quality-matched rerun at Q10 -- a published result (-> results/q10)
	@"$(ROOT)/run.sh" -profile $(PROFILE) \
	    --samplesheet assets/samplesheets/all.csv --mode both \
	    --min_qscore 10 --outdir results/q10 $(NF_ARGS)

# Whole 4-line records only: a plain `head -n` is fine at 40000 but silently
# truncates a record the moment DEMO_READS is set to a value that is not a
# multiple of the record length. Runs inside the analysis container so the
# subsample does not depend on whatever python the host happens to have.
define DEMO_SCRIPT
import os, sys
src, out, want = os.environ["SRC"], os.environ["OUT"], int(os.environ["NREADS"])
n = 0
with open(src) as fi, open(out, "w") as fo:
    while n < want:
        rec = [fi.readline() for _ in range(4)]
        if not rec[0]:
            break
        if not rec[0].startswith("@") or not rec[2].startswith("+"):
            sys.exit("error: record %d of %s is not a well-formed 4-line FASTQ "
                     "record; refusing to write a broken subsample" % (n + 1, src))
        fo.writelines(rec)
        n += 1
print("[demo-data] wrote {:,} reads -> {}".format(n, out))
if n < want:
    sys.stderr.write("warning: {} held only {:,} reads (asked for {:,})\n".format(src, n, want))
endef
export DEMO_SCRIPT

demo-data: ## Regenerate the smoke-test subsample into data/test/
	@test -f "$(ROOT)/$(DEMO_SOURCE)" \
	  || { echo "error: $(DEMO_SOURCE) not found under $(ROOT)" >&2; exit 1; }
	@mkdir -p "$(ROOT)/data/test"
	@printf '%s\n' "$$DEMO_SCRIPT" | docker run --rm -i \
	    -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/work -w /work \
	    -e SRC="$(DEMO_SOURCE)" -e OUT="$(DEMO_OUT)" -e NREADS="$(DEMO_READS)" \
	    "$(ANALYSIS_IMAGE)" python3 -

TEST_BAM        ?= results/test_s2/alignments/test_s2.qname.bam
TEST_CONTIG_MAP ?= results/references/lowinput_s2/contig_map.tsv

check: ## Assert consensus accounting and human masking (needs `make test` first)
	@test -f "$(ROOT)/$(TEST_BAM)" \
	  || { echo "error: $(TEST_BAM) not found. Run 'make test' first." >&2; exit 1; }
	@docker run --rm -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 tests/consensus_accounting.py \
	        --bam "$(TEST_BAM)" --contig-map "$(TEST_CONTIG_MAP)"
	@echo "==> human masking: decision logic"
	@docker run --rm -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 tests/mask_human_logic.py
	@echo "==> human masking: committed read-level fixture"
	@# Needs no mapping and no network; see assets/testdata/README.md.
	@docker run --rm -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 tests/human_masking.py

# Exit 2 means "pendings, no errors" and is tolerated; exit 1 (a real
# inconsistency) still fails the target.
# The evidence behind params.min_aln_frac. Reads a finished run's per-read
# assignments, so it needs no BAM pass and no re-mapping. Regenerate this rather
# than quoting the numbers from memory if the threshold is ever questioned.
threshold: ## Coverage distribution behind the attribution floor (--min_aln_frac)
	@docker run --rm -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/attribution_threshold.py \
	        --assignments 'results/*/competitive/*.assignments.tsv.gz' \
	        --outdir results/summary --chosen $(MIN_ALN_FRAC)

# Only meaningful after a run with --mask_human true, which is what writes the
# per-sample human_stats.json files this reads.
masksummary: ## Per-replicate human masking statistics (Supplementary table)
	@docker run --rm -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/masking_summary.py \
	        --stats 'results/*/human/*.human_stats.json' \
	        --out results/summary/human_masking.tsv

measurements: ## Check assets/measurements.tsv for gaps and inconsistencies
	@docker run --rm -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/check_measurements.py || test $$? -eq 2

verify: ## Assert every display item satisfies docs/display-items.md
	@docker run --rm -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/verify_display_items.py

# `attribution` is a prerequisite, not just documentation: plot_coverage.py reads
# results/summary/coverage_attribution.tsv and exits if it is missing. It is
# cheap (reads summary tables, no BAM pass), so always regenerating it is the
# right trade against a reader hitting a missing-input error on a fresh tree.
coverage: attribution ## Per-replicate coverage figure (Fig. S2)
	@docker run --rm -u "$$(id -u):$$(id -g)" -e MPLCONFIGDIR=/tmp/mpl \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/plot_coverage.py \
	        --summaries results/*/coverage/*.coverage_summary.tsv \
	        --profiles  results/*/coverage/*.coverage_profile.tsv \
	        --attribution results/summary/coverage_attribution.tsv \
	        --outdir results/summary

modedelta: ## Competitive vs sequential assignment figure (Fig. S1)
	@docker run --rm -u "$$(id -u):$$(id -g)" -e MPLCONFIGDIR=/tmp/mpl \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/plot_mode_delta.py \
	        --per-organism results/summary/per_organism.tsv \
	        --outdir results/summary

estcontrol: ## Estimated no-adaptive-sampling control from Mojarro et al. 2019
	@docker run --rm -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/comparison/estimated_control.py --outdir results/comparison

seqsummary: ## Per-replicate yield, read length and read quality (needs finished runs + FASTQs)
	@docker run --rm -u "$$(id -u):$$(id -g)" -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/sequencing_summary.py --out results/summary/sequencing_summary.tsv

runmeta: ## Per-run acquisition id, model and dates, read from the FASTQs (slow: full pass over ~100 GB)
	@docker run --rm -u "$$(id -u):$$(id -g)" -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/run_metadata.py --out results/summary/run_metadata.tsv

assigneddepth: ## Per-base depth of the reads assignment AWARDED to each organism (slow: one FASTQ pass per replicate)
	@bash "$(ROOT)/bin/assigned_depth.sh"

# `assigneddepth` is a prerequisite because pool_coverage.py --depth-kind assigned
# has nothing to pool without it. Safe to force: bin/assigned_depth.sh skips any
# replicate whose output already exists, so on a warm tree this is a no-op loop
# rather than another pass over 110 GB of FASTQ.
poolcov: assigneddepth ## Pooled-across-replicates coverage summary and Figure S3
	@docker run --rm -u "$$(id -u):$$(id -g)" -e MPLCONFIGDIR=/tmp/mpl \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/pool_coverage.py --depth-kind assigned \
	        --out-summary results/summary/pooled_coverage_summary.tsv \
	        --out-profile results/summary/pooled_coverage_profile.tsv
	@docker run --rm -u "$$(id -u):$$(id -g)" -e MPLCONFIGDIR=/tmp/mpl \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/pool_coverage.py --depth-kind alignment \
	        --out-summary results/summary/pooled_alignment_summary.tsv \
	        --out-profile results/summary/pooled_alignment_profile.tsv
	@docker run --rm -u "$$(id -u):$$(id -g)" -e MPLCONFIGDIR=/tmp/mpl \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/plot_pooled_coverage.py \
	        --summary results/summary/pooled_coverage_summary.tsv \
	        --profile results/summary/pooled_coverage_profile.tsv \
	        --alignment-summary results/summary/pooled_alignment_summary.tsv \
	        --outdir  results/summary

attribution: ## Per-replicate alignment vs attributable depth, and the 1x threshold
	@docker run --rm -u "$$(id -u):$$(id -g)" -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/coverage_attribution.py --out results/summary/coverage_attribution.tsv

# The one analysis step outside the Nextflow DAG: breseq over the contaminant
# reads of a finished run, in its own container. It produces the divergence table
# in Supplementary Section S3, so a reproduction that skips it is incomplete.
# It had no target until now, which made it the only step a reader had to know
# existed rather than discover from `make help`. Idempotent: a replicate whose
# output.gd exists is skipped. Slow -- breseq runs under amd64 emulation.
divergence: ## Contaminant divergence from stock MG1655 (Supplementary Section S3; slow, emulated)
	@bash "$(ROOT)/bin/contaminant_divergence.sh"

comparison: ## Prior-work comparison figure (Fig. 3) into results/comparison/
	@docker run --rm -u "$$(id -u):$$(id -g)" -e MPLCONFIGDIR=/tmp/mpl \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/comparison/plot_comparison.py --results-dir results

versions: ## Record the software the built images actually contain
	@"$(ROOT)/bin/software_versions.sh" --out "$(ROOT)/results/summary/software_versions.tsv"

clean: ## Remove Nextflow scratch (work/, .nextflow*); leaves results/ alone
	@rm -rf "$(ROOT)/work" "$(ROOT)/.nextflow"
	@rm -f "$(ROOT)"/.nextflow.log*
	@echo "==> removed work/ and .nextflow* (results/ and refs/ untouched)"
