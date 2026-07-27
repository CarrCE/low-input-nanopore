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

# Container engine profile passed to Nextflow; override with `make s2 PROFILE=singularity`.
PROFILE ?= docker
# Extra arguments appended to every pipeline run, e.g. `make test NF_ARGS=-resume`.
NF_ARGS ?=

# Smoke-test subsample (see assets/samplesheets/test.csv).
DEMO_SOURCE ?= data/lowinput_s2_r1.fastq
DEMO_OUT    ?= data/test/test_s2.fastq
DEMO_READS  ?= 40000

.DEFAULT_GOAL := help

.PHONY: help images test check measurements seqsummary versions s1 s2 demo-data clean

help: ## Show this help
	@printf 'low-input-nanopore -- make targets\n\n'
	@grep -hE '^[a-zA-Z0-9_-]+:.*##' "$(firstword $(MAKEFILE_LIST))" \
	  | awk 'BEGIN {FS = ":.*##[ ]*"} {printf "  %-11s %s\n", $$1, $$2}'
	@printf '\nVariables: PROFILE=%s  NF_ARGS=%s\n' '$(PROFILE)' '$(NF_ARGS)'
	@printf 'Images:    %s  %s\n' '$(TOOLS_IMAGE)' '$(ANALYSIS_IMAGE)'

images: ## Build the two local Docker images (tools + analysis)
	@echo "==> building $(TOOLS_IMAGE)"
	@docker build -t "$(TOOLS_IMAGE)" "$(ROOT)/docker/tools"
	@echo "==> building $(ANALYSIS_IMAGE)"
	@docker build -t "$(ANALYSIS_IMAGE)" "$(ROOT)/docker/analysis"
	@echo "==> done: $(TOOLS_IMAGE) $(ANALYSIS_IMAGE)"

test: ## Run the smoke test (40k reads, test profile)
	@"$(ROOT)/run.sh" -profile $(PROFILE),test $(NF_ARGS)

s1: ## Run lowinput_s1 (D6311 log-distributed DNA, 3 replicates)
	@"$(ROOT)/run.sh" -profile $(PROFILE) \
	    --samplesheet assets/samplesheets/lowinput_s1.csv $(NF_ARGS)

s2: ## Run lowinput_s2 (D6321 low-microbial-load cells, r0-r3)
	@"$(ROOT)/run.sh" -profile $(PROFILE) \
	    --samplesheet assets/samplesheets/lowinput_s2.csv $(NF_ARGS)

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

check: ## Assert consensus subtraction preserves read accounting (needs `make test` first)
	@test -f "$(ROOT)/$(TEST_BAM)" \
	  || { echo "error: $(TEST_BAM) not found. Run 'make test' first." >&2; exit 1; }
	@docker run --rm -u "$$(id -u):$$(id -g)" \
	    -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 tests/consensus_accounting.py \
	        --bam "$(TEST_BAM)" --contig-map "$(TEST_CONTIG_MAP)"

measurements: ## Check assets/measurements.tsv for gaps and inconsistencies
	@python3 "$(ROOT)/bin/check_measurements.py" || test $$? -eq 2

seqsummary: ## Per-replicate yield, read length and read quality (needs finished runs + FASTQs)
	@docker run --rm -u "$$(id -u):$$(id -g)" -v "$(ROOT)":/repo -w /repo "$(ANALYSIS_IMAGE)" \
	    python3 bin/sequencing_summary.py --out results/summary/sequencing_summary.tsv

versions: ## Record the software the built images actually contain
	@"$(ROOT)/bin/software_versions.sh" --out "$(ROOT)/results/summary/software_versions.tsv"

clean: ## Remove Nextflow scratch (work/, .nextflow*); leaves results/ alone
	@rm -rf "$(ROOT)/work" "$(ROOT)/.nextflow"
	@rm -f "$(ROOT)"/.nextflow.log*
	@echo "==> removed work/ and .nextflow* (results/ and refs/ untouched)"
