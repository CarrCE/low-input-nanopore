#!/usr/bin/env bash
#
# The carrier-depletion model against the seven runs (Supplementary Section S6,
# Figure S1 and Table S11 of the manuscript; see docs/carrier-model.md).
#
# Inputs
#   - the deposited, human-masked reads, looked for in this order:
#       $READS/<sample>.fastq.gz or $READS/<sample>.masked.fastq.gz, if READS is set
#       sra/<sample>.fastq.gz             the --fetch_from_sra download cache
#       results/<sample>/human/<sample>.masked.fastq.gz   a --mask_human run's output
#     The samplesheet's local `fastq` is deliberately NOT used: in the lab's own
#     checkout it holds the reads before human masking, and this analysis is
#     defined on the deposit.
#   - a finished `make all` run: results/<sample>/competitive/*.readlengths.tsv.gz,
#     results/summary/per_sample.tsv and per_organism.tsv, and
#     results/references/<set>/combined.fasta (lambda is taken from it)
#   - the MinKNOW instrument records (doi:10.5281/zenodo.23907136), one folder per
#     replicate: $INSTRUMENT/<sample>/{summary_,AS_decisions_,AS_timings_,
#     pore_activity_}<sample>.*.gz. `make instrumentdata` fetches and checks them.
#
# Outputs
#   results/carrier_model/        analysis tables (decomposition, cos ends,
#                                 decisions, escapes, recapture, time course)
#   results/carrier_model/work/   per-replicate intermediates (reusable)
#   results/summary/carrier-model.{pdf,png,csv,json}       the figure
#   results/summary/carrier-model-table.{csv,json}         the table
#
# Per-replicate steps are skipped when their output exists, so a rerun after a
# change to the analysis repeats only the cheap steps. Delete
# results/carrier_model/work to start over.
#
# Usage:  bash bin/carrier_model/run_carrier_model.sh
#         READS=... RESULTS=... INSTRUMENT=... THREADS=8 bash bin/carrier_model/run_carrier_model.sh
#
# Requirements: bash, docker, and the images from `make images`.

set -euo pipefail

REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
RESULTS="${RESULTS:-${REPO}/results}"
INSTRUMENT="${INSTRUMENT:-${REPO}/data/instrument}"
THREADS="${THREADS:-8}"
TOOLS_IMAGE="${TOOLS_IMAGE:-low-input-nanopore/tools:0.1.0}"
ANALYSIS_IMAGE="${ANALYSIS_IMAGE:-low-input-nanopore/analysis:0.1.0}"
SAMPLES=(lowinput_s1_r1 lowinput_s1_r2 lowinput_s1_r3
         lowinput_s2_r0 lowinput_s2_r1 lowinput_s2_r2 lowinput_s2_r3)
LAMBDA_CONTIG="NC_001416.1"

log() { printf '[carrier-model] %s\n' "$*"; }
die() { printf '[carrier-model] error: %s\n' "$*" >&2; exit 1; }

command -v docker >/dev/null 2>&1 || die "'docker' is required but not on PATH"
RESULTS="$(cd "${RESULTS}" 2>/dev/null && pwd -P)" || die "no results directory; run 'make all' first"
[ -d "${INSTRUMENT}" ] || die "no instrument records at ${INSTRUMENT}; run 'make instrumentdata' or set INSTRUMENT=..."
INSTRUMENT="$(cd "${INSTRUMENT}" && pwd -P)"

OUT="${RESULTS}/carrier_model"
WORK="${OUT}/work"
FIG="${RESULTS}/summary"
mkdir -p "${WORK}" "${FIG}"

# Containers see the repository at /repo, results at /results and the
# instrument records at /instrument, wherever those live on the host.
run_py () {
    docker run --rm -u "$(id -u):$(id -g)" -e MPLCONFIGDIR=/tmp/mpl \
        -v "${REPO}":/repo:ro -v "${RESULTS}":/results -v "${INSTRUMENT}":/instrument:ro \
        -w /repo "${ANALYSIS_IMAGE}" python3 "$@"
}

# The masked reads for one replicate, as a host path; empty if none is found.
masked_reads () {  # masked_reads <sample_id>
    local s="$1" c
    for c in ${READS:+"${READS}/${s}.fastq.gz" "${READS}/${s}.masked.fastq.gz"} \
             "${REPO}/${NXF_SRADIR:-sra}/${s}.fastq.gz" \
             "${RESULTS}/${s}/human/${s}.masked.fastq.gz"; do
        [ -f "${c}" ] && { printf '%s\n' "${c}"; return; }
    done
}

# Header-driven samplesheet lookup, as in bin/assigned_depth.sh.
sheet_col () {  # sheet_col <sample_id> <column_name>
    awk -F, -v r="$1" -v want="$2" '
        /^#/                { next }
        $1 == "sample_id"   { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $1 == r             { if (want in col) print $(col[want]); exit }
    ' "${REPO}"/assets/samplesheets/all.csv
}

for s in "${SAMPLES[@]}"; do
    rl="${RESULTS}/${s}/competitive/${s}.readlengths.tsv.gz"
    [ -f "${rl}" ] || die "${s}: no ${rl}; run 'make all' first"
    [ -d "${INSTRUMENT}/${s}" ] || die "${s}: no instrument records in ${INSTRUMENT}/${s}"

    # 1. one pass over the deposited reads
    if [ ! -s "${WORK}/${s}.reads.npz" ] || [ ! -s "${WORK}/${s}.read_ids.txt.gz" ]; then
        fastq="$(masked_reads "${s}")"
        [ -n "${fastq}" ] || die "${s}: no masked reads found (set READS=, or run with --fetch_from_sra)"
        log "${s}: reading ${fastq}"
        docker run --rm -u "$(id -u):$(id -g)" \
            -v "${REPO}":/repo:ro -v "${RESULTS}":/results -v "$(dirname "${fastq}")":/reads:ro \
            -w /repo "${ANALYSIS_IMAGE}" python3 bin/carrier_model/extract_reads.py "${s}" \
            "/reads/$(basename "${fastq}")" \
            "/results/${s}/competitive/${s}.readlengths.tsv.gz" /results/carrier_model/work
    fi

    # 2. map the sampled lambda-derived read ends to lambda, with the pipeline's minimap2
    if [ ! -s "${WORK}/${s}.ends.paf.gz" ]; then
        set_name="$(basename "$(sheet_col "${s}" reference_set)" .tsv)"
        combined="references/${set_name}/combined.fasta"
        [ -f "${RESULTS}/${combined}" ] || die "${s}: no ${RESULTS}/${combined}"
        log "${s}: mapping read ends to lambda"
        docker run --rm -u "$(id -u):$(id -g)" -v "${RESULTS}":/results -w /results "${TOOLS_IMAGE}" \
            bash -c "set -euo pipefail
                samtools faidx '${combined}' ${LAMBDA_CONTIG} > carrier_model/work/lambda.fa
                minimap2 -c -x map-ont --secondary=no -t ${THREADS} carrier_model/work/lambda.fa \
                    carrier_model/work/${s}.lambda_ends.fa 2> carrier_model/work/${s}.minimap2.log \
                  | gzip -c > carrier_model/work/${s}.ends.paf.gz.tmp
                mv carrier_model/work/${s}.ends.paf.gz.tmp carrier_model/work/${s}.ends.paf.gz
                minimap2 --version > carrier_model/work/minimap2.version"
    fi

    # 3. instrument records onto the deposited reads
    if [ ! -s "${WORK}/${s}.instrument.npz" ]; then
        log "${s}: joining instrument records"
        run_py bin/carrier_model/join_instrument.py "${s}" \
            "/results/${s}/competitive/${s}.readlengths.tsv.gz" \
            "/results/carrier_model/work/${s}.read_ids.txt.gz" "/instrument/${s}" \
            /results/carrier_model/work > "${WORK}/${s}.join.log"
    fi
done

log "cos ends"
run_py bin/carrier_model/cos_ends.py /results/carrier_model/work /results/carrier_model "${SAMPLES[@]}" \
    > "${OUT}/cos_ends.log"
log "instrument decisions, escapes, recapture, pores"
run_py bin/carrier_model/analyze_instrument.py /results/carrier_model/work /results/carrier_model /instrument
run_py bin/carrier_model/recapture_ends.py /results /results/carrier_model/work /results/carrier_model
run_py bin/carrier_model/decision_quality.py /results /results/carrier_model/work /results/carrier_model /instrument
log "model, figure and table"
run_py bin/carrier_model/analyze.py /repo /results /results/carrier_model/work \
    /results/carrier_model /results/summary
log "done: ${FIG#${REPO}/}/carrier-model.pdf, ${FIG#${REPO}/}/carrier-model-table.csv, ${OUT#${REPO}/}/"
