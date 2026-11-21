#!/usr/bin/env bash
#
# Does a lambda read's start type change how the live mapper saw it?
# (Supplementary Section S6; see bin/carrier_model/decision_edge.py for what is
# being tested and why, including the verdict rule, which is fixed in advance.)
#
# Tests whether cos-start lambda-derived reads mapped worse at the adaptive
# sampling decision than break-start ones. If they did, the rise of f_cos with
# read length is a readout of the live mapper at the reference edge rather than
# evidence about carrier intactness.
#
# Inputs
#   - a finished `make carriermodel` workdir:
#     results/carrier_model/work/<sample>.{read_ids.txt.gz,reads.npz,
#     instrument.npz,lambda_ids.tsv,ends.paf.gz}
#   - the MinKNOW instrument records (doi:10.5281/zenodo.23907136), one folder
#     per replicate holding AS_timings_<sample>.*.gz. `make instrumentdata`
#     fetches and checks them.
# Nothing is remapped and no read is rebasecalled: this reads the existing
# workdir and the instrument records only.
#
# Usage:  bash bin/carrier_model/run_decision_edge.sh
#         RESULTS=... INSTRUMENT=... OUT=... bash bin/carrier_model/run_decision_edge.sh
#
# Requirements: bash, docker, and the images from `make images`.

set -euo pipefail

REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
RESULTS="${RESULTS:-${REPO}/results}"
INSTRUMENT="${INSTRUMENT:-${REPO}/data/instrument}"
ANALYSIS_IMAGE="${ANALYSIS_IMAGE:-low-input-nanopore/analysis:0.1.0}"
SAMPLES=(lowinput_s1_r1 lowinput_s1_r2 lowinput_s1_r3
         lowinput_s2_r0 lowinput_s2_r1 lowinput_s2_r2 lowinput_s2_r3)

log() { printf '[decisionedge] %s\n' "$*"; }
die() { printf '[decisionedge] error: %s\n' "$*" >&2; exit 1; }

command -v docker >/dev/null 2>&1 || die "'docker' is required but not on PATH"
RESULTS="$(cd "${RESULTS}" 2>/dev/null && pwd -P)" || die "no results directory; run 'make all' first"
[ -d "${INSTRUMENT}" ] || die "no instrument records at ${INSTRUMENT}; run 'make instrumentdata' or set INSTRUMENT=..."
INSTRUMENT="$(cd "${INSTRUMENT}" && pwd -P)"

WORK="${RESULTS}/carrier_model/work"
OUT="${OUT:-${RESULTS}/carrier_model}"
mkdir -p "${OUT}"
OUT="$(cd "${OUT}" && pwd -P)"

for s in "${SAMPLES[@]}"; do
    for f in read_ids.txt.gz reads.npz instrument.npz lambda_ids.tsv ends.paf.gz; do
        [ -s "${WORK}/${s}.${f}" ] || die "${s}: no ${WORK}/${s}.${f}.
  This reads a finished carrier-model workdir; run 'make carriermodel' first."
    done
    [ -d "${INSTRUMENT}/${s}" ] || die "${s}: no instrument records in ${INSTRUMENT}/${s}"
done

log "comparing cos-start with break-start reads at the live decision"
docker run --rm -u "$(id -u):$(id -g)" \
    -v "${REPO}":/repo:ro -v "${RESULTS}":/results -v "${INSTRUMENT}":/instrument:ro \
    -w /repo "${ANALYSIS_IMAGE}" \
    python3 bin/carrier_model/decision_edge.py \
    /results/carrier_model/work /instrument "${OUT/#${RESULTS}//results}" \
    "${SAMPLES[@]}"
log "done: ${OUT#${REPO}/}/decision_edge.json"
