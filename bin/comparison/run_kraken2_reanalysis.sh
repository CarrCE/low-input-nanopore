#!/usr/bin/env bash
#
# Headless re-run of the Kraken2 classification behind the `kraken2_q1` and
# `kraken2_q10` rows in prior_studies.tsv.
#
# The original numbers were produced through the EPI2ME desktop application,
# which records its parameters only inside the run directory it creates. This
# script is the command-line equivalent, written so the provenance lives in
# version control instead: pinned workflow revision, pinned database, pinned
# taxonomy, explicit quality threshold.
#
#   workflow    epi2me-labs/wf-metagenomics v2.14.1
#               commit a57ff73c22b77c2754b7910cd8d24ab6056ed8cc
#   classifier  kraken2
#   database    Kraken2 PlusPF-8, 2024-12-28 build
#   taxonomy    NCBI new_taxdump 2025-01-01
#
# The database and taxonomy are passed as verified local paths rather than as
# `--database_set PlusPF-8`. Both routes resolve to the same two files on this
# revision -- v2.14.1 hard-codes those URLs in its `database_sets` map, which is
# how the pin in kraken2_db.manifest.tsv was derived -- but the local-path route
# is checksum-verified by fetch_kraken2_db.sh first, so a silently reissued
# upstream build cannot change the result without the fetch step failing.
#
# Usage:
#     bash bin/comparison/run_kraken2_reanalysis.sh --fastq DIR --out DIR [options]
#
#     --fastq DIR       directory of FASTQ files, one sample per subdirectory or
#                       a flat directory of files (wf-metagenomics conventions)
#     --out DIR         output directory
#     --min-qual N      read mean-quality lower limit (default 10)
#     --threads N       threads per task (default 4)
#     --profile NAME    nextflow profile (default standard, i.e. docker)
#     --db DIR          override the unpacked database directory
#     --taxonomy ZIP    override the taxonomy zip
#
# Note on --min-qual: wf-metagenomics `--min_read_qual` filters on the mean
# **Phred** of the read, whereas this repository's own filter uses the ONT
# convention of averaging in error-probability space. The two are not the same
# threshold; see docs/ and docs/comparison.md. The value is passed through
# unchanged so it matches what the prior-study rows were generated with.
#
# Requirements: bash, nextflow, docker, git; the databases fetched by
# bin/comparison/fetch_kraken2_db.sh.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"

WF_REPO="epi2me-labs/wf-metagenomics"
WF_TAG="v2.14.1"
WF_COMMIT="a57ff73c22b77c2754b7910cd8d24ab6056ed8cc"

DB_DEFAULT="${REPO_ROOT}/data/kraken2_db/pluspf8_20241228"
TAX_DEFAULT="${REPO_ROOT}/data/kraken2_db/new_taxdump_2025-01-01.zip"

FASTQ=""
OUT=""
MIN_QUAL=10
THREADS=4
PROFILE="standard"
DB="${DB_DEFAULT}"
TAXONOMY="${TAX_DEFAULT}"

log() { printf '[kraken2-reanalysis] %s\n' "$*"; }
die() { printf '[kraken2-reanalysis] error: %s\n' "$*" >&2; exit 1; }

while [ $# -gt 0 ]; do
    case "$1" in
        --fastq)     FASTQ="$2";    shift 2 ;;
        --out)       OUT="$2";      shift 2 ;;
        --min-qual)  MIN_QUAL="$2"; shift 2 ;;
        --threads)   THREADS="$2";  shift 2 ;;
        --profile)   PROFILE="$2";  shift 2 ;;
        --db)        DB="$2";       shift 2 ;;
        --taxonomy)  TAXONOMY="$2"; shift 2 ;;
        -h|--help)   sed -n '2,45p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)           die "unknown argument: $1" ;;
    esac
done

[ -n "${FASTQ}" ] || die "--fastq is required"
[ -n "${OUT}" ]   || die "--out is required"
[ -d "${FASTQ}" ] || die "--fastq is not a directory: ${FASTQ}"

command -v nextflow >/dev/null 2>&1 || die "'nextflow' is required but not on PATH"

# ---------------------------------------------------------------------------
# Databases must exist and must look like what the manifest describes. A missing
# bracken distribution is the failure worth naming explicitly: wf-metagenomics
# accepts a custom database directory only if it carries one, and the error it
# raises otherwise is opaque.
# ---------------------------------------------------------------------------
[ -d "${DB}" ] \
    || die "database directory not found: ${DB}
       Run: bash bin/comparison/fetch_kraken2_db.sh"
[ -s "${DB}/hash.k2d" ] || die "not a Kraken2 database (no hash.k2d): ${DB}"
ls "${DB}"/database*mers.kmer_distrib >/dev/null 2>&1 \
    || die "database has no bracken distribution (database*mers.kmer_distrib) in ${DB}"
[ -s "${TAXONOMY}" ] \
    || die "taxonomy not found: ${TAXONOMY}
       Run: bash bin/comparison/fetch_kraken2_db.sh"

# ---------------------------------------------------------------------------
# Pin the workflow revision, then confirm it. `nextflow pull -r <tag>` follows a
# tag, and a tag can be moved upstream; checking the resolved commit against the
# recorded SHA turns that from a silent substitution into a stop.
# ---------------------------------------------------------------------------
log "pulling ${WF_REPO} ${WF_TAG}"
nextflow pull "${WF_REPO}" -r "${WF_TAG}" >/dev/null

ASSET_DIR="${NXF_ASSETS:-${NXF_HOME:-${HOME}/.nextflow}/assets}/${WF_REPO}"
if [ -d "${ASSET_DIR}/.git" ] && command -v git >/dev/null 2>&1; then
    resolved="$(git -C "${ASSET_DIR}" rev-parse "${WF_TAG}^{commit}" 2>/dev/null || echo "")"
    if [ -n "${resolved}" ]; then
        [ "${resolved}" = "${WF_COMMIT}" ] \
            || die "workflow revision mismatch: ${WF_TAG} resolves to ${resolved}, expected ${WF_COMMIT}. The upstream tag has moved; do not use this run for the published comparison."
        log "workflow commit verified (${WF_COMMIT})"
    else
        log "warning: could not resolve ${WF_TAG} in ${ASSET_DIR}; commit not verified"
    fi
else
    log "warning: ${ASSET_DIR} is not a git checkout; commit not verified"
fi

mkdir -p "${OUT}"

log "database  ${DB}"
log "taxonomy  ${TAXONOMY}"
log "min_read_qual ${MIN_QUAL}, threads ${THREADS}, profile ${PROFILE}"

# --analyse_unclassified is left off deliberately: the prior-study rows count
# classified hits only, and turning it on changes the denominator.
nextflow run "${WF_REPO}" \
    -r "${WF_TAG}" \
    -profile "${PROFILE}" \
    --fastq "${FASTQ}" \
    --classifier kraken2 \
    --database "${DB}" \
    --taxonomy "${TAXONOMY}" \
    --min_read_qual "${MIN_QUAL}" \
    --include_read_assignments \
    --threads "${THREADS}" \
    --out_dir "${OUT}" \
    -w "${OUT}/work"

# ---------------------------------------------------------------------------
# Record what ran alongside the results. The workflow writes its own params.json,
# but not the database identity -- that is the thing this whole exercise exists
# to keep, so it is written explicitly.
# ---------------------------------------------------------------------------
cat > "${OUT}/reanalysis_provenance.txt" <<EOF
workflow        ${WF_REPO} ${WF_TAG}
commit          ${WF_COMMIT}
classifier      kraken2
database        ${DB}
taxonomy        ${TAXONOMY}
min_read_qual   ${MIN_QUAL}
fastq           ${FASTQ}
manifest        assets/comparison/kraken2_db.manifest.tsv
EOF

log "done. Provenance written to ${OUT}/reanalysis_provenance.txt"
log "abundance table:      ${OUT}/abundance_table_species.tsv"
log "per-read assignments: ${OUT}/reads_assignments/<alias>_lineages.kraken2.assignments.tsv"
