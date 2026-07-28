#!/usr/bin/env bash
#
# Fetch the 10 fastq_pass files of Basapathi Raghavendra et al. 2023 and stage
# them under data/raghavendra_2023/ with the aliases used by
# assets/samplesheets/raghavendra_2023.csv.
#
#   Basapathi Raghavendra J, Zorzano M-P, Kumaresan D, Martin-Torres J.
#   "DNA sequencing at the picogram level to investigate life on Mars and Earth."
#   Sci Rep 13:15277 (2023). doi:10.1038/s41598-023-42170-6
#   Data: https://zenodo.org/records/8208597  (doi:10.5281/zenodo.8208597)
#
# The deposit is a single 771 MiB zip ("MinION low detectability.zip") holding
# fast5 + fastq for every run. Only the 10 fastq_pass files listed below are
# extracted; the rest of the archive is left alone.
#
# Usage:
#     bash comparison/fetch_raghavendra.sh            # download if needed, extract
#     KEEP_ARCHIVE=1 bash comparison/fetch_raghavendra.sh   # keep the zip afterwards
#     ZIP=/path/to/local.zip bash comparison/fetch_raghavendra.sh  # use a local copy
#
# Idempotent: re-running when all 10 FASTQs are already present and valid is a
# no-op and does not touch the network.
#
# Requirements: bash, curl, unzip, gzip, and one of md5sum / md5 (optional).

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths are derived from the script location so this works from any cwd, and
# every expansion is quoted because the repository path contains spaces.
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
OUT_DIR="${REPO_ROOT}/data/raghavendra_2023"
CACHE_DIR="${OUT_DIR}/.archive"

ZENODO_RECORD="https://zenodo.org/records/8208597"
ZIP_URL="https://zenodo.org/records/8208597/files/MinION%20low%20detectability.zip?download=1"
ZIP_NAME="MinION low detectability.zip"
ZIP_MD5="0a426791b2ecf4f9819f71effd9bbaf0"   # from the Zenodo record metadata
ZIP_BYTES="808001111"

# ---------------------------------------------------------------------------
# alias <TAB> basename-inside-the-archive
#
# Matched inside the zip as "*/fastq_pass/<basename>" rather than by full path:
# the archive's directory names are awkward ("10pg yeast_2pg ec.oli") and the
# basenames are already unique, the "_pass_" infix distinguishing them from the
# fastq_fail files.
# ---------------------------------------------------------------------------
FILES="\
Ec_R1	ANP471_pass_4d26e51d_78259c2a_0.fastq.gz
Ec_R2	ANP315_pass_cd07c6d0_7c99cdb7_0.fastq.gz
Ec_R3	ANR193_pass_5181599a_4b0ad359_0.fastq.gz
YSC_R1	ANR915_pass_9147095f_57891227_0.fastq.gz
YSC_R2	ANQ439_pass_14b5c719_82187b38_0.fastq.gz
YSC_R3	ANP378_pass_301b7512_72027704_0.fastq.gz
Mix1_R1	APM812_pass_b4817b78_5e181096_0.fastq.gz
Mix1_R2	AQD920_pass_60a8e58c_0522ae15_0.fastq.gz
Mix2_R1	APU620_pass_3e0bfa22_1e2b9d5d_0.fastq.gz
Mix2_R2	APN068_pass_0a186932_3420adab_0.fastq.gz"

EXPECTED_N=10

log() { printf '[fetch-raghavendra] %s\n' "$*"; }
die() { printf '[fetch-raghavendra] error: %s\n' "$*" >&2; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || die "'$1' is required but not on PATH"; }
need curl
need unzip
need gzip

# ---------------------------------------------------------------------------
# Already done?  A FASTQ counts as present only if gzip can validate it, so a
# truncated extraction from an interrupted run is redownloaded rather than fed
# to the pipeline.
# ---------------------------------------------------------------------------
have_all=1
while IFS="$(printf '\t')" read -r alias_name member; do
    [ -n "${alias_name}" ] || continue
    target="${OUT_DIR}/${alias_name}.fastq.gz"
    if [ ! -s "${target}" ] || ! gzip -t "${target}" 2>/dev/null; then
        have_all=0
        break
    fi
done <<EOF
${FILES}
EOF

if [ "${have_all}" -eq 1 ]; then
    log "all ${EXPECTED_N} FASTQs already present and valid in ${OUT_DIR} -- nothing to do"
    exit 0
fi

mkdir -p "${OUT_DIR}" "${CACHE_DIR}"

# ---------------------------------------------------------------------------
# Locate or download the archive.  ZIP=... lets a caller point at a copy they
# already have; otherwise it is downloaded into the cache directory and reused
# on subsequent runs.  curl -C - resumes a partial download (Zenodo serves
# ranged requests).
# ---------------------------------------------------------------------------
ZIP="${ZIP:-${CACHE_DIR}/${ZIP_NAME}}"

if [ -s "${ZIP}" ]; then
    log "using existing archive: ${ZIP}"
else
    log "downloading ${ZIP_NAME} (~771 MiB) from ${ZENODO_RECORD}"
    log "this is a large one-off download; it is cached in ${CACHE_DIR}"
    curl -fL --retry 5 --retry-delay 5 --retry-connrefused -C - \
         -o "${ZIP}" "${ZIP_URL}" \
      || die "download failed; rerun to resume, or set ZIP=/path/to/${ZIP_NAME}"
fi

# ---- integrity -------------------------------------------------------------
actual_bytes="$(wc -c < "${ZIP}" | tr -d ' ')"
if [ "${actual_bytes}" != "${ZIP_BYTES}" ]; then
    die "archive size mismatch: expected ${ZIP_BYTES} bytes, got ${actual_bytes}. Delete '${ZIP}' and rerun."
fi

md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1; then md5 -q "$1"
    else echo ""; fi
}
actual_md5="$(md5_of "${ZIP}")"
if [ -n "${actual_md5}" ]; then
    [ "${actual_md5}" = "${ZIP_MD5}" ] \
        || die "archive MD5 mismatch: expected ${ZIP_MD5}, got ${actual_md5}. Delete '${ZIP}' and rerun."
    log "archive MD5 verified (${ZIP_MD5})"
else
    log "warning: no md5sum/md5 on PATH; size checked but checksum not verified"
fi

# ---------------------------------------------------------------------------
# Extract the 10 members, flattened (-j), into a scratch dir, then move each to
# its alias.  Extracting to scratch first keeps a failed run from leaving
# half-named files in data/raghavendra_2023/.
# ---------------------------------------------------------------------------
WORK_DIR="$(mktemp -d "${CACHE_DIR}/extract.XXXXXX")"
cleanup() { rm -rf "${WORK_DIR}"; }
trap cleanup EXIT

n_ok=0
while IFS="$(printf '\t')" read -r alias_name member; do
    [ -n "${alias_name}" ] || continue
    target="${OUT_DIR}/${alias_name}.fastq.gz"

    if [ -s "${target}" ] && gzip -t "${target}" 2>/dev/null; then
        log "${alias_name}: already present, skipping"
        n_ok=$(( n_ok + 1 ))
        continue
    fi

    # -o overwrite, -j flatten, -q quiet.  The pattern is quoted so the shell
    # leaves the '*' for unzip to expand against the archive.
    unzip -o -j -q "${ZIP}" "*/fastq_pass/${member}" -d "${WORK_DIR}" \
        || die "could not extract '${member}' from '${ZIP}'"

    [ -s "${WORK_DIR}/${member}" ] \
        || die "'${member}' not found in the archive (expected at */fastq_pass/${member})"
    gzip -t "${WORK_DIR}/${member}" \
        || die "'${member}' extracted but is not a valid gzip stream"

    mv -f "${WORK_DIR}/${member}" "${target}"
    log "${alias_name}: $(basename -- "${member}") -> $(basename -- "${target}")"
    n_ok=$(( n_ok + 1 ))
done <<EOF
${FILES}
EOF

# ---------------------------------------------------------------------------
# Verify: exactly EXPECTED_N valid FASTQs, and no strays in the output dir.
# ---------------------------------------------------------------------------
[ "${n_ok}" -eq "${EXPECTED_N}" ] \
    || die "staged ${n_ok} files but expected ${EXPECTED_N}"

n_present=0
while IFS="$(printf '\t')" read -r alias_name member; do
    [ -n "${alias_name}" ] || continue
    target="${OUT_DIR}/${alias_name}.fastq.gz"
    [ -s "${target}" ] || die "missing after extraction: ${target}"
    gzip -t "${target}" || die "corrupt after extraction: ${target}"
    n_present=$(( n_present + 1 ))
done <<EOF
${FILES}
EOF

[ "${n_present}" -eq "${EXPECTED_N}" ] \
    || die "verification found ${n_present} FASTQs, expected ${EXPECTED_N}"

n_found="$(find "${OUT_DIR}" -maxdepth 1 -name '*.fastq.gz' -type f | wc -l | tr -d ' ')"
if [ "${n_found}" != "${EXPECTED_N}" ]; then
    log "warning: ${OUT_DIR} holds ${n_found} *.fastq.gz files, expected exactly ${EXPECTED_N}"
    log "warning: unexpected files will not be read by the samplesheet, but check for stale copies"
fi

if [ "${KEEP_ARCHIVE:-0}" = "1" ]; then
    log "keeping archive at ${ZIP} (KEEP_ARCHIVE=1)"
else
    if [ "${ZIP}" = "${CACHE_DIR}/${ZIP_NAME}" ]; then
        rm -f "${ZIP}"
        log "removed cached archive (set KEEP_ARCHIVE=1 to keep it)"
    fi
fi

log "OK: ${EXPECTED_N} FASTQs staged in ${OUT_DIR}"
log "next: ./run.sh -profile docker --samplesheet assets/samplesheets/raghavendra_2023.csv"
