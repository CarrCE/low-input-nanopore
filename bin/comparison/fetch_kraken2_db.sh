#!/usr/bin/env bash
#
# Download and verify the pinned classification databases used by the
# prior-study Kraken2 reanalysis, as recorded in
# comparison/kraken2_db.manifest.tsv.
#
#   Kraken2 PlusPF-8, 2024-12-28 build   ~5.5 GiB compressed, ~7.5 GiB unpacked
#   NCBI new_taxdump, 2025-01-01          ~133 MiB
#
# Both are the exact artefacts that epi2me-labs/wf-metagenomics v2.14.1 resolves
# `--database_set PlusPF-8` to. Fetching them here, verifying them, and then
# passing the local paths to the workflow (see run_kraken2_reanalysis.sh) means
# the classification cannot silently drift onto a different database build.
#
# Usage:
#     bash comparison/fetch_kraken2_db.sh              # download if needed, verify, unpack
#     DEST=/path/to/dbs bash comparison/fetch_kraken2_db.sh
#     KEEP_ARCHIVE=1 bash comparison/fetch_kraken2_db.sh    # keep the 5.5 GiB tarball
#
# The default destination is data/kraken2_db/, matching where the read data and
# fetch_raghavendra.sh already put their downloads. That is ~7.5 GiB unpacked;
# if this clone lives in a cloud-synced folder, point DEST
# somewhere unsynced and pass the same path to run_kraken2_reanalysis.sh --db.
#
# Idempotent: re-running once everything is present and verified touches neither
# the network nor the unpacked files.
#
# Requirements: bash, curl, tar, unzip, and one of md5sum / md5.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"

MANIFEST="${SCRIPT_DIR}/kraken2_db.manifest.tsv"
DEST="${DEST:-${REPO_ROOT}/data/kraken2_db}"
CACHE_DIR="${DEST}/.archive"

DB_DIR="${DEST}/pluspf8_20241228"
TAXONOMY_ZIP="${DEST}/new_taxdump_2025-01-01.zip"

log() { printf '[fetch-kraken2-db] %s\n' "$*"; }
die() { printf '[fetch-kraken2-db] error: %s\n' "$*" >&2; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || die "'$1' is required but not on PATH"; }
need curl
need tar
need unzip

[ -f "${MANIFEST}" ] || die "manifest not found: ${MANIFEST}"

md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1; then md5 -q "$1"
    else echo ""; fi
}

# ---------------------------------------------------------------------------
# Read one field of one manifest row. The manifest is tab-separated with '#'
# comments; awk is given the key and the 1-based column so the field order lives
# in exactly one place (the header) rather than being duplicated here.
# ---------------------------------------------------------------------------
field() {
    awk -F'\t' -v key="$1" -v col="$2" '
        /^#/ { next }
        $1 == key { print $col; found = 1; exit }
        END { if (!found) exit 3 }
    ' "${MANIFEST}" || die "manifest has no entry for key '$1'"
}

# ---------------------------------------------------------------------------
# Download <url> to <path>, then check size and MD5 against the manifest. A file
# that fails either check is left in place for inspection but the script stops:
# silently redownloading would mask a URL that has started serving something
# else, which is the exact failure this pin exists to catch.
# ---------------------------------------------------------------------------
fetch_verified() {
    local key="$1" path="$2" label="$3"
    local url bytes want_md5 actual_bytes actual_md5

    url="$(field "${key}" 3)"
    bytes="$(field "${key}" 4)"
    want_md5="$(field "${key}" 5)"

    if [ -s "${path}" ]; then
        log "${label}: archive already present at ${path}"
    else
        log "${label}: downloading from ${url}"
        curl -fL --progress-bar --retry 5 --retry-delay 5 --retry-connrefused -C - \
             -o "${path}" "${url}" \
          || die "${label}: download failed; rerun to resume"
    fi

    actual_bytes="$(wc -c < "${path}" | tr -d ' ')"
    [ "${actual_bytes}" = "${bytes}" ] \
        || die "${label}: size mismatch -- expected ${bytes} bytes, got ${actual_bytes}. Delete '${path}' and rerun."

    actual_md5="$(md5_of "${path}")"
    if [ -n "${actual_md5}" ]; then
        [ "${actual_md5}" = "${want_md5}" ] \
            || die "${label}: MD5 mismatch -- expected ${want_md5}, got ${actual_md5}. The URL is no longer serving the pinned build; do not use it."
        log "${label}: MD5 verified (${want_md5})"
    else
        log "${label}: warning -- no md5sum/md5 on PATH; size checked but checksum not verified"
    fi
}

mkdir -p "${DEST}" "${CACHE_DIR}"

# ---- Kraken2 PlusPF-8 ------------------------------------------------------
# hash.k2d is the large one and is what kraken2 memory-maps; its presence is the
# marker for "already unpacked".
if [ -s "${DB_DIR}/hash.k2d" ] && [ -s "${DB_DIR}/taxo.k2d" ] && [ -s "${DB_DIR}/opts.k2d" ]; then
    log "PlusPF-8: already unpacked in ${DB_DIR}"
else
    DB_TAR="${CACHE_DIR}/k2_pluspf_08gb_20241228.tar.gz"
    log "PlusPF-8: this is a ~5.5 GiB download, cached in ${CACHE_DIR}"
    fetch_verified kraken2_pluspf8 "${DB_TAR}" "PlusPF-8"

    # Unpack into a scratch dir and move into place only on success, so an
    # interrupted extraction never leaves a half-populated database that the
    # "already unpacked" check above would accept on the next run.
    TMP_DIR="$(mktemp -d "${CACHE_DIR}/unpack.XXXXXX")"
    trap 'rm -rf "${TMP_DIR}"' EXIT
    log "PlusPF-8: unpacking (~7.5 GiB)"
    tar -xzf "${DB_TAR}" -C "${TMP_DIR}"
    [ -s "${TMP_DIR}/hash.k2d" ] || die "PlusPF-8: hash.k2d missing after unpack"
    rm -rf "${DB_DIR}"
    mv "${TMP_DIR}" "${DB_DIR}"
    trap - EXIT
    chmod -R u+rwX "${DB_DIR}"

    if [ "${KEEP_ARCHIVE:-0}" = "1" ]; then
        log "PlusPF-8: keeping ${DB_TAR} (KEEP_ARCHIVE=1)"
    else
        rm -f "${DB_TAR}"
        log "PlusPF-8: removed the tarball; set KEEP_ARCHIVE=1 to keep it next time"
    fi
fi

# ---- NCBI taxonomy ---------------------------------------------------------
# wf-metagenomics takes the taxonomy as a zip and unpacks it itself, so this one
# stays archived.
fetch_verified ncbi_taxdump "${TAXONOMY_ZIP}" "taxdump 2025-01-01"

# ---- report ----------------------------------------------------------------
cat <<EOF

[fetch-kraken2-db] ready.

  database  ${DB_DIR}
  taxonomy  ${TAXONOMY_ZIP}

Pass these to the workflow with:

  bash comparison/run_kraken2_reanalysis.sh --fastq <dir> --out <dir> --min-qual 10

which reads the same manifest and refuses to run against anything else.
EOF
