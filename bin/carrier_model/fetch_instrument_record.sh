#!/usr/bin/env bash
#
# Download the MinKNOW instrument-data record (doi:10.5281/zenodo.23907136) into
# data/instrument/<sample>/, the layout bin/carrier_model/run_carrier_model.sh
# reads. Zenodo stores the files flat; each name carries its replicate, so the
# folder is taken from the manifest. Every file is checked against the MD5 in
# assets/instrument_files.tsv; a mismatch is deleted and reported. Files that
# are already present with the right MD5 are not downloaded again.
#
# Usage: bash bin/carrier_model/fetch_instrument_record.sh [outdir]
# Requirements: bash, curl, md5sum or md5.

set -euo pipefail

REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
OUT="${1:-${REPO}/data/instrument}"
RECORD="${ZENODO_RECORD:-23907136}"
MANIFEST="${REPO}/assets/instrument_files.tsv"

log() { printf '[instrument-record] %s\n' "$*"; }
md5_of () { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }

n=0; bad=0
while IFS=$'\t' read -r sample file _bytes _rows want _rest; do
    name="$(basename "${file}")"
    dest="${OUT}/${sample}/${name}"
    mkdir -p "${OUT}/${sample}"
    if [ -f "${dest}" ] && [ "$(md5_of "${dest}")" = "${want}" ]; then
        log "have   ${sample}/${name}"; n=$((n + 1)); continue
    fi
    log "fetch  ${sample}/${name}"
    curl -fL --retry 5 --retry-delay 10 -C - -o "${dest}.part" \
        "https://zenodo.org/records/${RECORD}/files/${name}?download=1"
    if [ "$(md5_of "${dest}.part")" = "${want}" ]; then
        mv "${dest}.part" "${dest}"; n=$((n + 1))
    else
        rm -f "${dest}.part"; bad=$((bad + 1))
        log "error: ${name} does not match its MD5 in assets/instrument_files.tsv"
    fi
done < <(grep -v '^#' "${MANIFEST}" | tail -n +2)

log "${n} files verified in ${OUT}${bad:+; ${bad} failed}"
[ "${bad}" -eq 0 ]
