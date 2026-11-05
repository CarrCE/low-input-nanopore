#!/usr/bin/env bash
#
# Fetch the 9 deposited FASTQs of Zorzano et al. 2025 into data/zorzano_2025/
# and check each against the MD5s below.
#
#   Zorzano M-P et al. "Fragmented deoxyribonucleic acid in rocks and its
#   significance for studying life on Mars and Earth." Communications Earth &
#   Environment (2025).
#   Data: https://github.com/FISABIO-BioinformaticsService/CommsEarth_Environment-fDNArocks
#         (RawData/), pinned to commit 08ff9f2c1df255be55c7451bc1486146066b70b1
#
# The authors deposited their reads in a GitHub repository, not an archive, so
# the clone is pinned to a commit and every file is checked. These are the files
# behind the `kraken2_q1` and `kraken2_q10` Zorzano rows in
# assets/comparison/prior_studies.tsv; they were checked on 4 November 2026 to
# be byte-identical to the files those rows were classified from.
#
# Usage:
#     bash bin/comparison/fetch_zorzano.sh          # clone if needed, verify
#     OUT_DIR=/path bash bin/comparison/fetch_zorzano.sh
#
# Then stage them for the reanalysis:
#     python3 bin/comparison/stage_prior_reads.py zorzano \
#         --src data/zorzano_2025/RawData --out <samples dir>
#
# Idempotent: when the clone exists at the pinned commit and every file
# verifies, nothing touches the network.
#
# Requirements: bash, git, and one of md5sum / md5.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
OUT_DIR="${OUT_DIR:-${REPO_ROOT}/data/zorzano_2025}"

URL="https://github.com/FISABIO-BioinformaticsService/CommsEarth_Environment-fDNArocks.git"
COMMIT="08ff9f2c1df255be55c7451bc1486146066b70b1"

# basename in RawData/ <TAB> MD5
FILES="\
ExtractedWaterControl.fastq.gz	b65a16ca2882da6914ff2c70dc669583
MicrobialitesNonIrradiated.fastq.gz	993b651ce1369c75e16291d2b27903bd
MicrobialitesIrradiated0.874MGy.fastq.gz	14a904e651947e452db4553d6c4f2bf7
MicrobialitesIrradiated10.5MGy.fastq.gz	e50b68dc2bf4340152b225adeac49052
StromatolitesNonIrradiated.fastq.gz	76f99bff042e4466396af16c8a11d484
StromatolitesIrradiated10.5MGy.fastq.gz	8dd1e322258cfc8a1ddff74b9b0a4be2
OxideStromatolitesNonIrradiated.fastq.gz	2fe8a1cec287b80691f16e795587a584
OxideStromatolitesIrradiated.fastq.gz	f2ee9229fe6131b451b51b88dff9196a
Carbonates.fastq.gz	58ef414c4d09ab3e15fc6e12032a1719"

log() { printf '[fetch-zorzano] %s\n' "$*"; }
die() { printf '[fetch-zorzano] error: %s\n' "$*" >&2; exit 1; }
command -v git >/dev/null 2>&1 || die "'git' is required but not on PATH"
md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1; then md5 -q "$1"
    else die "need md5sum or md5"; fi
}

if [ ! -d "${OUT_DIR}/.git" ]; then
    log "cloning ${URL}"
    git clone --quiet "${URL}" "${OUT_DIR}"
fi
head="$(git -C "${OUT_DIR}" rev-parse HEAD)"
if [ "${head}" != "${COMMIT}" ]; then
    log "checking out pinned commit ${COMMIT}"
    git -C "${OUT_DIR}" fetch --quiet origin
    git -C "${OUT_DIR}" checkout --quiet "${COMMIT}" \
        || die "pinned commit ${COMMIT} is not in ${URL}"
fi

bad=0
while IFS=$'\t' read -r name want; do
    f="${OUT_DIR}/RawData/${name}"
    if [ ! -f "${f}" ]; then
        log "missing ${name}"; bad=$((bad + 1)); continue
    fi
    got="$(md5_of "${f}")"
    if [ "${got}" = "${want}" ]; then
        log "ok      ${name}"
    else
        log "MISMATCH ${name}: ${got}, expected ${want}"; bad=$((bad + 1))
    fi
done <<< "${FILES}"
[ "${bad}" -eq 0 ] || die "${bad} file(s) missing or failed their MD5"
log "all 9 files verified in ${OUT_DIR}/RawData"
