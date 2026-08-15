#!/usr/bin/env bash
#
# Launcher for the low-input-nanopore pipeline.
#
# Why this exists: Nextflow cannot run from a project path containing spaces.
# It emits `export PATH=...` and an inner `bash <path>` without quoting, so a
# path like ".../My Research Data/..." breaks the task wrapper before any
# command runs. That is a Nextflow limitation, not something the
# pipeline can fix internally.
#
# When the repository sits at a path with spaces, this script transparently
# creates a space-free symlink outside the repo and launches Nextflow through
# it. Nextflow then generates space-free PATH and work-dir references, while
# still bind-mounting the real (spaced) location for staged inputs -- which it
# escapes correctly.
#
# Usage:
#   ./run.sh --samplesheet assets/samplesheets/lowinput_s2.csv
#   ./run.sh -profile docker,test
#
# Any arguments are passed through to `nextflow run`.

set -euo pipefail

REPO_REAL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

if [[ "$REPO_REAL" == *" "* ]]; then
    # Derive a stable, collision-free link name from the real path so repeated
    # runs reuse the same location and `-resume` keeps working.
    HASH="$(printf '%s' "$REPO_REAL" | shasum -a 256 | cut -c1-12)"
    LINK_BASE="${TMPDIR:-/tmp}"
    LINK_BASE="${LINK_BASE%/}"
    LINK="${LINK_BASE}/low-input-nanopore-${HASH}"

    if [[ "$LINK" == *" "* ]]; then
        echo "error: TMPDIR itself contains spaces (${LINK_BASE})." >&2
        echo "       Set TMPDIR to a space-free directory and re-run." >&2
        exit 1
    fi

    ln -sfn "$REPO_REAL" "$LINK"

    echo "note: project path contains spaces; launching via ${LINK}"
    echo "      (Nextflow cannot handle spaces in the project path)"

    # The work directory deliberately defaults OUTSIDE the repository. A clone
    # often lives in a cloud-synced folder, and Nextflow's work/ holds tens of
    # GB of intermediate BAMs -- putting it in the repo would push all of that
    # through the sync client. Override with NXF_WORK if you want it
    # somewhere specific (e.g. a large external scratch disk).
    WORKDIR="${NXF_WORK:-${LINK}-work}"
    mkdir -p "$WORKDIR"
    echo "      work dir: ${WORKDIR}"
    cd "$LINK"
    exec nextflow run "${LINK}/main.nf" -w "$WORKDIR" "$@"
fi

# Same rationale as above: keep heavy intermediates out of a synced repo.
TMPBASE="${TMPDIR:-/tmp}"
TMPBASE="${TMPBASE%/}"
WORKDIR="${NXF_WORK:-${TMPBASE}/low-input-nanopore-work}"
mkdir -p "$WORKDIR"
echo "note: work dir: ${WORKDIR}"
cd "$REPO_REAL"
exec nextflow run "${REPO_REAL}/main.nf" -w "$WORKDIR" "$@"
