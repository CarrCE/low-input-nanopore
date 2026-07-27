#!/usr/bin/env bash
# Record the software this pipeline actually ran with.
#
# The manuscript needs a software table, and the honest way to build one is to
# ask the images rather than to transcribe the Dockerfiles. A build ARG records
# what was requested; only the image records what was installed. The two differ
# whenever an upstream pin is a moving pointer -- which, in this repository, is
# true in two places (see the `pin` column).
#
# Nothing is pushed to a registry: both project images are built locally from
# pinned source, so they have no registry digest. The local image ID (a digest
# over the image config) is recorded instead. It identifies the image on the
# machine that produced these results; it is not a portable content address, and
# rebuilding on another host will produce a different ID from identical inputs.
# What makes the build reproducible is the pinned source in docker/, not the ID.
#
# Usage:
#     bin/software_versions.sh                                   # to stdout
#     bin/software_versions.sh --out results/summary/software_versions.tsv
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out=""
[[ "${1:-}" == "--out" ]] && out="${2:?--out needs a path}"

# `container *=` and not merely `container`: the breseq block also sets
# `containerOptions`, and a looser pattern captures both lines.
img_for() { awk -F"'" -v l="withLabel: $1" '$0 ~ l,/}/{ if ($0 ~ /container *=/) print $2 }' "$root/conf/base.config"; }
tools_img=$(img_for tools)
anal_img=$(img_for analysis)
breseq_img=$(img_for breseq)

# Ask the image; empty if the image is not built on this host.
inspect() { docker image inspect --format '{{.Id}}' "$1" 2>/dev/null || true; }
run()     { docker run --rm --entrypoint sh "$@" 2>/dev/null || true; }

emit() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5"; }

{
  emit component version image pin note

  # ---- host -------------------------------------------------------------
  nf=$( (nextflow -version 2>/dev/null || true) | awk '/^ *version /{print $2; exit}')
  emit Nextflow "${nf:-not installed}" host ">=23.10.0" \
       "workflow engine; minimum declared in nextflow.config"
  emit Docker "$( (docker --version 2>/dev/null || true) | awk '{gsub(/,/,""); print $3}')" \
       host any "container runtime; any OCI runtime Nextflow supports will do"

  # ---- tools image ------------------------------------------------------
  v=$(run "$tools_img" -c 'minimap2 --version')
  emit minimap2 "${v:-unknown}" "$tools_img" "2.28" "built from source tag v2.28"
  v=$(run "$tools_img" -c 'samtools --version | head -1' | awk '{print $2}')
  emit samtools "${v:-unknown}" "$tools_img" "1.21" "built from source tag 1.21"
  v=$(run "$tools_img" -c 'samtools --version | sed -n 2p' | awk '{print $3}')
  emit htslib "${v:-unknown}" "$tools_img" "1.21" "built from source tag 1.21"
  v=$(run "$tools_img" -c 'seqkit version' | awk '{print $2}')
  emit seqkit "${v:-unknown}" "$tools_img" "2.8.2" "prebuilt release binary"
  v=$(run "$tools_img" -c 'datasets --version' | awk '{print $3}')
  emit "NCBI datasets" "${v:-unknown}" "$tools_img" "v2 (rolling)" \
       "ROLLING: the Dockerfile fetches NCBI's /command-line/v2/ path, which serves the current v2 build. The version column is what this image received."
  v=$(run "$tools_img" -c 'cat /etc/debian_version')
  emit "Debian (tools base)" "${v:-unknown}" "$tools_img" "bookworm-20241111-slim" "immutable base tag"
  emit "tools image ID" "$(inspect "$tools_img")" "$tools_img" n/a \
       "local build; no registry digest exists (never pushed)"

  # ---- analysis image ---------------------------------------------------
  v=$(run "$anal_img" -c 'python3 -V' | awk '{print $2}')
  emit Python "${v:-unknown}" "$anal_img" "3.12-slim-bookworm (rolling patch)" \
       "ROLLING: the base tag pins the minor series, not the patch. The version column is what this image received."
  while read -r pkg; do
      n=${pkg%%==*}; ver=${pkg##*==}
      [[ -z "$n" ]] && continue
      emit "$n" "$ver" "$anal_img" "$ver" "pinned in docker/analysis/requirements.txt"
  done < <(run "$anal_img" -c 'pip list --format=freeze' \
           | grep -iE '^(numpy|pandas|matplotlib|openpyxl|pysam|scipy)==' | sort)
  emit "analysis image ID" "$(inspect "$anal_img")" "$anal_img" n/a \
       "local build; no registry digest exists (never pushed)"

  # ---- breseq (upstream image, so a real digest does exist) -------------
  v=$(docker run --rm --platform linux/amd64 --entrypoint sh "$breseq_img" \
        -c 'breseq --version' 2>/dev/null | awk '/^breseq/{print $2; exit}')
  emit breseq "${v:-unknown}" "$breseq_img" "0.40.1--h3be2455_0" \
       "upstream biocontainer, immutable tag; optional path, not used for the reported results"
  d=$(docker image inspect --format '{{index .RepoDigests 0}}' "$breseq_img" 2>/dev/null || true)
  emit "breseq image digest" "${d##*@}" "$breseq_img" n/a "registry digest (this image is pulled, not built)"
} > "${out:-/dev/stdout}"

[[ -n "$out" ]] && echo "[software-versions] -> $out" >&2 || true
