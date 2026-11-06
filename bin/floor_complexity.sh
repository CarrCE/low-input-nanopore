#!/usr/bin/env bash
#
# `make floorcomplexity`: the sequence complexity of the community reads the
# attribution floor excludes, and the simple-repeat content of the reference
# genomes (bin/floor_complexity.py; Supplementary Section S2).
#
# Needs a finished `make all` and `make nofloor`, and the DEPOSITED reads. The
# deposited reads are human-masked. The unmasked originals in data/ are never
# read: analyses use the masked reads only. For each sample the reads are
# looked for, in order, at
#   $READS/<sample>.fastq.gz or $READS/<sample>.masked.fastq.gz, if READS is set
#   $SRADIR/<sample>.fastq.gz          the --fetch_from_sra download cache
#   $RESULTS/<sample>/human/<sample>.masked.fastq.gz
#
# Environment (defaults in brackets):
#   RESULTS  the run with the floor on [results]
#   NOFLOOR  the same run with the floor off [$RESULTS/nofloor]
#   READS, SRADIR [sra], REFDIR [refs]
#   OUT      [$RESULTS/summary]
#
# The reference genomes scanned are those of lowinput_s1 and lowinput_s2, and
# of the reference-fragmentation test's complete set if `make fragmentation`
# has fetched it.

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "${ROOT}"
abs() { (cd "$1" && pwd -P); }
die() { printf '[floorcomplexity] error: %s\n' "$*" >&2; exit 1; }

RESULTS="${RESULTS:-results}"
[ -d "${RESULTS}" ] || die "${RESULTS} not found. Run 'make all' first."
RESULTS="$(abs "${RESULTS}")"
NOFLOOR="${NOFLOOR:-${RESULTS}/nofloor}"
[ -d "${NOFLOOR}" ] || die "${NOFLOOR} not found. Run 'make nofloor' first."
NOFLOOR="$(abs "${NOFLOOR}")"
SRADIR="${SRADIR:-sra}"
REFDIR="$(abs "${REFDIR:-refs}")"
OUT="${OUT:-${RESULTS}/summary}"
mkdir -p "${OUT}"
OUT="$(abs "${OUT}")"
ANALYSIS_IMAGE="${ANALYSIS_IMAGE:-low-input-nanopore/analysis:0.1.0}"

mounts=(-v "${ROOT}":/repo:ro -v "${RESULTS}":/results:ro -v "${NOFLOOR}":/nofloor:ro
        -v "${REFDIR}":/refs:ro -v "${OUT}":/out)
fastq=()
i=0
for a in "${RESULTS}"/*/competitive/*.assignments.tsv.gz; do
    s="$(basename "${a}" .assignments.tsv.gz)"
    f=""
    for c in ${READS:+"${READS}/${s}.fastq.gz" "${READS}/${s}.masked.fastq.gz"} \
             "${SRADIR}/${s}.fastq.gz" "${RESULTS}/${s}/human/${s}.masked.fastq.gz"; do
        [ -f "${c}" ] && { f="${c}"; break; }
    done
    [ -n "${f}" ] || die "${s}: no deposited (masked) reads found; set READS= or SRADIR="
    d="$(abs "$(dirname "${f}")")"
    i=$((i + 1))
    mounts+=(-v "${d}:/reads${i}:ro")
    fastq+=("${s}=/reads${i}/$(basename "${f}")")
    printf '[floorcomplexity] %s: %s\n' "${s}" "${d}/$(basename "${f}")"
done

refsets=(lowinput_s1=assets/references/lowinput_s1.tsv lowinput_s2=assets/references/lowinput_s2.tsv)
[ -d "${REFDIR}/lowinput_s1_complete/genomes" ] \
    && refsets+=(lowinput_s1_complete=assets/fragmentation/lowinput_s1_complete.tsv)

docker run --rm -u "$(id -u):$(id -g)" "${mounts[@]}" -w /repo "${ANALYSIS_IMAGE}" \
    python3 bin/floor_complexity.py --with /results --without /nofloor \
    --fastq "${fastq[@]}" --refset "${refsets[@]}" --refdir /refs --out /out
printf '[floorcomplexity] done: %s/floor_complexity.json\n' "${OUT}"
