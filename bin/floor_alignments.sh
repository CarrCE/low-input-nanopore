#!/usr/bin/env bash
#
# `make flooralignments`: where the attribution floor's excluded reads align
# (bin/floor_alignments.py; Supplementary Section S2).
#
# Needs a finished `make all` and `make nofloor`. Two steps:
#
#   1. choose the reads and extract their alignments out of each replicate's
#      alignment BAM into <out>/floor_alignments/<sample>.bam -- a few MB from
#      a 7 GB BAM, cached, so a second run costs nothing;
#   2. measure where those alignments land.
#
# Step 1 keeps the evidence behind the floor paragraph after the full BAMs are
# gone: they are ~7 GB per replicate and nothing else reads them, so a finished
# replication can delete them and these extracts still answer the question.
#
# Environment (defaults in brackets):
#   RESULTS  the run with the floor on [results]
#   NOFLOOR  the same run with the floor off [$RESULTS/nofloor]
#   REFDIR   the fetched genomes [refs]
#   OUT      [$RESULTS/summary]
#   CONTROLS community reads >= 10 kb kept per replicate, for comparison [400]
#   SEED     [1]

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "${ROOT}"
die() { printf '[flooralignments] error: %s\n' "$*" >&2; exit 1; }
log() { printf '[flooralignments] %s\n' "$*"; }
abs() { (cd "$1" 2>/dev/null && pwd -P) || die "no such directory: $1"; }

RESULTS="${RESULTS:-results}"
[ -d "${RESULTS}" ] || die "${RESULTS} not found. Run 'make all' first."
RESULTS="$(abs "${RESULTS}")"
NOFLOOR="${NOFLOOR:-${RESULTS}/nofloor}"
[ -d "${NOFLOOR}" ] || die "${NOFLOOR} not found. Run 'make nofloor' first."
NOFLOOR="$(abs "${NOFLOOR}")"
REFDIR="$(abs "${REFDIR:-refs}")"
OUT="${OUT:-${RESULTS}/summary}"
mkdir -p "${OUT}/floor_alignments"
OUT="$(abs "${OUT}")"
ALN="${OUT}/floor_alignments"
CONTROLS="${CONTROLS:-400}"
SEED="${SEED:-1}"
TOOLS_IMAGE="${TOOLS_IMAGE:-low-input-nanopore/tools:0.1.0}"
ANALYSIS_IMAGE="${ANALYSIS_IMAGE:-low-input-nanopore/analysis:0.1.0}"

# Each replicate's reference set, named by its TSV in the samplesheet its run
# used: assets/references/lowinput_s1.tsv -> lowinput_s1, which is also the
# directory the pipeline fetches that set's genomes into.
refset_of() {  # refset_of <sample_id>
    awk -F, -v r="$1" '
        /^#/              { next }
        $1 == "sample_id" { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $1 == r           { print $(col["reference_set"]); exit }
    ' "${ROOT}"/assets/samplesheets/all.csv
}

samples=()
maps=(); sets=()
for a in "${RESULTS}"/*/competitive/*.assignments.tsv.gz; do
    s="$(basename "${a}" .assignments.tsv.gz)"
    rs="$(basename "$(refset_of "${s}")" .tsv)"
    [ -n "${rs}" ] || die "${s}: no reference_set in assets/samplesheets/all.csv"
    cm="${RESULTS}/references/${rs}/contig_map.tsv"
    [ -f "${cm}" ] || die "${s}: no contig map at ${cm}"
    samples+=("${s}")
    maps+=("${s}=/results/references/${rs}/contig_map.tsv")
    sets+=("${s}=${rs}")
done
[ ${#samples[@]} -gt 0 ] || die "no finished replicates in ${RESULTS}"

run_py() {
    docker run --rm -u "$(id -u):$(id -g)" \
        -v "${ROOT}":/repo:ro -v "${RESULTS}":/results:ro -v "${NOFLOOR}":/nofloor:ro \
        -v "${REFDIR}":/refs:ro -v "${OUT}":/out -w /repo "${ANALYSIS_IMAGE}" \
        python3 bin/floor_alignments.py --with /results --without /nofloor \
        --contig-map "${maps[@]}" --refset "${sets[@]}" --refdir /refs \
        --controls "${CONTROLS}" --seed "${SEED}" "$@"
}

# 1. the read ids, then their alignments
if [ ! -s "${ALN}/${samples[0]}.excluded.txt" ]; then
    log "choosing reads -> ${ALN}"
    run_py --ids /out/floor_alignments --write-ids
fi
bams=()
for s in "${samples[@]}"; do
    out="${ALN}/${s}.bam"
    if [ ! -s "${out}" ]; then
        src="${RESULTS}/${s}/alignments/${s}.qname.bam"
        [ -f "${src}" ] || die "${s}: no ${src}. The alignment BAMs are needed once,
  to extract these reads; after that ${ALN}/${s}.bam stands in for them."
        cat "${ALN}/${s}.excluded.txt" "${ALN}/${s}.control.txt" > "${ALN}/${s}.want.txt"
        log "${s}: extracting $(wc -l < "${ALN}/${s}.want.txt" | tr -d ' ') reads"
        docker run --rm -u "$(id -u):$(id -g)" \
            -v "${RESULTS}":/results:ro -v "${OUT}":/out -w /out "${TOOLS_IMAGE}" \
            samtools view -b -N "/out/floor_alignments/${s}.want.txt" \
            -o "/out/floor_alignments/${s}.bam" "/results/${s}/alignments/${s}.qname.bam"
        rm -f "${ALN}/${s}.want.txt"
    fi
    bams+=("${s}=/out/floor_alignments/${s}.bam")
done

# 2. measure
run_py --ids /out/floor_alignments --bam "${bams[@]}" --out /out
log "done: ${OUT}/floor_alignments.json"
