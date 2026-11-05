#!/usr/bin/env bash
#
# The reference-fragmentation test behind `make fragmentation`
# (docs/reference-fragmentation.md): does a contig-level reference change read
# assignment? Three arms on the D6311 replicates:
#
#   A  the published run (contig-level assemblies of the exact yeast strains),
#      taken from results/ and results/nofloor -- run `make all` and
#      `make nofloor` first
#   B  complete genomes of other strains (assets/fragmentation/lowinput_s1_complete.tsv)
#   C  B's genomes cut where A's assemblies break (lowinput_s1_fragmented.tsv)
#
# B and C are each run with the attribution floor on and off. Everything is
# written under results/fragmentation/; the comparison is
# results/fragmentation/summary/fragmentation.json.
#
# The pipeline downloads genomes by accession and cannot download a cut one, so
# C's cut genomes are placed in its reference cache (<refdir>/lowinput_s1_fragmented)
# before C runs; FETCH_GENOMES then finds them and is skipped.
#
# Environment:
#   NF_ARGS   passed to every pipeline run (e.g. --fetch_from_sra --sradir DIR)
#   REFDIR    the pipeline's --refdir (default: <repo>/refs)
#   PROFILE   nextflow profile (default docker)
#
# Requirements: a finished `make all` and `make nofloor`; docker; the images
# from `make images`.

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "${ROOT}"
PROFILE="${PROFILE:-docker}"
REFDIR="${REFDIR:-${ROOT}/refs}"
NF_ARGS="${NF_ARGS:-}"
TOOLS_IMAGE="${TOOLS_IMAGE:-low-input-nanopore/tools:0.1.0}"
ANALYSIS_IMAGE="${ANALYSIS_IMAGE:-low-input-nanopore/analysis:0.1.0}"
OUT="results/fragmentation"

log() { printf '[fragmentation] %s\n' "$*"; }
die() { printf '[fragmentation] error: %s\n' "$*" >&2; exit 1; }

[ -f results/summary/per_organism.tsv ] || die "results/ has no finished run. Run 'make all' first."
[ -f results/nofloor/summary/per_organism.tsv ] || die "results/nofloor not found. Run 'make nofloor' first."
mkdir -p "${OUT}/build" "${REFDIR}"

# Samplesheets: the published D6311 sheet with only the reference set changed.
for arm in complete fragmented; do
    sed "s#assets/references/lowinput_s1.tsv#assets/fragmentation/lowinput_s1_${arm}.tsv#" \
        assets/samplesheets/lowinput_s1.csv > "${OUT}/samplesheet_${arm}.csv"
    grep -q "lowinput_s1_${arm}.tsv" "${OUT}/samplesheet_${arm}.csv" \
        || die "could not point the samplesheet at the ${arm} reference set"
done

run() {  # outdir samplesheet min_aln_frac
    log "pipeline -> $1 (--min_aln_frac $3)"
    # shellcheck disable=SC2086
    "${ROOT}/run.sh" -profile "${PROFILE}" --samplesheet "$2" --mode competitive \
        --min_aln_frac "$3" --refdir "${REFDIR}" --outdir "$1" -resume ${NF_ARGS}
}

# ---- arm B -----------------------------------------------------------------
run "${OUT}/complete" "${OUT}/samplesheet_complete.csv" 0.10

# ---- cut B's yeasts where A's assemblies break -------------------------------
acc() {  # tsv organism
    awk -F'\t' -v o="$2" '!/^#/ && $1 == o {print $2}' "$1"
}
fna() {  # refset accession
    find "${REFDIR}/$1/genomes/ncbi_dataset/data/$2" -name '*.fna' | head -1
}
[ -d "${REFDIR}/lowinput_s1/genomes" ] \
    || die "${REFDIR}/lowinput_s1 not found: arm A's genomes are needed (the cache of 'make all'; pass REFDIR if it ran with --refdir)."
CDIR="${REFDIR}/lowinput_s1_fragmented"
if [ ! -d "${CDIR}/genomes" ]; then
    cp -R "${REFDIR}/lowinput_s1_complete/genomes" "${OUT}/build/genomes.tmp"
    for org in "Saccharomyces cerevisiae" "Cryptococcus neoformans"; do
        a="$(acc assets/references/lowinput_s1.tsv "${org}")"
        b="$(acc assets/fragmentation/lowinput_s1_complete.tsv "${org}")"
        [ -n "${a}" ] && [ -n "${b}" ] || die "no accession for ${org}"
        afna="$(fna lowinput_s1 "${a}")"; bfna="$(fna lowinput_s1_complete "${b}")"
        [ -n "${afna}" ] && [ -n "${bfna}" ] || die "genome FASTA missing for ${org}"
        log "aligning ${org}: ${a} (contigs) to ${b} (complete)"
        docker run --rm -u "$(id -u):$(id -g)" -v "${REFDIR}":/refs -v "${ROOT}":/repo -w /repo \
            "${TOOLS_IMAGE}" minimap2 -x asm20 -c --secondary=no -t 8 \
            "/refs/${bfna#"${REFDIR}"/}" "/refs/${afna#"${REFDIR}"/}" \
            > "${OUT}/build/${a}.on.${b}.paf" 2> "${OUT}/build/${a}.on.${b}.log"
        target="${OUT}/build/genomes.tmp/ncbi_dataset/data/${b}/$(basename "${bfna}")"
        docker run --rm -u "$(id -u):$(id -g)" -v "${REFDIR}":/refs -v "${ROOT}":/repo -w /repo \
            "${ANALYSIS_IMAGE}" python3 bin/fragment_reference.py \
            --complete "/refs/${bfna#"${REFDIR}"/}" --paf "${OUT}/build/${a}.on.${b}.paf" \
            --out "${target}" --cuts "${OUT}/build/cuts_${b}.tsv"
    done
    mkdir -p "${CDIR}"
    mv "${OUT}/build/genomes.tmp" "${CDIR}/genomes"
else
    log "${CDIR} exists; using the cut genomes already there"
fi

# ---- arm C, then both arms with the floor off --------------------------------
run "${OUT}/fragmented"         "${OUT}/samplesheet_fragmented.csv" 0.10
run "${OUT}/complete_nofloor"   "${OUT}/samplesheet_complete.csv"   0
run "${OUT}/fragmented_nofloor" "${OUT}/samplesheet_fragmented.csv" 0

# ---- compare -----------------------------------------------------------------
docker run --rm -u "$(id -u):$(id -g)" -v "${ROOT}":/repo -w /repo "${ANALYSIS_IMAGE}" \
    python3 bin/fragmentation_compare.py \
    --a results --a-nofloor results/nofloor \
    --b "${OUT}/complete" --b-nofloor "${OUT}/complete_nofloor" \
    --c "${OUT}/fragmented" --c-nofloor "${OUT}/fragmented_nofloor" \
    --cuts "${OUT}"/build/cuts_*.tsv --out "${OUT}/summary"
log "done: ${OUT}/summary/fragmentation.json"
