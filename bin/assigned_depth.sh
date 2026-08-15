#!/usr/bin/env bash
#
# Per-base depth of the reads competitive assignment AWARDED to each community
# organism, as opposed to every read that happens to align there.
#
# Why this exists
# ---------------
# COVERAGE_PROFILE runs `samtools depth` over the primary alignments to community
# contigs. It does not filter to the reads that were awarded to the organism, so
# for a member sharing sequence with an abundant relative the depth is largely
# that relative's reads. Pooled over the replicates, E. coli B-1109 shows 158x of
# alignment depth of which 2.08x is attributable; E. faecalis shows 2.09x of
# which 0.0016x is. See docs/TODO.md item 3.
#
# The 1x interpretability threshold already uses attributable depth
# (bin/coverage_attribution.py), which is what decides whether an organism is
# reported as characterised. This script produces the other half: the depth
# PROFILE on that same basis, so the shape of the coverage can be read too.
#
# Why it does not re-run the pipeline
# -----------------------------------
# A finished run recorded, per read, which organism won. The community-assigned
# reads are a small fraction of the whole -- about 712,000 reads across all seven
# replicates, against roughly 60 million sequenced -- so they can be pulled out
# by ID and re-mapped in minutes. Re-running the pipeline to recover the same
# information would take hours and produce the same alignments.
#
# What it does
# ------------
#   1. read IDs and their awarded organism, from <rep>.assignments.tsv.gz
#      (role=sample rows only)
#   2. `seqkit grep` those reads out of the FASTQ
#   3. map with the SAME minimap2 invocation and the SAME combined reference the
#      pipeline used, so alignment scores and therefore positions are comparable
#   4. keep a primary alignment only where the contig's organism IS the organism
#      the read was awarded to -- this is the filter the pipeline lacks
#   5. `samtools depth -a` over that
#
# Step 4 is the whole point. A read awarded to E. coli B-1109 that aligns best to
# lambda is not evidence about B-1109's coverage, and a read awarded to Listeria
# whose primary lands on E. faecalis is not evidence about E. faecalis.
#
# Usage:
#     bash bin/assigned_depth.sh                    # every finished replicate
#     bash bin/assigned_depth.sh lowinput_s1_r1     # named replicates only
#
#     RESULTS=... THREADS=8 bash bin/assigned_depth.sh
#
# Idempotent: a replicate whose output already exists is skipped.
# Requirements: bash, docker, and the images from `make images`.

set -euo pipefail

REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
RESULTS="${RESULTS:-${REPO}/results}"
THREADS="${THREADS:-8}"
TOOLS_IMAGE="${TOOLS_IMAGE:-low-input-nanopore/tools:0.1.0}"
ANALYSIS_IMAGE="${ANALYSIS_IMAGE:-low-input-nanopore/analysis:0.1.0}"

log() { printf '[assigned-depth] %s\n' "$*"; }
die() { printf '[assigned-depth] error: %s\n' "$*" >&2; exit 1; }

command -v docker >/dev/null 2>&1 || die "'docker' is required but not on PATH"

if [ "$#" -gt 0 ]; then
    REPLICATES=("$@")
else
    REPLICATES=()
    for d in "${RESULTS}"/*/; do
        r="$(basename "${d}")"
        [ "${r}" = "test_s2" ] && continue          # smoke test, not an experiment
        [ -f "${d}/competitive/${r}.assignments.tsv.gz" ] && REPLICATES+=("${r}")
    done
fi
[ "${#REPLICATES[@]}" -gt 0 ] || die "no finished replicates under ${RESULTS}"
log "replicates: ${REPLICATES[*]}"

# Look a column up BY NAME. These used to be positional ($6 for reference_set,
# $4 for fastq), which broke the moment a column was inserted before them:
# adding `biosample_accession` shifted reference_set from 6 to 7, and the
# positional read then returned an accession, which `set_name` happily used as a
# directory name. Nothing errored -- it looked for a reference set that could not
# exist. Header-driven lookup cannot fail that way.
sheet_col () {  # sheet_col <sample_id> <column_name>
    awk -F, -v r="$1" -v want="$2" '
        /^#/                { next }
        $1 == "sample_id"   { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $1 == r             { if (want in col) print $(col[want]); exit }
    ' "${REPO}"/assets/samplesheets/*.csv 2>/dev/null | head -1
}

for rep in "${REPLICATES[@]}"; do
    # Reference set from the samplesheet, not from the sample name: they can
    # disagree, and guessing silently looks up the wrong reference.
    set_name="$(sheet_col "${rep}" reference_set)"
    set_name="$(basename "${set_name%.tsv}")"
    set_name="${set_name:-${rep%_r*}}"

    fastq="$(sheet_col "${rep}" fastq)"
    # Under --fetch_from_sra the samplesheet's local path is unused and may be
    # blank; the reads live in the download cache under the sample id. Fall back
    # to it so this step works from a clone that never had local FASTQs.
    if [ -z "${fastq}" ] || [ ! -f "${REPO}/${fastq}" ]; then
        cached="${NXF_SRADIR:-sra}/${rep}.fastq.gz"
        if [ -f "${REPO}/${cached}" ]; then
            log "${rep}: no local FASTQ; using the download cache ${cached}"
            fastq="${cached}"
        fi
    fi
    [ -n "${fastq}" ] || die "${rep}: no fastq in any samplesheet and none in the download cache"

    outdir="${RESULTS}/${rep}/coverage"
    out="${outdir}/${rep}.assigned_depth.tsv.gz"
    assignments="${RESULTS}/${rep}/competitive/${rep}.assignments.tsv.gz"
    contig_map="${RESULTS}/references/${set_name}/contig_map.tsv"
    combined="${RESULTS}/references/${set_name}/combined.fasta"

    if [ -s "${out}" ]; then log "${rep}: already done"; continue; fi
    [ -f "${assignments}" ] || die "${rep}: no assignments at ${assignments}"
    [ -f "${combined}" ]    || die "${rep}: no combined reference at ${combined}"

    mkdir -p "${outdir}"
    work="$(mktemp -d)"
    trap 'rm -rf "${work}"' EXIT

    log "${rep}: collecting community-assigned read IDs"
    docker run --rm -i -v "${REPO}":/repo -v "${work}":/work -w /repo "${ANALYSIS_IMAGE}" \
        python3 - "${assignments#${REPO}/}" <<'PY'
import gzip, sys
path = sys.argv[1]
n = 0
with gzip.open(path, "rt") as fh, open("/work/ids.txt", "w") as ids, \
     open("/work/read_org.tsv", "w") as pairs:
    header = fh.readline().rstrip("\n").split("\t")
    i_id, i_org = header.index("read_id"), header.index("organism")
    i_role, i_call = header.index("role"), header.index("call")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        # Only reads AWARDED to a community organism. Ambiguous calls are
        # excluded deliberately: a tie is a statement that the read cannot be
        # attributed, and counting it toward either organism's coverage would
        # manufacture attribution the assignment step declined to make.
        if f[i_role] == "sample" and f[i_call] == "assigned":
            ids.write(f[i_id] + "\n")
            pairs.write(f"{f[i_id]}\t{f[i_org]}\n")
            n += 1
print(f"[assigned-depth]   {n:,} community-assigned reads")
PY

    log "${rep}: extracting them from $(basename "${fastq}")"
    docker run --rm -v "${REPO}":/repo -v "${work}":/work -w /repo "${TOOLS_IMAGE}" \
        bash -lc "seqkit grep -f /work/ids.txt '${fastq}' > /work/reads.fastq"

    # seqkit grep exits 0 when it matches nothing. If the samplesheet resolves
    # to the wrong replicate, or the reads were re-basecalled since assignment,
    # the profile below is silently shallow and feeds pooled attributable depth,
    # Figure S3 and Table S5 with no other symptom.
    n_ids=$(wc -l < "${work}/ids.txt")
    got=$(( $(wc -l < "${work}/reads.fastq") / 4 ))
    if [ "${got}" -ne "${n_ids}" ]; then
        echo "error: ${rep}: asked seqkit for ${n_ids} reads, recovered ${got}." \
             "The FASTQ and the assignment table disagree; refusing to write a" \
             "depth profile from an incomplete extraction." >&2
        exit 1
    fi

    log "${rep}: mapping against ${set_name} and filtering to the awarded organism"
    docker run --rm -v "${REPO}":/repo -v "${work}":/work -w /repo "${TOOLS_IMAGE}" \
        bash -lc "
        set -euo pipefail
        # Same invocation as MAP_COMPETITIVE, so scores and positions match.
        minimap2 -ax map-ont -t ${THREADS} -N 10 --secondary=yes \
            '${combined#${REPO}/}' /work/reads.fastq 2>/work/minimap2.log \
        | awk -F'\t' -v OFS='\t' '
            NR == FNR                     { org[\$1] = \$2; next }
            FILENAME ~ /contig_map/       { if (FNR > 1 && \$3 == \"sample\")
                                                cmap[\$1] = \$2; next }
            /^@SQ/                        { split(\$2, a, \":\")
                                            if (a[2] in cmap) print
                                            next }
            /^@/                          { print; next }
            {
              # primary only, and only where the contig belongs to the organism
              # this read was actually awarded to
              if (int(\$2/256) % 2 || int(\$2/2048) % 2) next   # 0x100|0x800
              if (!(\$3 in cmap)) next
              if (!(\$1 in org)) next
              if (cmap[\$3] != org[\$1]) next
              print
            }
          ' /work/read_org.tsv '${contig_map#${REPO}/}' - \
        | samtools sort -@ ${THREADS} -m 1G -o /work/sorted.bam -
        samtools index /work/sorted.bam
        samtools depth -a -@ ${THREADS} /work/sorted.bam | gzip -c > /work/depth.tsv.gz
        "

    mv "${work}/depth.tsv.gz" "${out}"
    rm -rf "${work}"; trap - EXIT
    log "${rep}: -> ${out#${REPO}/}"
done

log "done"
