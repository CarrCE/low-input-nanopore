#!/usr/bin/env bash
#
# How far is the carrier-derived contaminant from the reference we assign it
# against? Runs breseq over the contaminant reads of one or more finished
# replicates and reports mapping rate, depth and variant counts per replicate.
#
# Why this exists
# ---------------
# The original lowinput_s1 analysis did not assign reads to stock E. coli K-12
# MG1655; it built a breseq consensus of the E. coli actually present in the
# carrier prep and subtracted reads matching that. This pipeline assigns against
# the stock reference instead, so it owes the reader evidence that the stock
# reference is close enough for that to be sound. This script is that evidence.
#
# It is a validation check, not part of the analysis path. Nothing in results/
# depends on it. The equivalent in-pipeline route is --breseq_consensus, which
# re-maps; this is the cheap route for a run that has already finished.
#
# Why it does not re-map
# ----------------------
# A finished run already recorded, per read, which organism won the competitive
# assignment. Pulling the contaminant read IDs out of <sample>.assignments.tsv.gz
# takes seconds, and `seqkit grep` extracts them from the FASTQ in seconds more.
# Re-mapping a 16 GB FASTQ to answer this would take hours and would produce the
# same seed set. Do not re-run the pipeline for this.
#
# Note on comparability across datasets: the seed set is whatever competitive
# assignment awarded to the contaminant, so its character depends on the
# community. Where a community member is itself an E. coli (lowinput_s1), reads
# that could be either strain land in the ambiguous class and never reach
# breseq, making the seed purer. Where nothing competes (lowinput_s2), every
# enterobacterial read is assigned to the contaminant and the seed is broader.
# Variant counts are therefore comparable within a dataset, and only loosely
# between datasets. See docs/TODO.md item 2.
#
# Usage:
#     bash bin/contaminant_divergence.sh                 # every replicate in results/
#     bash bin/contaminant_divergence.sh lowinput_s1_r1  # named replicates only
#
#     OUTDIR=... RESULTS=... THREADS=8 bash bin/contaminant_divergence.sh
#
# Idempotent: a replicate whose output.gd already exists is skipped.
#
# Requirements: bash, docker, and the images from `make images`.

set -euo pipefail

REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
RESULTS="${RESULTS:-${REPO}/results}"
OUTDIR="${OUTDIR:-${REPO}/results/contaminant_divergence}"
THREADS="${THREADS:-8}"

TOOLS_IMAGE="${TOOLS_IMAGE:-low-input-nanopore/tools:0.1.0}"
ANALYSIS_IMAGE="${ANALYSIS_IMAGE:-low-input-nanopore/analysis:0.1.0}"
# Pinned in conf/base.config; amd64-only and must run as root, hence the
# --platform and HOME overrides here.
BRESEQ_IMAGE="${BRESEQ_IMAGE:-quay.io/biocontainers/breseq:0.40.1--h3be2455_0}"

log() { printf '[contaminant-divergence] %s\n' "$*"; }
die() { printf '[contaminant-divergence] error: %s\n' "$*" >&2; exit 1; }

command -v docker >/dev/null 2>&1 || die "'docker' is required but not on PATH"

if [ "$#" -gt 0 ]; then
    REPLICATES=("$@")
else
    REPLICATES=()
    for d in "${RESULTS}"/*/; do
        r="$(basename "${d}")"
        [ -f "${d}/competitive/${r}.assignments.tsv.gz" ] && REPLICATES+=("${r}")
    done
fi
[ "${#REPLICATES[@]}" -gt 0 ] || die "no finished replicates found under ${RESULTS}"

mkdir -p "${OUTDIR}"
log "replicates: ${REPLICATES[*]}"

for rep in "${REPLICATES[@]}"; do
    # The reference set comes from the samplesheet rather than from the sample
    # name. They usually agree, but not always -- the smoke-test sample
    # `test_s2` runs against the lowinput_s2 reference set -- and guessing from
    # the name silently looks up the wrong reference or none at all.
    set_name="$(awk -F, -v r="${rep}" '
        $1 == r { n = split($6, p, "/"); sub(/\.tsv$/, "", p[n]); print p[n]; exit }
    ' "${REPO}"/assets/samplesheets/*.csv 2>/dev/null | head -1)"
    set_name="${set_name:-${rep%_r*}}"

    out="${OUTDIR}/${rep}"
    assignments="${RESULTS}/${rep}/competitive/${rep}.assignments.tsv.gz"
    contig_map="${RESULTS}/references/${set_name}/contig_map.tsv"
    combined="${RESULTS}/references/${set_name}/combined.fasta"

    if [ -s "${out}/out/output/output.gd" ]; then
        log "${rep}: already done"
        continue
    fi
    [ -f "${assignments}" ] || die "${rep}: no assignments at ${assignments}"
    [ -f "${contig_map}" ]  || die "${rep}: no contig map at ${contig_map}"

    mkdir -p "${out}"
    log "${rep}: extracting contaminant reads"

    # The contaminant organism name comes from the contig map rather than being
    # hardcoded, so this works for any reference set that declares one.
    docker run --rm -v "${REPO}":/repo -v "${OUTDIR}":/out -w /repo "${ANALYSIS_IMAGE}" \
        python3 - "${assignments#${REPO}/}" "${contig_map#${REPO}/}" "${rep}" <<'PY'
import gzip, sys
assignments, contig_map, rep = sys.argv[1:4]
orgs = {f[1] for f in (l.rstrip("\n").split("\t") for l in open(contig_map))
        if len(f) >= 3 and f[2] == "contaminant"}
if not orgs:
    sys.exit(f"error: {contig_map} declares no contaminant organism")
n = 0
with gzip.open(assignments, "rt") as fh, open(f"/out/{rep}/ids.txt", "w") as out:
    next(fh)
    for line in fh:
        f = line.split("\t", 3)
        if f[1] in orgs and f[2] == "assigned":
            out.write(f[0] + "\n")
            n += 1
print(f"  {n} contaminant read IDs ({', '.join(sorted(orgs))})")
if n == 0:
    sys.exit("error: no reads were assigned to the contaminant")
PY

    log "${rep}: pulling reads and contaminant reference"
    docker run --rm -v "${REPO}":/repo -v "${OUTDIR}":/out -w /repo "${TOOLS_IMAGE}" bash -c "
        set -e
        # Column 4 is fastq; column 2 is experiment. Reading \$2 here made the
        # test below compare a non-empty experiment name against a file that
        # cannot exist, so the data/<rep>.fastq fallback never fired and this
        # script could not find any FASTQ at all.
        fq=\$(awk -F, -v r='${rep}' '\$1 == r {print \$4}' /repo/assets/samplesheets/*.csv | head -1)
        [ -n \"\$fq\" ] || { echo 'error: ${rep} not found in any samplesheet' >&2; exit 1; }
        [ -f \"\$fq\" ] || { echo \"error: FASTQ not found for ${rep}: \$fq\" >&2; exit 1; }
        seqkit grep -f /out/${rep}/ids.txt \"\$fq\" > /out/${rep}/contam.fastq 2>/dev/null
        awk -F'\t' 'NR>1 && \$3==\"contaminant\" {print \$1}' \
            '${contig_map#${REPO}/}' > /out/${rep}/ctgs.txt
        samtools faidx '${combined#${REPO}/}'
        xargs samtools faidx '${combined#${REPO}/}' \
            < /out/${rep}/ctgs.txt > /out/${rep}/contaminant.fasta
    "

    log "${rep}: running breseq (nanopore mode)"
    docker run --rm --platform linux/amd64 -e HOME=/tmp -v "${OUTDIR}/${rep}":/w -w /w "${BRESEQ_IMAGE}" bash -c "
        breseq -x --no-junction-prediction -j ${THREADS} -n ${rep} \
            -r contaminant.fasta -o out contam.fastq > breseq.log 2>&1
        gdtools APPLY -r contaminant.fasta -f FASTA -o consensus.fasta \
            out/output/output.gd >> breseq.log 2>&1
    " || { tail -30 "${OUTDIR}/${rep}/breseq.log" >&2; die "${rep}: breseq failed"; }

    log "${rep}: done"
done

# ---- summary ---------------------------------------------------------------
summary="${OUTDIR}/summary.tsv"
{
    printf 'replicate\tcontaminant_reads\tmapped_pct\tdepth\tsnps\tindels\tstructural\n'
    for rep in "${REPLICATES[@]}"; do
        gd="${OUTDIR}/${rep}/out/output/output.gd"
        [ -s "${gd}" ] || continue
        ref_bp=$(grep -v '^>' "${OUTDIR}/${rep}/contaminant.fasta" | tr -d '\n' | wc -c | tr -d ' ')
        awk -v rep="${rep}" -v g="${ref_bp}" '
            /^#=INPUT-READS/     {ir=$2}
            /^#=CONVERTED-READS/ {cr=$2}
            /^#=MAPPED-READS/    {mr=$2}
            /^#=MAPPED-BASES/    {mb=$2}
            $1=="SNP"                     {s++}
            $1 ~ /^(SUB|DEL|INS)$/        {i++}
            $1 ~ /^(MOB|AMP|CON|INV)$/    {v++}
            END {printf "%s\t%d\t%.1f\t%.1f\t%d\t%d\t%d\n", rep, ir, 100*mr/cr, mb/g, s+0, i+0, v+0}
        ' "${gd}"
    done
} > "${summary}"

log "summary -> ${summary}"
column -t -s "$(printf '\t')" "${summary}"
