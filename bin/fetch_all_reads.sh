#!/usr/bin/env bash
#
# Download every sample's reads from ENA into the read cache, on the host.
#
# The same downloads happen inside the pipeline with --fetch_from_sra, but there
# each one is a Nextflow task: it inherits the pipeline's time limit and runs in
# a container, and stopping one means stopping the run. Here they are plain
# processes: stop with Ctrl-C and run again, and every partial file resumes
# where it stopped. The pipeline then finds the files in the cache and skips its
# own download step:
#
#   make fetchreads
#   nextflow run . -profile docker --samplesheet assets/samplesheets/all.csv \
#       --fetch_from_sra --mode both
#
# Each file goes through bin/fetch_ena_reads.sh, so it gets the same checks:
# the archive's md5, repair of ENA's injected error text, and the analysed-file
# md5 where that applies. A file is downloaded as <sample>.fastq.gz.part and
# renamed only once it passes, so an interrupted download is never mistaken for
# a finished one.
#
# ENA limits throughput per connection, not per client, so downloading every
# sample at once is several times faster than one at a time.
#
# Usage:  bash bin/fetch_all_reads.sh
# Env:    SRADIR      cache directory (default: sra/, the pipeline's --sradir default)
#         SAMPLESHEET samplesheet to read (default: assets/samplesheets/all.csv)
#         JOBS        downloads at once (default: all samples)
# Needs:  bash, curl, awk, and md5sum or md5. No Docker.

set -euo pipefail

REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
SRADIR="${SRADIR:-${REPO}/sra}"
SHEET="${SAMPLESHEET:-${REPO}/assets/samplesheets/all.csv}"
MD5_TABLE="${REPO}/assets/deposited_files.tsv"
mkdir -p "${SRADIR}/logs"

log() { printf '[fetchreads] %s\n' "$*"; }

# sample_id <TAB> accession, by column NAME (sra_accession wins, as in main.nf).
rows="$(awk -F, '
    /^#/ { next }
    $1 == "sample_id" { for (i = 1; i <= NF; i++) col[$i] = i; next }
    NF > 1 {
        acc = ("sra_accession" in col) ? $(col["sra_accession"]) : ""
        if (acc == "" && ("biosample_accession" in col)) acc = $(col["biosample_accession"])
        print $1 "\t" acc
    }' "$SHEET")"
[ -n "$rows" ] || { echo "error: no samples in ${SHEET}" >&2; exit 1; }
JOBS="${JOBS:-$(printf '%s\n' "$rows" | grep -c .)}"

fetch_one () {  # fetch_one <sample_id> <accession>
    local s="$1" acc="$2" out="${SRADIR}/$1.fastq.gz"
    if [ -s "$out" ]; then log "have   ${s}"; return 0; fi
    if [ -z "$acc" ]; then log "error  ${s}: no accession in ${SHEET}"; return 1; fi
    log "fetch  ${s} (${acc}); log in ${SRADIR}/logs/${s}.log"
    if bash "${REPO}/bin/fetch_ena_reads.sh" --accession "$acc" --sample-id "$s" \
            --out "${out}.part" --md5-table "$MD5_TABLE" > "${SRADIR}/logs/${s}.log" 2>&1; then
        mv "${out}.part" "$out"
        log "done   ${s}"
    else
        log "FAILED ${s}; see ${SRADIR}/logs/${s}.log"
        return 1
    fi
}

pids=(); names=(); failed=0
while IFS=$'\t' read -r s acc; do
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$JOBS" ]; do sleep 2; done
    fetch_one "$s" "$acc" & pids+=($!); names+=("$s")
done <<< "$rows"
for i in "${!pids[@]}"; do wait "${pids[$i]}" || failed=$((failed + 1)); done

n="$(printf '%s\n' "$rows" | grep -c .)"
if [ "$failed" -gt 0 ]; then
    log "${failed} of ${n} samples failed; run again to resume them"
    exit 1
fi
log "all ${n} samples verified in ${SRADIR}"
