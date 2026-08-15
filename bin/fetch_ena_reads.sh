#!/usr/bin/env bash
#
# Fetch one sample's reads from the public archive, by accession.
#
# Resolution goes through ENA's portal API rather than the SRA toolkit for two
# reasons. First, the toolkit publishes no linux/arm64 build, and this study's
# reference platform is Apple Silicon; adding it would mean a second emulated
# image for a download. Second, ENA serves the *submitted* file -- the exact
# bytes we uploaded -- alongside the archive's regenerated copy, and only the
# submitted file can be checked against the md5 in assets/deposited_files.tsv.
# curl is already in the tools image, so this adds no dependency at all.
#
# The accession may be a run (SRR/ERR/DRR) or a sample (SAMN/SAMEA/ERS/SRS).
# A sample accession is what this repository ships, because run accessions do
# not exist until the submission is processed; ENA resolves the sample to its
# run at fetch time, so nothing here has to be edited when they are issued.
#
# Usage:
#   fetch_ena_reads.sh --accession SAMN62407365 --sample-id lowinput_s1_r1 \
#                      --out lowinput_s1_r1.fastq.gz \
#                      [--md5-table assets/deposited_files.tsv] \
#                      [--portal https://www.ebi.ac.uk/ena/portal/api] \
#                      [--report saved-filereport.tsv]
#
# --report substitutes a saved filereport for the live query, which is how
# tests/sra_fetch.sh exercises every branch below without a network.
#
# Exit status is non-zero, with a diagnosis on stderr, for every failure mode
# we can distinguish -- an unreleased submission and a typo'd accession look
# identical in the API's response (an empty table), so the message covers both.

set -euo pipefail

PORTAL="https://www.ebi.ac.uk/ena/portal/api"
ACCESSION=""
SAMPLE_ID=""
OUT=""
MD5_TABLE=""
SAVED_REPORT=""

while [ $# -gt 0 ]; do
    case "$1" in
        --accession) ACCESSION="$2"; shift 2 ;;
        --sample-id) SAMPLE_ID="$2"; shift 2 ;;
        --out)       OUT="$2";       shift 2 ;;
        --md5-table) MD5_TABLE="$2"; shift 2 ;;
        --portal)    PORTAL="$2";    shift 2 ;;
        --report)    SAVED_REPORT="$2"; shift 2 ;;
        -h|--help)   sed -n '2,31p' "$0"; exit 0 ;;
        *) echo "fetch_ena_reads.sh: unknown argument '$1'" >&2; exit 2 ;;
    esac
done

[ -n "$ACCESSION" ] || { echo "fetch_ena_reads.sh: --accession is required" >&2; exit 2; }
[ -n "$OUT" ]       || { echo "fetch_ena_reads.sh: --out is required"       >&2; exit 2; }
SAMPLE_ID="${SAMPLE_ID:-$ACCESSION}"

case "$ACCESSION" in
    SRR*|ERR*|DRR*|SAMN*|SAMEA*|SAMD*|ERS*|SRS*|DRS*|SRX*|ERX*|DRX*) ;;
    *) echo "fetch_ena_reads.sh: '$ACCESSION' is not a run, experiment or sample accession." >&2
       echo "  A BioProject (PRJNA...) covers every sample at once and cannot identify one." >&2
       exit 2 ;;
esac

# ---------------------------------------------------------------------------
# 1. Resolve the accession to a file listing
# ---------------------------------------------------------------------------
FIELDS="run_accession,submitted_ftp,submitted_md5,fastq_ftp,fastq_md5"
REPORT="$(mktemp)"
trap 'rm -f "$REPORT"' EXIT

if [ -n "$SAVED_REPORT" ]; then
    echo "==> using saved filereport ${SAVED_REPORT} for ${ACCESSION} (${SAMPLE_ID})" >&2
    cp "$SAVED_REPORT" "$REPORT"
else
    echo "==> resolving ${ACCESSION} (${SAMPLE_ID}) at ${PORTAL}" >&2
    curl -fsSL --retry 10 --retry-delay 5 --retry-all-errors \
         -o "$REPORT" \
         "${PORTAL}/filereport?accession=${ACCESSION}&result=read_run&fields=${FIELDS}&format=tsv" \
      || { echo "error: ENA portal query failed for ${ACCESSION}" >&2; exit 1; }
fi

# Row 1 is the header. Anything less means the accession resolved to nothing.
# awk rather than `wc -l`, which undercounts a file with no trailing newline.
# The parentheses are load-bearing: `print NR > 0` is a redirect to a file
# named "0", not a comparison.
n_rows="$(awk 'END {print (NR > 0 ? NR - 1 : 0)}' "$REPORT")"
if [ "$n_rows" -lt 1 ]; then
    cat >&2 <<EOF
error: ENA returned no runs for ${ACCESSION} (${SAMPLE_ID}).

  The two causes are indistinguishable from here:
    - the submission is not yet released, so no public record exists; or
    - the accession is wrong.

  Reads for this study are deposited under BioProject PRJNA1513130 and become
  public on release. Until then, run without --fetch_from_sra and point the
  samplesheet at local FASTQs.
EOF
    exit 1
fi
if [ "$n_rows" -gt 1 ]; then
    echo "error: ${ACCESSION} resolves to ${n_rows} runs; this pipeline expects one per sample." >&2
    echo "  Put the specific run accession in the samplesheet's sra_accession column:" >&2
    tail -n +2 "$REPORT" | cut -f1 | sed 's/^/    /' >&2
    exit 1
fi

# Field-by-field with awk, NOT `IFS=$'\t' read`. Tab is IFS whitespace, so read
# collapses runs of tabs into one delimiter: a run with no submitted file has an
# empty column 2, and every later field would silently shift left by one.
field () { tail -n +2 "$REPORT" | head -1 | awk -F'\t' -v i="$1" '{print $i}'; }
RUN="$(field 1)"
SUB_FTP="$(field 2)"
SUB_MD5="$(field 3)"
FQ_FTP="$(field 4)"
FQ_MD5="$(field 5)"

# ---------------------------------------------------------------------------
# 2. Choose which copy to download
# ---------------------------------------------------------------------------
# `submitted_ftp` is the file as uploaded and is byte-identical to what was
# analysed. `fastq_ftp` is regenerated by the archive: same sequence, different
# bytes, and read names may be rewritten. Prefer the former and say so, because
# a silent fallback would turn a failed md5 check into a mystery.
if [ -n "${SUB_FTP:-}" ]; then
    URLS="$SUB_FTP"; MD5S="$SUB_MD5"; SOURCE="submitted"
elif [ -n "${FQ_FTP:-}" ]; then
    URLS="$FQ_FTP";  MD5S="$FQ_MD5";  SOURCE="archive-generated"
    echo "warning: ${RUN} has no submitted file; falling back to the archive's" >&2
    echo "  regenerated FASTQ. Sequence is the same; bytes and read names may not be." >&2
else
    echo "error: ${RUN} lists no downloadable file. The run exists but its data are" >&2
    echo "  not yet available -- most likely still being processed after release." >&2
    exit 1
fi

case "$URLS" in
    *";"*) echo "error: ${RUN} has multiple files (${URLS}). These are single-end" >&2
           echo "  nanopore runs; a paired or split layout means the wrong record." >&2
           exit 1 ;;
esac

# ENA reports these as bare host/path ("ftp.sra.ebi.ac.uk/vol1/..."), which is
# served over https. Match ANY scheme, not just http: a file:// URL is how
# tests/sra_fetch.sh drives this offline, and prefixing that with https:// sends
# curl into its full retry schedule against a hostname that cannot exist.
URL="$URLS"
case "$URL" in *://*) ;; *) URL="https://${URL}" ;; esac

# ---------------------------------------------------------------------------
# 3. Download and verify
# ---------------------------------------------------------------------------
echo "==> ${SAMPLE_ID}: ${RUN} (${SOURCE}) <- ${URL}" >&2
# --no-progress-meter, not -s: a 6 GB transfer's progress bar is unreadable
# noise in .command.err, but errors must still be reported.
curl -fL -C - --no-progress-meter --retry 20 --retry-delay 5 --retry-all-errors \
     -o "$OUT" "$URL"

# md5sum on Linux (and so in the tools image); `md5 -q` on macOS, where this is
# run by hand and by tests/sra_fetch.sh.
md5_of () {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1
    else md5 -q "$1"; fi
}

verify_md5 () {
    local want="$1" label="$2"
    local got
    got="$(md5_of "$OUT")"
    if [ "$got" != "$want" ]; then
        echo "error: md5 mismatch against ${label} for ${SAMPLE_ID}" >&2
        echo "  expected ${want}" >&2
        echo "  got      ${got}" >&2
        rm -f "$OUT"
        return 1
    fi
    echo "    md5 ok against ${label}: ${got}" >&2
}

if [ -n "${MD5S:-}" ]; then
    verify_md5 "$MD5S" "the archive's own checksum"
else
    echo "warning: ENA reported no checksum for ${RUN}; transport integrity unverified." >&2
fi

# The stronger check: does this file hash to what we analysed? Only meaningful
# for the submitted copy, and only for samples this repository deposited.
if [ -n "$MD5_TABLE" ] && [ -f "$MD5_TABLE" ]; then
    want="$(awk -F'\t' -v s="$SAMPLE_ID" '$1 == s {print $4}' "$MD5_TABLE")"
    if [ -z "$want" ]; then
        echo "    ${SAMPLE_ID} is not in ${MD5_TABLE}; skipping the analysed-file check." >&2
    elif [ "$SOURCE" != "submitted" ]; then
        echo "warning: cannot confirm this is the analysed file. ${MD5_TABLE} records the" >&2
        echo "  submitted bytes, and the archive served its regenerated copy instead." >&2
    else
        verify_md5 "$want" "assets/deposited_files.tsv (the analysed file)"
    fi
fi

echo "==> ${SAMPLE_ID}: $(stat -c %s "$OUT" 2>/dev/null || stat -f %z "$OUT") bytes" >&2
