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

  Three causes are indistinguishable from here:
    - ENA has not yet mirrored the record, though it is public at NCBI;
    - the submission is not yet released anywhere; or
    - the accession is wrong.

  ENA receives INSDC submissions from NCBI on a delay of a few days. This
  study's runs (SRR40180147-SRR40180153, BioProject PRJNA1513130) were mirrored
  by 2 Nov 2026, so for them an empty result most likely means a mistyped
  accession or a passing portal fault. Check which case you are in:

      curl -s 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=run_accession&format=tsv'

  An empty result with a header row means ENA does not have it yet. Meanwhile,
  download the FASTQs from NCBI and point the samplesheet's fastq column at
  them; note that reads regenerated from the SRA archive differ byte-wise from
  the deposited files and will not match assets/deposited_files.tsv.
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
#
# Retries are this loop, not curl's --retry. On an in-process retry curl
# truncates the file back to where THAT invocation started ("Throwing away N
# bytes"), so a connection ENA dropped at 4 GB restarted the file from zero.
# Each pass here is a fresh curl whose -C - resumes from the file's current
# size. --speed-limit/--speed-time turn a stalled connection into a failure
# the loop can resume, rather than a hang until the task's time limit.
FETCH_ATTEMPTS="${FETCH_ATTEMPTS:-30}"
attempt=1
until curl -fL -C - --no-progress-meter --connect-timeout 30 \
           --speed-limit 1024 --speed-time 120 \
           -o "$OUT" "$URL"; do
    rc=$?
    if [ "$attempt" -ge "$FETCH_ATTEMPTS" ]; then
        echo "error: download of ${URL} failed ${attempt} times (last curl exit ${rc})" >&2
        exit 1
    fi
    # 33: the server refused the byte range, so the partial file is unusable.
    if [ "$rc" -eq 33 ]; then rm -f "$OUT"; fi
    have="$( [ -f "$OUT" ] && { stat -c %s "$OUT" 2>/dev/null || stat -f %z "$OUT"; } || echo 0)"
    echo "    curl exit ${rc}; resuming at byte ${have} (attempt $((attempt + 1)) of ${FETCH_ATTEMPTS})" >&2
    attempt=$((attempt + 1))
    sleep 10
done

# md5sum on Linux (and so in the tools image); `md5 -q` on macOS, where this is
# run by hand and by tests/sra_fetch.sh.
md5_of () {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1
    else md5 -q "$1"; fi
}

size_of () { stat -c %s "$1" 2>/dev/null || stat -f %z "$1"; }

# ENA's download server streams each file from a storage back end. When that
# back end drops mid-file, the server writes the back end's XML error into the
# response AS DATA and carries on from just after it, so the file arrives at
# exactly the right length with ~200 bytes of text where data should be, and
# nothing on the client side reports an error (seen three times on 2-3 Nov
# 2026, e.g. SRR40180153 at byte 758,224,317):
#
#   <Error><Code>ConnectionClosedException</Code><Message>Premature end of
#   Content-Length delimited message body (expected: ...; received: ...)
#   </Message><ErrorMessage/><RequestId/></Error>
#
# Each such blob replaces exactly its own length of real bytes, so fetching that
# byte range again and writing it back restores the file; the md5 check that
# follows is what decides whether the repair worked. Returns 0 if it repaired
# anything. A false positive would need this exact text inside gzip output.
repair_injected () {
    local pat='<Error><Code>[A-Za-z]+</Code><Message>[^<]*</Message><ErrorMessage/><RequestId/></Error>'
    local n=0 off text len part
    while IFS=: read -r off text; do
        [ -n "$off" ] || continue
        len="$(printf '%s' "$text" | LC_ALL=C wc -c | tr -d ' ')"
        part="$(mktemp)"
        if ! curl -fsS --retry 5 --retry-delay 5 --retry-all-errors \
                  -r "${off}-$((off + len - 1))" -o "$part" "$URL" \
           || [ "$(size_of "$part")" -ne "$len" ]; then
            echo "    could not refetch bytes ${off}-$((off + len - 1)) to repair them" >&2
            rm -f "$part"; return 1
        fi
        dd if="$part" of="$OUT" bs=1 seek="$off" conv=notrunc 2>/dev/null
        rm -f "$part"
        echo "    repaired ${len} bytes of server error text written into the download at byte ${off}" >&2
        n=$((n + 1))
    done < <(LC_ALL=C grep -aobE "$pat" "$OUT" || true)
    [ "$n" -gt 0 ]
}

verify_md5 () {
    local want="$1" label="$2" got="${3:-}"
    [ -n "$got" ] || got="$(md5_of "$OUT")"
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
    # Only a file that fails its checksum is scanned, so a good download costs
    # one md5 pass, not a scan as well.
    got="$(md5_of "$OUT")"
    if [ "$got" != "$MD5S" ] && repair_injected; then
        got="$(md5_of "$OUT")"
    fi
    verify_md5 "$MD5S" "the archive's own checksum" "$got"
else
    repair_injected || true
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
