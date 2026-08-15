#!/usr/bin/env bash
#
# Assert bin/fetch_ena_reads.sh over saved ENA filereports.
#
# Every branch is exercised without a network: --report substitutes a saved
# portal response, and the "download" is a file:// URL curl reads off disk. The
# point is not that curl works, it is that the *decisions* are right --- which
# copy is chosen, when a checksum is enforced, and whether a bad file survives.
#
# The two checks that earn their keep:
#
#   md5 mismatch     the output must be DELETED, not left on disk. A pipeline
#                    that leaves a corrupt FASTQ behind hands storeDir a cached
#                    "result" that never re-downloads.
#   archive fallback the analysed-file checksum must be SKIPPED with a warning,
#                    not failed. SRA regenerates fastq_ftp from its own archive:
#                    same sequence, different bytes. Failing there would make a
#                    correct download look like corruption.
#
# Runs anywhere with bash, awk, curl and an md5 tool. No image, no network.

set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FETCH="${HERE}/../bin/fetch_ena_reads.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

fail=0
check () {  # check <description> <condition-exit-status>
    if [ "$2" -eq 0 ]; then echo "  ok   $1"; else echo "  FAIL $1"; fail=$((fail + 1)); fi
}

# ---- fixtures -------------------------------------------------------------
printf '@r1\nACGT\n+\nIIII\n' > "${TMP}/reads.fastq.gz"   # content, not real gzip
GOOD_MD5="$(md5sum "${TMP}/reads.fastq.gz" 2>/dev/null | cut -d' ' -f1 || md5 -q "${TMP}/reads.fastq.gz")"
FILE_URL="file://${TMP}/reads.fastq.gz"

printf 'sample_id\tfilename\tbytes\tmd5\treads\tbases\n' > "${TMP}/md5table.tsv"
printf 'demo\tdemo.fastq.gz\t20\t%s\t1\t4\n' "$GOOD_MD5" >> "${TMP}/md5table.tsv"

hdr='run_accession\tsubmitted_ftp\tsubmitted_md5\tfastq_ftp\tfastq_md5\n'

printf "$hdr"                                     >  "${TMP}/submitted.tsv"
printf 'ERR1\t%s\t%s\t\t\n' "$FILE_URL" "$GOOD_MD5" >> "${TMP}/submitted.tsv"

printf "$hdr"                                     >  "${TMP}/archive.tsv"
printf 'ERR2\t\t\t%s\t%s\n' "$FILE_URL" "$GOOD_MD5" >> "${TMP}/archive.tsv"

printf "$hdr"                                     >  "${TMP}/badmd5.tsv"
printf 'ERR3\t%s\t%s\t\t\n' "$FILE_URL" "00000000000000000000000000000000" >> "${TMP}/badmd5.tsv"

printf "$hdr"                                     >  "${TMP}/empty.tsv"

printf "$hdr"                                     >  "${TMP}/two.tsv"
printf 'ERR4\t%s\t%s\t\t\n' "$FILE_URL" "$GOOD_MD5" >> "${TMP}/two.tsv"
printf 'ERR5\t%s\t%s\t\t\n' "$FILE_URL" "$GOOD_MD5" >> "${TMP}/two.tsv"

printf "$hdr"                                     >  "${TMP}/nofiles.tsv"
printf 'ERR6\t\t\t\t\n'                           >> "${TMP}/nofiles.tsv"

run () {  # run <report> <sample-id> [extra args...]
    local report="$1" sid="$2"; shift 2
    rm -f "${TMP}/out.gz"
    bash "$FETCH" --accession SAMN00000001 --sample-id "$sid" \
        --out "${TMP}/out.gz" --report "${TMP}/${report}" \
        --md5-table "${TMP}/md5table.tsv" "$@" > "${TMP}/log" 2>&1
}

# ---- 1. the submitted copy, verified against the analysed file -------------
run submitted.tsv demo; rc=$?
check "submitted copy: exits 0"   "$rc"
[ -s "${TMP}/out.gz" ]; check "submitted copy: file kept" $?
grep -q "(submitted)"                     "${TMP}/log"; check "submitted copy: reports its source"  $?
grep -q "the analysed file"               "${TMP}/log"; check "submitted copy: checks deposited_files.tsv" $?

# ---- 2. the archive's regenerated copy -------------------------------------
run archive.tsv demo; rc=$?
check "archive copy: exits 0"     "$rc"
grep -q "no submitted file"               "${TMP}/log"; check "archive copy: warns about the fallback"    $?
grep -q "cannot confirm this is the analysed file" "${TMP}/log"
check "archive copy: skips the analysed-file check rather than failing it" $?

# ---- 3. a corrupt download must not survive --------------------------------
run badmd5.tsv demo; rc=$?
[ "$rc" -ne 0 ];              check "md5 mismatch: exits non-zero" $?
[ ! -e "${TMP}/out.gz" ];    check "md5 mismatch: output deleted" $?

# ---- 4. nothing released yet -----------------------------------------------
run empty.tsv demo; rc=$?
[ "$rc" -ne 0 ];              check "no runs: exits non-zero" $?
grep -q "not yet released"                "${TMP}/log"; check "no runs: names the likely cause"          $?

# ---- 5. an accession that is not one sample --------------------------------
run two.tsv demo; rc=$?
[ "$rc" -ne 0 ];              check "two runs: exits non-zero" $?
grep -q "sra_accession"                   "${TMP}/log"; check "two runs: says how to disambiguate"       $?

# ---- 6. a released run whose files are not there yet -----------------------
run nofiles.tsv demo; rc=$?
[ "$rc" -ne 0 ];              check "no files: exits non-zero" $?

# ---- 7. a sample this repository never deposited ---------------------------
run submitted.tsv someone_elses; rc=$?
check "unknown sample: exits 0"   "$rc"
grep -q "skipping the analysed-file check" "${TMP}/log"
check "unknown sample: skips the check quietly"    $?

# ---- 8. a BioProject cannot identify one sample ----------------------------
bash "$FETCH" --accession PRJNA1513130 --out "${TMP}/out.gz" > "${TMP}/log" 2>&1; rc=$?
[ "$rc" -ne 0 ];              check "BioProject accession: rejected" $?
grep -q "cannot identify one"             "${TMP}/log"; check "BioProject accession: explains why"        $?

# ---- 9. the shipped accessions and checksums agree -------------------------
# A stale accession or a missing checksum row cannot be caught at run time: the
# fetch would succeed and silently deliver the wrong sample, or skip the check
# that proves the file is the analysed one. Both tables are small and hand-
# maintained, so assert their agreement here.
SHEET="${HERE}/../assets/samplesheets/all.csv"
DEP="${HERE}/../assets/deposited_files.tsv"

sheet_ids="$(awk -F, '!/^#/ && $1 != "sample_id" && NF > 1 {print $1}' "$SHEET" | sort)"
dep_ids="$(awk -F'\t' '!/^#/ && $1 != "sample_id" && NF > 1 {print $1}' "$DEP" | sort)"
[ "$sheet_ids" = "$dep_ids" ]
check "all.csv and deposited_files.tsv cover the same samples" $?

n_sheet="$(printf '%s\n' "$sheet_ids" | grep -c .)"
[ "$n_sheet" -eq 7 ]; check "all seven replicates are listed" $?

bios="$(awk -F, '!/^#/ && $1 != "sample_id" && NF > 1 {print $6}' "$SHEET")"
[ "$(printf '%s\n' "$bios" | grep -c '^SAMN[0-9][0-9]*$')" -eq 7 ]
check "every row carries a well-formed BioSample accession" $?
[ "$(printf '%s\n' "$bios" | sort -u | grep -c .)" -eq 7 ]
check "no two samples share a BioSample accession" $?

# Samplesheet columns must be read BY NAME, never by position. Inserting
# `biosample_accession` shifted reference_set from column 6 to 7 and silently
# fed an accession to bin/assigned_depth.sh as a directory name -- no error, just
# a lookup for a reference set that cannot exist. The helper below is what makes
# that impossible; assert it is still in use.
grep -q 'sheet_col ()' "${HERE}/../bin/assigned_depth.sh"
check "assigned_depth.sh looks samplesheet columns up by name" $?

refsets="$(awk -F, '!/^#/ && $1 != "sample_id" && NF > 1 {
             for (i = 1; i <= NF; i++) if (i == rs) print $i }
           $1 == "sample_id" { for (i = 1; i <= NF; i++) if ($i == "reference_set") rs = i }' "$SHEET")"
missing=0
for r in $refsets; do [ -f "${HERE}/../$r" ] || missing=$((missing + 1)); done
[ "$missing" -eq 0 ] && [ "$(printf '%s\n' "$refsets" | grep -c .)" -eq 7 ]
check "every reference_set named in all.csv resolves to a real file" $?

md5s="$(awk -F'\t' '!/^#/ && $1 != "sample_id" && NF > 1 {print $4}' "$DEP")"
[ "$(printf '%s\n' "$md5s" | grep -c '^[0-9a-f]\{32\}$')" -eq 7 ]
check "every deposited file has a well-formed md5" $?
[ "$(printf '%s\n' "$md5s" | sort -u | grep -c .)" -eq 7 ]
check "no two deposited files share an md5" $?

if [ "$fail" -gt 0 ]; then
    echo; echo "FAILED: ${fail} check(s)"; exit 1
fi
echo; echo "sra fetch: all checks passed"
