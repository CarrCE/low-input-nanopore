#!/usr/bin/env bash
#
# Assert bin/filter_by_qscore.awk's quality decisions on hand-built reads.
#
# The filter takes a read's mean quality from dorado's qs:f: header tag when it
# is there and recomputes it from the quality string when it is not -- as for
# every read in the archive's regenerated FASTQ. These cases pin both paths:
#
#   tag present          the tag decides, whatever the quality string says
#   no tag, TRIM=0       every base counts (fastcat / wf-metagenomics)
#   no tag, TRIM=60      the first 60 bases are left out, as dorado does
#   short read           a read no longer than TRIM uses every base
#   exactly Q10          with TRIM>0 the value is rounded as dorado writes it
#                        (4 decimals), so a read of Q10 bases passes MINQ=10
#                        exactly as its qs:f:10.000 tag would
#
# Runs anywhere with bash and awk.

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
FILTER="${HERE}/../bin/filter_by_qscore.awk"
pass=0; fail=0
q() { printf "%${2}s" "" | tr ' ' "$1"; }          # q CHAR N: N copies of CHAR
rec() { printf '@%s%s\n%s\n+\n%s\n' "$1" "${2:+	$2}" "$(q A ${#3})" "$3"; }

check() {  # name expected(kept|dropped) trim fastq-record
    local name="$1" want="$2" trim="$3" fq="$4" got
    if [ -n "$trim" ]; then
        got=$(printf '%s' "$fq" | awk -v MINQ=10 -v TRIM="$trim" -f "$FILTER" 2>/dev/null | wc -l | tr -d ' ')
    else
        got=$(printf '%s' "$fq" | awk -v MINQ=10 -f "$FILTER" 2>/dev/null | wc -l | tr -d ' ')
    fi
    [ "$got" -gt 0 ] && got=kept || got=dropped
    if [ "$got" = "$want" ]; then echo "  ok   $name"; pass=$((pass + 1))
    else echo "  FAIL $name: $got, expected $want"; fail=$((fail + 1)); fi
}

# '#' is Q2, '5' is Q20, '+' is Q10.
low60_high40="$(q '#' 60)$(q 5 40)"

check "tag present: tag says Q25, bases say Q2"  kept    ""  "$(rec r1 'qs:f:25.0000' "$(q '#' 100)")"
check "tag present: tag says Q5, bases say Q20"  dropped 60  "$(rec r2 'qs:f:5.0000'  "$(q 5 100)")"
check "no tag, TRIM unset: 60 x Q2 + 40 x Q20"   dropped ""  "$(rec r3 '' "$low60_high40")"
check "no tag, TRIM=0: same read"                dropped 0   "$(rec r4 '' "$low60_high40")"
check "no tag, TRIM=60: same read"               kept    60  "$(rec r5 '' "$low60_high40")"
check "no tag, TRIM=60, 50 bp of Q2"             dropped 60  "$(rec r6 '' "$(q '#' 50)")"
check "no tag, TRIM=60, 50 bp of Q20"            kept    60  "$(rec r7 '' "$(q 5 50)")"
check "no tag, TRIM=60, 60 x Q2 + 2 x Q10"       kept    60  "$(rec r8 '' "$(q '#' 60)$(q + 2)")"
check "tag present, qs:f:10.0000"                kept    60  "$(rec r9 'qs:f:10.0000' "$(q '#' 62)")"

echo "${pass} passed, ${fail} failed"
[ "$fail" -eq 0 ]
