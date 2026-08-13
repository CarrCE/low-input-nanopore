#!/usr/bin/env python3
"""Assert the --min-aln-frac attribution floor over real alignments.

The floor stops an organism claiming a read it barely explains: a read is
credited to its winner with its FULL length, so without it a 70 kb read touching
a genome over 200 bp booked 70 kb against that genome. See
bin/attribution_threshold.py for the distribution the default came from.

This runs the real assigner over the smoke-test BAM at three settings and checks
the properties that must hold whatever the data happens to contain:

  floor 0.0   no read is ever called low_coverage -- the floor is genuinely off,
              so a run without it is unchanged from before the feature existed
  floor 0.9   the mechanism actually fires, and every read it fires on really
              is below the floor
  any floor   every read still lands in exactly one class; the floor moves reads
              between classes and never loses one
  floor 0.1   nothing assigned survives below the floor

The 0.9 case matters: a test that only ran the default would pass just as well
against code that never fired at all, since a clean 40,000-read subsample need
not contain a single weakly-attributed read.

Usage (wired into `make check`, which supplies the paths):
  attribution_floor.py --bam results/test_s2/alignments/test_s2.qname.bam \
                       --contig-map results/references/lowinput_s2/contig_map.tsv
"""

import argparse
import gzip
import os
import subprocess
import sys
import tempfile

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")


def run_assigner(bam, contig_map, frac, workdir):
    prefix = os.path.join(workdir, f"f{frac}")
    subprocess.run(
        [sys.executable, os.path.join(BIN, "assign_reads.py"), bam,
         "--contig-map", contig_map, "--prefix", prefix,
         "--min-aln-frac", str(frac), "--mode", "competitive"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    rows = []
    with gzip.open(prefix + ".assignments.tsv.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append((f[ix["read_id"]], f[ix["call"]],
                         int(f[ix["read_length"]]), int(f[ix["aligned_bases"]])))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--contig-map", required=True)
    args = ap.parse_args()

    failures = []

    def check(cond, msg):
        print(("  ok   " if cond else "  FAIL ") + msg)
        if not cond:
            failures.append(msg)

    with tempfile.TemporaryDirectory() as tmp:
        off = run_assigner(args.bam, args.contig_map, 0.0, tmp)
        dflt = run_assigner(args.bam, args.contig_map, 0.10, tmp)
        hard = run_assigner(args.bam, args.contig_map, 0.9, tmp)

    print(f"assignments read: {len(off):,}")

    check(not [r for r in off if r[1] == "low_coverage"],
          "floor 0.0 calls no read low_coverage")

    fired = [r for r in hard if r[1] == "low_coverage"]
    check(len(fired) > 0,
          f"floor 0.9 fires on at least one read (fired on {len(fired):,})")
    check(all(a < 0.9 * L for _, _, L, a in fired),
          "every read floor 0.9 rejected is genuinely below 0.9 coverage")

    check(len(off) == len(dflt) == len(hard),
          "read count is identical at every floor (reads move class, none are lost)")
    check({r[0] for r in off} == {r[0] for r in hard},
          "the same read ids appear at every floor")

    low = [r for r in dflt if r[1] == "low_coverage"]
    kept = [r for r in dflt if r[1] == "assigned"]
    check(all(a >= 0.10 * L for _, _, L, a in kept),
          "nothing assigned at floor 0.1 sits below the floor")
    print(f"  info floor 0.1 moved {len(low):,} reads to low_coverage")

    # A floor can only ever move reads OUT of assigned, never into it.
    off_assigned = {r[0] for r in off if r[1] == "assigned"}
    dflt_assigned = {r[0] for r in dflt if r[1] == "assigned"}
    check(dflt_assigned <= off_assigned,
          "the floor only removes reads from assigned, never adds them")

    if failures:
        print(f"\nFAILED: {len(failures)} check(s)")
        sys.exit(1)
    print("\nattribution floor: all checks passed")


if __name__ == "__main__":
    main()
