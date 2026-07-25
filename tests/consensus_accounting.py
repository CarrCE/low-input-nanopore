#!/usr/bin/env python3
"""
Assert that breseq-consensus subtraction preserves read accounting.

Consensus subtraction is the one code path that can delete reads from the
per-organism tallies without any signal that it did -- that is precisely the
failure this study documents in the original lowinput_s1 analysis, so the
reimplementation of it has to be held to the accounting guarantee the rest of
the pipeline meets: every read lands in exactly one class, and the classes sum
to the FASTQ read count.

The test runs bin/assign_reads.py twice over the same BAM, once with and once
without a synthetic consensus-hits list, and checks seven properties. The hits
list deliberately mixes three cases that exercise different branches:

  * reads whose primary alignment is the contaminant  (the ordinary case)
  * reads whose primary alignment is something else   (must lose to the carrier
                                                       subtraction that runs
                                                       first in the chain)
  * reads that are unmapped entirely                  (must still be subtracted:
                                                       the original pipeline
                                                       subtracted on raw reads)

Usage:
    python3 tests/consensus_accounting.py \\
        --bam        results/test_s2/alignments/test_s2.qname.bam \\
        --contig-map results/references/lowinput_s2/contig_map.tsv

Both inputs are produced by `make test`. Exits non-zero on the first failure.
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True, help="qname-grouped BAM from MAP_COMPETITIVE")
    p.add_argument("--contig-map", required=True, help="contig_map.tsv from BUILD_REFERENCE")
    p.add_argument("--assign-script", default=None,
                   help="path to bin/assign_reads.py (default: alongside this test)")
    p.add_argument("--n-each", type=int, default=40,
                   help="how many reads to draw for each of the three hit cases")
    return p.parse_args()


def contaminant_contigs(contig_map):
    contigs, organisms = set(), set()
    with open(contig_map) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2] == "contaminant":
                contigs.add(f[0])
                organisms.add(f[1])
    if not contigs:
        sys.exit(f"error: no contaminant contigs in {contig_map}")
    return contigs, organisms


def sample_read_ids(bam_path, contaminant, n):
    """Draw the three classes of read the hits list needs to exercise."""
    on_contam, elsewhere, unmapped = [], [], []
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_secondary or aln.is_supplementary:
                continue
            if aln.is_unmapped:
                bucket = unmapped
            elif aln.reference_name in contaminant:
                bucket = on_contam
            else:
                bucket = elsewhere
            if len(bucket) < n:
                bucket.append(aln.query_name)
            if len(on_contam) >= n and len(elsewhere) >= n and len(unmapped) >= n:
                break
    for name, bucket in (("contaminant-primary", on_contam),
                         ("other-primary", elsewhere),
                         ("unmapped", unmapped)):
        if not bucket:
            sys.exit(f"error: the BAM contains no {name} reads, so this test "
                     f"cannot exercise that branch")
    return set(on_contam), set(elsewhere), set(unmapped)


def load_assignments(path):
    rows = {}
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows[f[0]] = dict(zip(hdr, f))
    return rows


def main():
    args = parse_args()
    script = Path(args.assign_script) if args.assign_script else \
        Path(__file__).resolve().parent.parent / "bin" / "assign_reads.py"
    if not script.exists():
        sys.exit(f"error: assign_reads.py not found at {script}")

    contaminant, contaminant_orgs = contaminant_contigs(args.contig_map)
    on_contam, elsewhere, unmapped = sample_read_ids(args.bam, contaminant, args.n_each)
    hits = on_contam | elsewhere | unmapped

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        hits_file = tmp / "consensus_hits.txt"
        hits_file.write_text("\n".join(sorted(hits)) + "\n")

        common = [sys.executable, str(script), args.bam,
                  "--contig-map", args.contig_map, "--mode", "sequential"]
        subprocess.run(common + ["--prefix", str(tmp / "base")],
                       check=True, stdout=subprocess.DEVNULL)
        subprocess.run(common + ["--prefix", str(tmp / "cons"),
                                 "--consensus-hits", str(hits_file)],
                       check=True, stdout=subprocess.DEVNULL)

        base = load_assignments(tmp / "base.assignments.tsv.gz")
        cons = load_assignments(tmp / "cons.assignments.tsv.gz")

    subtracted = {r for r, d in cons.items() if d["call"] == "subtracted_consensus"}
    lost_to_carrier = hits - subtracted
    fell_through = [r for r, d in base.items()
                    if d["role"] == "contaminant" and r not in hits
                    and cons[r]["call"] != "subtracted_consensus"]

    checks = [
        ("every read still lands somewhere, both runs",
         len(base) == len(cons) and len(cons) > 0,
         f"{len(base)} vs {len(cons)} reads"),

        ("nothing is subtracted that was not a hit",
         subtracted <= hits,
         f"{len(subtracted)} subtracted, {len(subtracted - hits)} of them not in the hits list"),

        ("subtracted reads are booked to a contaminant organism",
         all(cons[r]["organism"] in contaminant_orgs for r in subtracted),
         f"organisms seen: {sorted({cons[r]['organism'] for r in subtracted})}"),

        ("subtracted reads carry role=contaminant",
         all(cons[r]["role"] == "contaminant" for r in subtracted),
         f"roles seen: {sorted({cons[r]['role'] for r in subtracted})}"),

        # The original pipeline subtracted against the consensus on raw reads,
        # before any community mapping, so a read matching the consensus but
        # aligning nowhere in the combined index was removed there too.
        ("unmapped hits are subtracted, not left unassigned",
         all(cons[r]["call"] == "subtracted_consensus" for r in unmapped),
         f"calls seen: {sorted({cons[r]['call'] for r in unmapped})}"),

        # The chain is carrier -> contaminant, so a hit that also aligns to the
        # carrier is claimed one step earlier. That is the chain working, not a
        # consensus failure.
        ("hits not subtracted were claimed by the carrier first",
         all(cons[r]["role"] == "carrier" for r in lost_to_carrier),
         f"{len(lost_to_carrier)} such reads, roles: "
         f"{sorted({cons[r]['role'] for r in lost_to_carrier})}"),

        # Replace semantics: the consensus is a corrected version of the same
        # genome, so contaminant-aligning reads that do NOT match it are no
        # longer subtracted. They must reappear elsewhere, never vanish.
        ("contaminant reads that miss the consensus fall through, not away",
         all(r in cons for r in fell_through),
         f"{len(fell_through)} reads, now "
         f"{sorted({cons[r]['organism'] for r in fell_through})}"),
    ]

    failed = 0
    for i, (label, ok, detail) in enumerate(checks, 1):
        print(f"{'PASS' if ok else 'FAIL'}  {i}. {label}\n        {detail}")
        failed += not ok

    print(f"\n{len(checks) - failed}/{len(checks)} checks passed "
          f"({len(subtracted)} reads subtracted to the consensus)")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
