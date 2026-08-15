#!/usr/bin/env python3
"""
Per-run acquisition metadata, read from the reads themselves.

Supplementary Table S9 reports a run identifier, a date and a basecalling model
for each library. Those were originally taken by eye from the first FASTQ
header, which is wrong for the date: **a dorado-emitted FASTQ is not sorted by
acquisition time**. For lowinput_s1_r1 the first record in the file starts at
2025-04-11T02:47Z while the earliest read of the run starts at
2025-04-09T22:34Z -- a two-day error in a published table, from a file whose
first line looked authoritative.

So the date has to be a minimum over every read, not a lookup. That means a full
pass over each FASTQ (~100 GB across the seven runs, tens of minutes); it is run
once and its output committed, not run per analysis.

Fields per run, all derived:

    reads          records seen
    run_id         MinKNOW acquisition id, from the RG tag
    model          basecalling model, from the RG tag
    first_read     earliest read start time (`st` tag), UTC -- the reported date
    last_read      latest read start time, UTC
    duration_h     last - first. Note this is the span of *retained* reads and
                   so is a lower bound on the acquisition window: a run that
                   produced nothing in its final hours reads as shorter than it
                   was.

The RG tag is asserted constant within a file rather than taken from the first
record, on the same principle that produced the date bug.

Usage:
    python3 bin/run_metadata.py --out results/summary/run_metadata.tsv
    python3 bin/run_metadata.py lowinput_s1_r1
"""
from __future__ import annotations

import argparse
import gzip
import sys
from datetime import datetime
from pathlib import Path

COLS = ["sample_id", "reads", "run_id", "model", "first_read", "last_read", "duration_h"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("replicates", nargs="*", help="default: every sample in the samplesheets")
    p.add_argument("--samplesheets", default="assets/samplesheets")
    p.add_argument("--out", default=None)
    return p.parse_args()


def samples(sheet_dir, sradir="sra"):
    """sample_id -> fastq path, from the samplesheets rather than a guess.

    Falls back to the --fetch_from_sra download cache (`<sample_id>.fastq.gz`),
    so this works from a clone whose reads came from the archive rather than
    from the samplesheet's local `fastq` column.
    """
    out = {}
    for p in sorted(Path(sheet_dir).glob("*.csv")):
        if p.name in ("test.csv", "raghavendra_2023.csv", "all.csv"):
            continue
        rows = [l for l in p.read_text().splitlines()
                if l.strip() and not l.startswith("#")]
        cols = rows[0].split(",")
        for line in rows[1:]:
            f = dict(zip(cols, line.split(",")))
            sid = f.get("sample_id")
            if not sid:
                continue
            local = Path(f["fastq"]) if f.get("fastq") else None
            cached = Path(sradir) / f"{sid}.fastq.gz"
            if local and local.is_file():
                out[sid] = local
            elif cached.is_file():
                out[sid] = cached
            elif local:
                out[sid] = local
    return out


def scan(path):
    n = 0
    first = last = None
    rgs = set()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4:
                continue
            n += 1
            for field in line.rstrip("\n").split("\t")[1:]:
                if field.startswith("st:Z:"):
                    t = field[5:]
                    if first is None or t < first:
                        first = t
                    if last is None or t > last:
                        last = t
                elif field.startswith("RG:Z:"):
                    rgs.add(field[5:])
    if len(rgs) > 1:
        sys.exit(f"error: {path} mixes {len(rgs)} read groups: {sorted(rgs)}. This "
                 f"script assumes one acquisition per file.")
    return n, first, last, (rgs.pop() if rgs else "")


def main():
    args = parse_args()
    root = Path(__file__).resolve().parent.parent
    sheets = root / args.samplesheets
    known = samples(sheets)
    reps = args.replicates or sorted(known)

    rows = []
    for sid in reps:
        fq = known.get(sid)
        if fq is None:
            sys.exit(f"error: {sid} is in no samplesheet")
        if not fq.is_absolute():
            fq = root / fq
        if not fq.is_file():
            sys.exit(f"error: {fq} not found")

        print(f"[run-metadata] {sid}: scanning {fq.name}", file=sys.stderr)
        n, first, last, rg = scan(fq)

        # RG is "<acquisition-id>_<model>"; the model itself contains
        # underscores, so split once from the left only.
        run_id, _, model = rg.partition("_")
        dur = ""
        if first and last:
            fmt = "%Y-%m-%dT%H:%M:%S.%f%z"
            dur = round((datetime.strptime(last, fmt)
                         - datetime.strptime(first, fmt)).total_seconds() / 3600, 1)

        rows.append(dict(sample_id=sid, reads=n, run_id=run_id, model=model,
                         first_read=first or "", last_read=last or "",
                         duration_h=dur))
        print(f"  {n:,} reads  {str(first)[:10]} -> {str(last)[:10]}  {dur} h",
              file=sys.stderr)

    text = "\t".join(COLS) + "\n" + "\n".join(
        "\t".join(str(r[c]) for c in COLS) for r in rows) + "\n"
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text)
        print(f"[run-metadata] -> {args.out}", file=sys.stderr)
    print(text)


if __name__ == "__main__":
    main()
