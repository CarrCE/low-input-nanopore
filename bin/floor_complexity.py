#!/usr/bin/env python3
"""What the attribution floor excludes, by sequence: the complexity of the
community reads the 10% floor excludes, and the simple-repeat content of the
reference genomes (Supplementary Section S2). Run by `make floorcomplexity`.

The excluded reads are the community reads the floor removes: awarded to a
community organism (role=sample) with the floor off, and called low_coverage
with it on. This is the set `floor_effect.py` counts.

Reads are the deposited, human-masked reads. In masked reads, N marks sequence
the human screen removed (Supplementary Section S3). Each read is measured on
its unmasked bases joined together (its "core"), and only when the core holds
at least --min-core bases:

  * n_frac        share of the read that is N
  * zlib          zlib (level 6) compressed size of the core over its length;
                  ordinary sequence gives about 0.3, a simple repeat far less
  * low_entropy   share of the core's non-overlapping 64-bp windows whose
                  trinucleotide entropy is below 3 bits (of a possible 6)
  * tandem_frac   share of the core in tandem runs of at least 24 bp of a
                  unit of 1-3 bp, such as poly(A) or (CA)n
  * unit          the unit holding most of those bases, reported in its
                  canonical form (smallest rotation of it or its reverse
                  complement), so (CA)n, (AC)n, (TG)n and (GT)n are all "AC"

A read is "dominated by a 1-3 bp repeat" when tandem_frac >= 0.5.

Controls are community reads of at least 10 kb assigned with the floor off and
not excluded, up to --controls per sample, drawn with a fixed seed.

For the reference genomes, the simple repeats counted are homopolymers and
CA or AT dinucleotide repeats (either strand, any phase) of at least 24 bp:
their total length per Mb of genome, and the longest.

Writes <out>/floor_complexity.json (every number, named, per experiment) and
<out>/floor_complexity_references.tsv.
"""
import argparse
import csv
import glob
import gzip
import json
import math
import os
import random
import re
import statistics
import sys
import zlib
from collections import Counter, defaultdict

import pysam

TANDEM = re.compile(r"([ACGT]{1,3}?)\1{7,}")
GENOME_REPEAT = re.compile(r"A{24,}|C{24,}|G{24,}|T{24,}|(?:AC){12,}|(?:CA){12,}"
                           r"|(?:GT){12,}|(?:TG){12,}|(?:AT){12,}|(?:TA){12,}")
COMP = str.maketrans("ACGT", "TGCA")


def experiment_of(sample):
    return sample.rsplit("_r", 1)[0]


def samples_in(root):
    return sorted(os.path.basename(os.path.dirname(os.path.dirname(p)))
                  for p in glob.glob(f"{root}/*/competitive/*.assignments.tsv.gz"))


def assignments(root, sample):
    """read_id -> (call, role, read_length), in file order."""
    out = {}
    with gzip.open(f"{root}/{sample}/competitive/{sample}.assignments.tsv.gz", "rt") as f:
        head = next(f).rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(head)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            out[p[ix["read_id"]]] = (p[ix["call"]], p[ix["role"]], int(p[ix["read_length"]]))
    return out


def canonical(unit):
    rc = unit[::-1].translate(COMP)
    return min(s[i:] + s[:i] for s in (unit, rc) for i in range(len(s)))


def measure(seq, min_core):
    n = seq.count("N")
    d = {"length": len(seq), "n_frac": n / len(seq) if seq else 0.0}
    core = seq.replace("N", "")
    d["core"] = len(core)
    if len(core) < min_core:
        return d
    d["zlib"] = len(zlib.compress(core.encode(), 6)) / len(core)
    low = tot = 0
    for i in range(0, len(core) - 63, 64):
        w = core[i:i + 64]
        c = Counter(w[j:j + 3] for j in range(62))
        e = -sum(v / 62 * math.log2(v / 62) for v in c.values())
        tot += 1
        low += e < 3.0
    d["low_entropy"] = low / tot if tot else None
    units = Counter()
    for m in TANDEM.finditer(core):
        if m.end() - m.start() >= 24:
            units[canonical(m.group(1))] += m.end() - m.start()
    d["tandem_frac"] = sum(units.values()) / len(core)
    d["unit"] = units.most_common(1)[0][0] if units else None
    return d


def median(xs):
    xs = [x for x in xs if x is not None]
    return statistics.median(xs) if xs else None


def summarise(rows):
    meas = [r for r in rows if "zlib" in r]
    unmasked = [r for r in meas if r["n_frac"] == 0]
    dom = [r for r in meas if r["tandem_frac"] >= 0.5]

    def block(rs):
        return {"reads": len(rs),
                "median_length": median([r["length"] for r in rs]),
                "median_zlib": median([r["zlib"] for r in rs]),
                "median_low_entropy": median([r["low_entropy"] for r in rs]),
                "median_tandem_frac": median([r["tandem_frac"] for r in rs]),
                "dominated_by_1_3bp_repeat": sum(r["tandem_frac"] >= 0.5 for r in rs)}
    return {
        "reads": len(rows),
        "median_length": median([r["length"] for r in rows]),
        "with_any_n": sum(r["n_frac"] > 0 for r in rows),
        "more_than_half_n": sum(r["n_frac"] > 0.5 for r in rows),
        "median_n_frac": median([r["n_frac"] for r in rows]),
        "measurable": block(meas),
        "measurable_without_n": block(unmasked),
        "dominated_units": dict(Counter(r["unit"] for r in dom).most_common()),
    }


def reference_rows(refsets, refdir):
    rows = []
    for name, tsv in refsets:
        for r in csv.reader(open(tsv), delimiter="\t"):
            if not r or r[0].startswith("#") or r[0] == "organism":
                continue
            organism, acc = r[0], r[1]
            fna = sorted(glob.glob(f"{refdir}/{name}/genomes/ncbi_dataset/data/{acc}/*.fna"))
            if not fna:
                sys.exit(f"error: no genome for {acc} under {refdir}/{name}; "
                         "run the pipeline for this reference set first")
            seq = "".join(l.strip().upper() for l in open(fna[0]) if not l.startswith(">"))
            runs = [m.end() - m.start() for m in GENOME_REPEAT.finditer(seq)]
            rows.append({"reference_set": name, "organism": organism, "accession": acc,
                         "genome_bp": len(seq), "repeat_bp": sum(runs),
                         "repeat_bp_per_mb": round(sum(runs) / len(seq) * 1e6, 1),
                         "repeats": len(runs), "longest": max(runs) if runs else 0})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--with", dest="on", required=True, help="run with the floor on")
    ap.add_argument("--without", dest="off", required=True, help="the same run, floor off")
    ap.add_argument("--fastq", nargs="+", required=True, metavar="SAMPLE=PATH",
                    help="the deposited (masked) reads of each sample")
    ap.add_argument("--refset", nargs="+", default=[], metavar="NAME=TSV",
                    help="reference sets whose genomes to scan for simple repeats")
    ap.add_argument("--refdir", default="refs")
    ap.add_argument("--controls", type=int, default=400)
    ap.add_argument("--min-core", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    fastq = dict(x.split("=", 1) for x in a.fastq)
    samples = samples_in(a.on)
    if samples_in(a.off) != samples:
        sys.exit(f"error: {a.off} does not hold the same samples as {a.on}")
    missing = [s for s in samples if s not in fastq]
    if missing:
        sys.exit(f"error: no reads given for {', '.join(missing)}")

    random.seed(a.seed)
    rows = {"excluded": defaultdict(list), "control": defaultdict(list)}
    for s in samples:
        on, off = assignments(a.on, s), assignments(a.off, s)
        if on.keys() != off.keys():
            sys.exit(f"error: {s}: the runs with and without the floor hold different reads")
        excl = {r for r, v in on.items()
                if v[0] == "low_coverage" and off[r][0] == "assigned" and off[r][1] == "sample"}
        pool = [r for r, v in off.items()
                if v[0] == "assigned" and v[1] == "sample" and v[2] >= 10_000 and r not in excl]
        ctrl = set(random.sample(pool, min(len(pool), a.controls)))
        want = excl | ctrl
        found = 0
        with pysam.FastxFile(fastq[s]) as f:
            for rec in f:
                if rec.name in want:
                    d = measure(rec.sequence.upper(), a.min_core)
                    d["sample"] = s
                    rows["excluded" if rec.name in excl else "control"][experiment_of(s)].append(d)
                    found += 1
        if found != len(want):
            sys.exit(f"error: {s}: found {found} of {len(want)} reads in {fastq[s]}; "
                     "are these the reads the runs were made from?")
        print(f"{s}: {len(excl)} excluded, {len(ctrl)} controls", file=sys.stderr)

    refsets = [tuple(x.split("=", 1)) for x in a.refset]
    refs = reference_rows(refsets, a.refdir) if refsets else []
    result = {
        "definitions": {"min_core": a.min_core, "controls_per_sample": a.controls,
                        "seed": a.seed, "control_min_length": 10_000,
                        "low_entropy": "64-bp windows, trinucleotide entropy < 3 bits",
                        "tandem": "runs >= 24 bp of a 1-3 bp unit",
                        "dominated": "tandem_frac >= 0.5",
                        "genome_repeats": "homopolymers and CA/AT dinucleotide repeats >= 24 bp"},
        "experiments": {e: {"excluded": summarise(rows["excluded"][e]),
                            "control": summarise(rows["control"][e])}
                        for e in sorted(set(rows["excluded"]) | set(rows["control"]))},
        "references": refs,
    }
    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/floor_complexity.json", "w") as f:
        json.dump(result, f, indent=1)
    with open(f"{a.out}/floor_complexity_references.tsv", "w") as f:
        cols = ["reference_set", "organism", "accession", "genome_bp", "repeat_bp",
                "repeat_bp_per_mb", "repeats", "longest"]
        f.write("\t".join(cols) + "\n")
        for r in refs:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(json.dumps(result["experiments"], indent=1))


if __name__ == "__main__":
    main()
