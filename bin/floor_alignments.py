#!/usr/bin/env python3
"""Where the attribution floor's excluded reads align (Supplementary Section S2).
Run by `make flooralignments`.

`floor_complexity.py` says what the excluded reads are made of: mostly
low-complexity sequence, about half of it 1-3 bp tandem repeats, and in these
panels only the two yeast genomes hold simple repeats at all. That leaves the
mechanism one step short, because it does not show that the reads align to
those repeats. This script closes that step from the alignments themselves.

Two steps, because the second needs a BAM holding only these reads:

  --write-ids DIR   choose the reads and write DIR/<sample>.{excluded,control}.txt.
                    Excluded: every read the floor excludes (low_coverage with
                    the floor on). Control: community reads of at least 10 kb
                    that the floor keeps, up to --controls per sample, drawn
                    with --seed. bin/floor_alignments.sh then extracts their
                    alignments with samtools.
  (no --write-ids)  measure, reading those same id files, so the selection is
                    recorded rather than recomputed.

For each read it takes every alignment minimap2 reported, the organism of each
alignment's contig, and how much of the reference span the alignment covers
inside a simple repeat: homopolymers and CA or AT dinucleotide repeats of at
least 24 bp, the definition floor_complexity.py applies to whole genomes.

Reports, per experiment, for the excluded reads of each role and for the
controls:
  * reads, how many align at all, and alignments per read;
  * the share of aligned reference bases inside a simple repeat, excluded
    against control -- the number that decides the mechanism;
  * aligned and simple-repeat reference bases by organism;
  * how concentrated they are: the distinct repeat tracts carrying an aligned
    excluded read, and the tracts the most reads align to.

Writes <out>/floor_alignments.json and <out>/floor_alignments_loci.tsv.
"""
import argparse
import csv
import glob
import gzip
import json
import os
import random
import statistics
import sys
from collections import Counter, defaultdict

import pysam

from floor_complexity import GENOME_REPEAT, experiment_of

ROLES = ("sample", "carrier", "contaminant")


def contig_map(path):
    """contig -> (organism, role); also the accessions the reference set holds."""
    out, accs = {}, set()
    for r in csv.DictReader(open(path), delimiter="\t"):
        out[r["contig"]] = (r["organism"], r["role"])
        accs.add(r["accession"])
    return out, sorted(accs)


def repeat_intervals(refdir, refset, accessions):
    """contig -> sorted, disjoint (start, end) simple repeats, 0-based half-open."""
    out = {}
    for acc in accessions:
        fna = sorted(glob.glob(f"{refdir}/{refset}/genomes/ncbi_dataset/data/{acc}/*.fna"))
        if not fna:
            sys.exit(f"error: no genome for {acc} under {refdir}/{refset}")
        name, seq = None, []
        for line in open(fna[0]):
            if line.startswith(">"):
                if name is not None:
                    out[name] = find_repeats("".join(seq))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip().upper())
        if name is not None:
            out[name] = find_repeats("".join(seq))
    return out


def find_repeats(seq):
    return [(m.start(), m.end()) for m in GENOME_REPEAT.finditer(seq)]


def overlap(intervals, start, end):
    """Bases of [start, end) covered by the sorted, disjoint intervals."""
    n = 0
    for s, e in intervals:
        if e <= start:
            continue
        if s >= end:
            break
        n += min(e, end) - max(s, start)
    return n


def assignments(root, sample):
    out = {}
    with gzip.open(f"{root}/{sample}/competitive/{sample}.assignments.tsv.gz", "rt") as f:
        head = next(f).rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(head)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            out[p[ix["read_id"]]] = (p[ix["call"]], p[ix["role"]], p[ix["organism"]],
                                     int(p[ix["read_length"]]))
    return out


def select(on, off, controls, rng):
    """(excluded read ids, control read ids) for one sample."""
    excluded = [r for r, v in on.items() if v[0] == "low_coverage"]
    pool = [r for r, v in off.items()
            if v[0] == "assigned" and v[1] == "sample" and v[3] >= 10_000
            and on.get(r, ("",))[0] == "assigned"]
    return excluded, rng.sample(pool, min(len(pool), controls))


def median(xs):
    xs = [x for x in xs if x is not None]
    return statistics.median(xs) if xs else None


def read_rows(bam, want, cmap, repeats):
    """read id -> its alignments, aligned reference bases, and repeat overlap."""
    rows = {}
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as f:
        for a in f:
            if a.is_unmapped or a.query_name not in want:
                continue
            org = cmap.get(a.reference_name, ("unknown", "unknown"))[0]
            ivs = repeats.get(a.reference_name, [])
            span = a.reference_end - a.reference_start
            rep = overlap(ivs, a.reference_start, a.reference_end)
            d = rows.setdefault(a.query_name, {
                "alignments": 0, "ref_bases": 0, "repeat_bases": 0, "longest": 0,
                "by_org": Counter(), "repeat_by_org": Counter(), "loci": set()})
            d["alignments"] += 1
            d["ref_bases"] += span
            d["repeat_bases"] += rep
            d["longest"] = max(d["longest"], span)
            d["by_org"][org] += span
            d["repeat_by_org"][org] += rep
            if rep:
                for s, e in ivs:
                    if s < a.reference_end and e > a.reference_start:
                        d["loci"].add((a.reference_name, s, e))
    return rows


def summarise(rows, ids):
    aligned = [rows[r] for r in ids if r in rows]
    shares = [d["repeat_bases"] / d["ref_bases"] for d in aligned if d["ref_bases"]]
    by_org, repeat_by_org = Counter(), Counter()
    for d in aligned:
        by_org.update(d["by_org"])
        repeat_by_org.update(d["repeat_by_org"])
    return {
        "reads": len(ids),
        "with_an_alignment": len(aligned),
        "median_alignments_per_read": median([d["alignments"] for d in aligned]),
        "median_longest_alignment": median([d["longest"] for d in aligned]),
        "median_repeat_share_of_aligned_reference": median(shares),
        "reads_over_half_repeat": sum(s > 0.5 for s in shares),
        "aligned_reference_bases": sum(by_org.values()),
        "repeat_reference_bases": sum(repeat_by_org.values()),
        "aligned_bases_by_organism": dict(by_org.most_common(8)),
        "repeat_bases_by_organism": dict((k, v) for k, v in repeat_by_org.most_common(8) if v),
        # The share for one organism's own alignments. Compare excluded with
        # control here, not overall: the two sets align to different organisms,
        # and most of the panel holds no simple repeat for either to land in.
        "repeat_share_by_organism": {k: round(repeat_by_org[k] / v, 6)
                                     for k, v in by_org.most_common(8) if v},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--with", dest="on", required=True, help="run with the floor on")
    ap.add_argument("--without", dest="off", required=True, help="the same run, floor off")
    ap.add_argument("--contig-map", nargs="+", required=True, metavar="SAMPLE=PATH")
    ap.add_argument("--refset", nargs="+", required=True, metavar="SAMPLE=NAME")
    ap.add_argument("--ids", required=True, metavar="DIR",
                    help="where the id files are (or, with --write-ids, go)")
    ap.add_argument("--write-ids", action="store_true",
                    help="only choose the reads and write the id files")
    ap.add_argument("--bam", nargs="*", default=[], metavar="SAMPLE=PATH",
                    help="the extracted alignments; required unless --write-ids")
    ap.add_argument("--refdir", default="refs")
    ap.add_argument("--controls", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out")
    a = ap.parse_args()

    maps = dict(x.split("=", 1) for x in a.contig_map)
    sets = dict(x.split("=", 1) for x in a.refset)
    samples = sorted(maps)

    if a.write_ids:
        os.makedirs(a.ids, exist_ok=True)
        rng = random.Random(a.seed)
        for s in samples:
            excl, ctrl = select(assignments(a.on, s), assignments(a.off, s), a.controls, rng)
            for kind, ids in (("excluded", excl), ("control", ctrl)):
                with open(f"{a.ids}/{s}.{kind}.txt", "w") as f:
                    f.write("".join(r + "\n" for r in ids))
            print(f"{s}: {len(excl)} excluded, {len(ctrl)} controls")
        return

    if not a.bam or not a.out:
        sys.exit("error: --bam and --out are required unless --write-ids is given")
    bams = dict(x.split("=", 1) for x in a.bam)

    def ids_of(sample, kind):
        with open(f"{a.ids}/{sample}.{kind}.txt") as f:
            return [l.strip() for l in f if l.strip()]

    rows_all, loci, repeats_cache = {}, Counter(), {}
    picked = {"excluded": defaultdict(lambda: defaultdict(list)),
              "control": defaultdict(lambda: defaultdict(list))}
    for s in samples:
        cmap, accs = contig_map(maps[s])
        refset = sets[s]
        if refset not in repeats_cache:
            repeats_cache[refset] = repeat_intervals(a.refdir, refset, accs)
            n = sum(len(v) for v in repeats_cache[refset].values())
            print(f"{refset}: {n} simple repeats in the reference", file=sys.stderr)
        repeats = repeats_cache[refset]
        off = assignments(a.off, s)
        excl, ctrl = ids_of(s, "excluded"), ids_of(s, "control")
        rows = read_rows(bams[s], set(excl) | set(ctrl), cmap, repeats)
        rows_all[s] = rows
        e = experiment_of(s)
        for r in excl:
            role = off.get(r, ("", "unknown"))[1]
            picked["excluded"][e][role if role in ROLES else "other"].append(r)
            for c, st, en in rows.get(r, {}).get("loci", ()):
                loci[(e, c, st, en)] += 1
        picked["control"][e]["sample"].extend(ctrl)
        print(f"{s}: {len(excl)} excluded ({sum(r in rows for r in excl)} align), "
              f"{len(ctrl)} controls", file=sys.stderr)

    result = {"definitions": {
        "simple_repeat": "homopolymers and CA/AT dinucleotide repeats >= 24 bp,"
                         " as in floor_complexity.py",
        "controls": f"community reads >= 10 kb the floor keeps, up to {a.controls}"
                    f" per sample, seed {a.seed}",
        "repeat_share": "simple-repeat bases of the aligned reference span, over that span",
    }, "experiments": {}}
    for kind in ("excluded", "control"):
        for e, by_role in picked[kind].items():
            rows = {}
            for s in samples:
                if experiment_of(s) == e:
                    rows.update(rows_all[s])
            dest = result["experiments"].setdefault(e, {}).setdefault(kind, {})
            for role, ids in by_role.items():
                dest[role] = summarise(rows, ids)

    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/floor_alignments.json", "w") as f:
        json.dump(result, f, indent=1)
    with open(f"{a.out}/floor_alignments_loci.tsv", "w") as f:
        f.write("experiment\tcontig\tstart\tend\tlength\treads\n")
        for (e, c, st, en), n in loci.most_common():
            f.write(f"{e}\t{c}\t{st}\t{en}\t{en - st}\t{n}\n")
    print(json.dumps(result["experiments"], indent=1))
    print(f"{len(loci)} distinct simple-repeat tracts carry an aligned excluded read")


if __name__ == "__main__":
    main()
