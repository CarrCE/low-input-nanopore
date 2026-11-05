#!/usr/bin/env python3
"""Compare the three arms of the reference-fragmentation test (`make fragmentation`,
docs/reference-fragmentation.md) on the D6311 replicates.

  A  the published run: contig-level assemblies of the exact yeast strains
  B  complete genomes of other strains (S288C, H99)
  C  B's genomes cut where A's assemblies break

B against C isolates fragmentation (same sequence, only the breaks differ). A
against B is the strain swap, and mixes fragmentation with strain divergence.
Each arm is given with the floor on and, for the floor's own effect, off.

Reports, competitive mode:
  * reads and bases per organism, per replicate and pooled, in each arm;
  * every read whose organism or call differs between B and C, by transition;
  * yeast reads whose primary alignment in B crosses one of C's cuts, and where
    C sent them (needs B's BAMs and the cut positions);
  * community bases the floor removes in each arm, by organism;
  * per experiment, base-weighted enrichment and bases per femtogram.

Writes <out>/fragmentation.json and <out>/fragmentation_by_organism.tsv.
"""
import argparse
import csv
import gzip
import json
import os
from collections import Counter, defaultdict

SAMPLES = ["lowinput_s1_r1", "lowinput_s1_r2", "lowinput_s1_r3"]
YEASTS = ["Saccharomyces cerevisiae", "Cryptococcus neoformans"]


def per_org(root):
    d = defaultdict(Counter)
    for r in csv.DictReader(open(f"{root}/summary/per_organism.tsv"), delimiter="\t"):
        if r["mode"] == "competitive" and r["sample_id"] in SAMPLES:
            d[r["organism"]]["reads"] += int(r["reads"])
            d[r["organism"]]["bases"] += int(r["bases"])
            d[(r["sample_id"], r["organism"])]["reads"] += int(r["reads"])
            d[(r["sample_id"], r["organism"])]["bases"] += int(r["bases"])
            d[r["organism"]]["role"] = r["role"]
    return d


def experiment(root):
    for r in csv.DictReader(open(f"{root}/summary/experiment_summary.tsv"), delimiter="\t"):
        if r["mode"] == "competitive" and r["experiment"] == "lowinput_s1":
            return {"enrichment_bases_mean": float(r["enrichment_bases_mean"]),
                    "bases_per_fg_mean": float(r["bases_per_fg_mean"])}


def calls(root, s):
    out = {}
    with gzip.open(f"{root}/{s}/competitive/{s}.assignments.tsv.gz", "rt") as f:
        next(f)
        for line in f:
            p = line.split("\t")
            out[p[0]] = (p[1], p[2])
    return out


def floor_loss(on, off):
    """Community bases per organism with the floor off minus on."""
    loss = {}
    for k, v in off.items():
        if isinstance(k, str) and v.get("role") == "sample":
            loss[k] = v["bases"] - on.get(k, Counter())["bases"]
    return loss


def cut_crossers(b_root, cuts_files, b_calls, c_calls):
    """Yeast reads whose B primary alignment spans a cut, and their fate in C."""
    import pysam
    cuts = defaultdict(list)
    for path in cuts_files:
        for row in csv.DictReader(open(path), delimiter="\t"):
            cuts[row["chrom"]].append(int(row["pos"]))
    fate = Counter()
    for s in SAMPLES:
        bam = f"{b_root}/{s}/alignments/{s}.qname.bam"
        with pysam.AlignmentFile(bam, "rb", check_sq=False) as f:
            for a in f:
                if a.is_unmapped or a.is_secondary or a.is_supplementary:
                    continue
                cs = cuts.get(a.reference_name)
                if not cs or not any(a.reference_start < x < a.reference_end for x in cs):
                    continue
                b = b_calls[s].get(a.query_name)
                c = c_calls[s].get(a.query_name)
                fate[(b[0] if b else "?", c[0] if c else "?")] += 1
    return [{"organism_B": k[0], "organism_C": k[1], "reads": v} for k, v in fate.most_common()]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    for arm in "ABC":
        ap.add_argument(f"--{arm.lower()}", required=True, help=f"arm {arm}, floor on")
        ap.add_argument(f"--{arm.lower()}-nofloor", required=True, help=f"arm {arm}, floor off")
    ap.add_argument("--cuts", nargs="+", required=True, help="cut-position TSVs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    roots = {"A": a.a, "B": a.b, "C": a.c}
    roots0 = {"A": a.a_nofloor, "B": a.b_nofloor, "C": a.c_nofloor}

    po = {k: per_org(v) for k, v in roots.items()}
    po0 = {k: per_org(v) for k, v in roots0.items()}

    b_calls = {s: calls(a.b, s) for s in SAMPLES}
    c_calls = {s: calls(a.c, s) for s in SAMPLES}
    trans, n_reads = Counter(), 0
    for s in SAMPLES:
        if b_calls[s].keys() != c_calls[s].keys():
            raise SystemExit(f"error: arms B and C hold different reads for {s}")
        n_reads += len(b_calls[s])
        for r, bv in b_calls[s].items():
            if bv != c_calls[s][r]:
                trans[(bv, c_calls[s][r])] += 1
    yeast_changed = sum(n for (bv, cv), n in trans.items()
                        if any(y in bv[0] or y in cv[0] for y in YEASTS))

    losses = {k: floor_loss(po[k], po0[k]) for k in "ABC"}
    floor = {}
    for k in "ABC":
        tot = sum(v for v in losses[k].values() if v > 0)
        y = sum(losses[k].get(o, 0) for o in YEASTS)
        floor[k] = {"community_bases_removed": tot,
                    "yeast_share_pct": 100.0 * y / tot if tot else 0.0,
                    "by_yeast": {o: losses[k].get(o, 0) for o in YEASTS}}

    result = {
        "samples": SAMPLES,
        "runs": {"floor_on": roots, "floor_off": roots0},
        "yeast_reads": {k: {o: {s: po[k].get((s, o), Counter())["reads"] for s in SAMPLES}
                            for o in YEASTS} for k in "ABC"},
        "yeast_bases": {k: {o: po[k].get(o, Counter())["bases"] for o in YEASTS} for k in "ABC"},
        "B_vs_C": {
            "reads": n_reads,
            "changed": sum(trans.values()),
            "changed_involving_a_yeast": yeast_changed,
            "transitions": [{"B": f"{bv[0]} [{bv[1]}]", "C": f"{cv[0]} [{cv[1]}]", "reads": n}
                            for (bv, cv), n in trans.most_common()],
        },
        "cut_crossing_reads": cut_crossers(a.b, a.cuts, b_calls, c_calls),
        "floor": floor,
        "experiment": {k: {"floor_on": experiment(roots[k]), "floor_off": experiment(roots0[k])}
                       for k in "ABC"},
    }
    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/fragmentation.json", "w") as f:
        json.dump(result, f, indent=1)
    orgs = sorted({o for k in "ABC" for o in po[k] if isinstance(o, str)})
    with open(f"{a.out}/fragmentation_by_organism.tsv", "w") as f:
        f.write("organism\t" + "\t".join(f"{k}_{m}" for k in "ABC" for m in ("reads", "bases")) + "\n")
        for o in orgs:
            f.write(o + "\t" + "\t".join(str(po[k].get(o, Counter())[m]) for k in "ABC"
                                         for m in ("reads", "bases")) + "\n")
    print(json.dumps({k: result[k] for k in ("yeast_reads", "floor", "experiment")}, indent=1))
    print(f"B vs C: {result['B_vs_C']['changed']} of {n_reads} reads change, "
          f"{yeast_changed} involving a yeast")
    print("reads crossing a cut (B -> C):", result["cut_crossing_reads"])


if __name__ == "__main__":
    main()
