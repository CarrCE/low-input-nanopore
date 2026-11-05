#!/usr/bin/env python3
"""Every number behind the attribution-floor paragraph (Supplementary Section S2),
regenerated from three finished runs instead of quoted from memory.

  --with     a run with the floor on (the published analysis, results/)
  --without  the same run with the floor off (--min_aln_frac 0, results/nofloor)
  --q10      a run at --min_qscore 10 (results/q10); only its set of read IDs is
             used, as the reads that pass the Q10 filter

All three must cover the same samples, from the same reads. Competitive mode.

What it reports, pooled over the samples unless it says otherwise:
  * community (role=sample) reads and bases before any floor;
  * the alignment coverage of those reads (aligned bases / read length): the
    share above 98% and below 2%;
  * the reads the floor excludes (call "low_coverage" with the floor on), all
    of them and by role, and for the community: count, share of community
    reads, bases as a share of community bases, and mean length, against the
    mean length of all reads. These are the reads the assigner actually
    excluded -- no organism aligns to 10% of them -- which is stricter than
    attribution_threshold.py's view, which looks at the winning organism only;
  * who loses them: each excluded read is booked to the organism it was
    awarded to with the floor off, and the losses are tabulated by organism,
    with the carrier's loss also as a share of its bases;
  * per sample, the share of excluded reads the Q10 filter would also remove;
  * for one sample and organism (--focus, default D6311 r1 and S. cerevisiae):
    the alignment coverage of its awarded reads at all qualities and at Q10,
    and the share of its aligned bases and awarded read length kept under the
    floor;
  * per experiment, base-weighted enrichment and bases per femtogram with the
    floor off and on, and the largest relative change in any read-weighted
    quantity.

Outputs <out>/floor_effect.json (every number, named) and
<out>/floor_effect_by_organism.tsv.
"""
import argparse
import csv
import glob
import gzip
import json
import os
import sys
from collections import Counter, defaultdict

READ_WEIGHTED = ["reads_per_fg_mean", "enrichment_reads_mean"]


def samples_in(root):
    return sorted(os.path.basename(os.path.dirname(os.path.dirname(p)))
                  for p in glob.glob(f"{root}/*/competitive/*.assignments.tsv.gz"))


def assignments(root, sample):
    """read_id -> (organism, call, role, read_length, aligned_bases)"""
    out = {}
    with gzip.open(f"{root}/{sample}/competitive/{sample}.assignments.tsv.gz", "rt") as f:
        head = next(f).rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(head)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            out[p[ix["read_id"]]] = (p[ix["organism"]], p[ix["call"]], p[ix["role"]],
                                     int(p[ix["read_length"]]), int(p[ix["aligned_bases"]]))
    return out


def read_ids(root, sample):
    with gzip.open(f"{root}/{sample}/competitive/{sample}.assignments.tsv.gz", "rt") as f:
        next(f)
        return {line.split("\t", 1)[0] for line in f}


def experiments(root):
    rows = {}
    for r in csv.DictReader(open(f"{root}/summary/experiment_summary.tsv"), delimiter="\t"):
        if r["mode"] == "competitive":
            rows[r["experiment"]] = r
    return rows


def pct(a, b):
    return 100.0 * a / b if b else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--with", dest="on", required=True)
    ap.add_argument("--without", dest="off", required=True)
    ap.add_argument("--q10", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--focus-sample", default="lowinput_s1_r1")
    ap.add_argument("--focus-organism", default="Saccharomyces cerevisiae")
    a = ap.parse_args()

    samples = samples_in(a.on)
    for root in (a.off, a.q10):
        if samples_in(root) != samples:
            sys.exit(f"error: {root} does not hold the same samples as {a.on}: "
                     f"{samples_in(root)} vs {samples}")

    comm_reads = comm_bases = 0
    cov_hi = cov_lo = 0
    all_reads = all_bases = 0
    excl_n = excl_bases = 0
    loss = defaultdict(Counter)          # organism -> reads, bases
    role_of = {}
    off_bases_by_org = Counter()
    q10_removed = {}
    focus = {}
    excl_by_role = defaultdict(Counter)  # sample / carrier / contaminant / ambiguous

    for s in samples:
        on, off = assignments(a.on, s), assignments(a.off, s)
        if on.keys() != off.keys():
            sys.exit(f"error: {s}: the runs with and without the floor hold different reads")
        q10 = read_ids(a.q10, s)
        excluded = [r for r, v in on.items() if v[1] == "low_coverage"]
        for r, (org, call, role, rlen, aln) in off.items():
            all_reads += 1; all_bases += rlen
            if call == "assigned":
                off_bases_by_org[org] += rlen
                role_of[org] = role
            if call == "assigned" and role == "sample":
                comm_reads += 1; comm_bases += rlen
                c = aln / rlen if rlen else 0.0
                cov_hi += c > 0.98
                cov_lo += c < 0.02
        for r in excluded:
            org, call, role, rlen, _ = off[r]
            # Ambiguous reads already carry "ambiguous:<a>|<b>" as their organism.
            key = org if call in ("assigned", "ambiguous") else f"{call}: {org}"
            role_of.setdefault(key, role if call == "assigned" else call)
            loss[key]["reads"] += 1; loss[key]["bases"] += rlen
            excl_n += 1; excl_bases += rlen
            excl_by_role[role_of[key]]["reads"] += 1
            excl_by_role[role_of[key]]["bases"] += rlen
        q10_removed[s] = pct(sum(r not in q10 for r in excluded), len(excluded))

        if s == a.focus_sample:
            fo = a.focus_organism
            sel = [(rlen, aln, r) for r, (org, call, _, rlen, aln) in off.items()
                   if org == fo and call == "assigned"]
            sel_q = [(rlen, aln) for rlen, aln, r in sel if r in q10]
            kept = [(rlen, aln) for r, (org, call, _, rlen, aln) in on.items()
                    if org == fo and call == "assigned"]
            focus = {
                "sample": s, "organism": fo,
                "coverage_all_reads_pct": pct(sum(x[1] for x in sel), sum(x[0] for x in sel)),
                "coverage_q10_reads_pct": pct(sum(x[1] for x in sel_q), sum(x[0] for x in sel_q)),
                "aligned_bases_kept_pct": pct(sum(x[1] for x in kept), sum(x[1] for x in sel)),
                "awarded_length_lost_pct": 100 - pct(sum(x[0] for x in kept), sum(x[0] for x in sel)),
            }

    comm_loss_bases = sum(v["bases"] for k, v in loss.items() if role_of.get(k) == "sample")
    by_org = []
    for k, v in sorted(loss.items(), key=lambda kv: -kv[1]["bases"]):
        by_org.append({
            "organism": k, "role": role_of.get(k, ""), "reads": v["reads"], "bases": v["bases"],
            "share_of_excluded_community_bases_pct":
                pct(v["bases"], comm_loss_bases) if role_of.get(k) == "sample" else None,
            "share_of_own_bases_pct": pct(v["bases"], off_bases_by_org[k]) if k in off_bases_by_org else None,
        })

    exp_off, exp_on = experiments(a.off), experiments(a.on)
    exps = {}
    for e in sorted(exp_on):
        o, n = exp_off[e], exp_on[e]
        exps[e] = {
            "enrichment_bases_without": float(o["enrichment_bases_mean"]),
            "enrichment_bases_with": float(n["enrichment_bases_mean"]),
            "bases_per_fg_without": float(o["bases_per_fg_mean"]),
            "bases_per_fg_with": float(n["bases_per_fg_mean"]),
            "enrichment_bases_change_pct":
                pct(float(n["enrichment_bases_mean"]) - float(o["enrichment_bases_mean"]),
                    float(o["enrichment_bases_mean"])),
            "max_read_weighted_change_pct": max(
                abs(pct(float(n[c]) - float(o[c]), float(o[c]))) for c in READ_WEIGHTED),
        }

    result = {
        "samples": samples,
        "runs": {"with_floor": a.on, "without_floor": a.off, "q10": a.q10},
        "community_before_floor": {"reads": comm_reads, "bases": comm_bases},
        "coverage_above_98pct_share_pct": pct(cov_hi, comm_reads),
        "coverage_below_2pct_share_pct": pct(cov_lo, comm_reads),
        "excluded": {
            "all_reads": excl_n,
            "all_bases": excl_bases,
            "mean_length_bp": excl_bases / excl_n if excl_n else 0,
            "by_role": {k: dict(v) for k, v in sorted(excl_by_role.items())},
            "community_reads": excl_by_role["sample"]["reads"],
            "community_bases": comm_loss_bases,
            "share_of_community_reads_pct": pct(excl_by_role["sample"]["reads"], comm_reads),
            "share_of_community_bases_pct": pct(comm_loss_bases, comm_bases),
            "community_mean_length_bp": (comm_loss_bases / excl_by_role["sample"]["reads"]
                                         if excl_by_role["sample"]["reads"] else 0),
        },
        "mean_read_length_all_reads_bp": all_bases / all_reads if all_reads else 0,
        "q10_removes_pct_of_excluded_by_sample": q10_removed,
        "focus": focus,
        "experiments": exps,
        "by_organism": by_org,
    }
    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/floor_effect.json", "w") as f:
        json.dump(result, f, indent=1)
    with open(f"{a.out}/floor_effect_by_organism.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["organism", "role", "reads", "bases",
                    "share_of_excluded_community_bases_pct", "share_of_own_bases_pct"])
        for r in by_org:
            w.writerow([r["organism"], r["role"], r["reads"], r["bases"],
                        "" if r["share_of_excluded_community_bases_pct"] is None
                        else f'{r["share_of_excluded_community_bases_pct"]:.2f}',
                        "" if r["share_of_own_bases_pct"] is None else f'{r["share_of_own_bases_pct"]:.3f}'])
    print(json.dumps({k: result[k] for k in ("community_before_floor", "excluded",
                                              "q10_removes_pct_of_excluded_by_sample",
                                              "focus", "experiments")}, indent=1))


if __name__ == "__main__":
    main()
