#!/usr/bin/env python3
"""Two tests that need where reads start on lambda.

1. Recapture. For consecutive reads in the same channel where BOTH are in the
   P=0.02 hash sample (so both have mapped 5' ends), does the second read begin
   at the same molecule end as the first (same strand, 5' within +/-20 bp)?
   Internal starts are used for the rate, since an internal break position is
   nearly unique to one molecule; cos starts coincide by chance. Chance level is
   estimated from pairs separated by a natural end (the molecule went through).

2. Why kept lambda-derived reads were kept. 5' unaligned clip and decision
   depth for long reads by decision (unblock vs sequence vs none).

Inputs: <workdir>/<s>.{reads,instrument,pairs}.npz, <workdir>/<s>.ends.paf.gz,
        <results>/<s>/competitive/<s>.readlengths.tsv.gz
Output: <outdir>/recapture_ends.tsv, <outdir>/keep_decisions_5prime.tsv
Usage: recapture_ends.py <results> <workdir> <outdir>
"""
import gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import SAMPLES  # noqa: E402

G = 48502
P = 0.02
TOL = 20
COS_WIN = (-60, 20)


def write(path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, float_format="%.6g")


def five_prime(work, s):
    """read_id -> (strand, ref coord of 5' end, distance from upstream cos, 5' clip)."""
    prim = {}
    with gzip.open(work / f"{s}.ends.paf.gz", "rt") as fh:
        for line in fh:
            if "tp:A:P" not in line:
                continue
            f = line.split("\t")
            q = f[0]
            if not (q.endswith("_whole") or q.endswith("_head")):
                continue
            rec = (int(f[1]), int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8]))
            if q not in prim or rec[2] - rec[1] > prim[q][2] - prim[q][1]:
                prim[q] = rec
    out = {}
    for q, (qlen, qs, qe, strand, ts, te) in prim.items():
        rid = q.rsplit("_", 1)[0]
        p5 = ts - qs if strand == "+" else te + qs
        dcos = p5 if strand == "+" else G - p5
        out[rid] = (strand, p5, dcos, qs)
    return out


def main(ror, work, out):
    ror, work, out = Path(ror), Path(work), Path(out)
    rc_rows, kp_rows = [], []
    for s in SAMPLES:
        ids = pd.read_csv(ror / s / "competitive" / f"{s}.readlengths.tsv.gz", sep="\t",
                          usecols=["read_id"])["read_id"].to_numpy()
        fp = five_prime(work, s)
        ins = np.load(work / f"{s}.instrument.npz"); act = ins["action"]; decb = ins["dec_bases"]
        d = np.load(work / f"{s}.reads.npz"); L = d["length"]
        pr = np.load(work / f"{s}.pairs.npz")
        prev, nxt, gap = pr["prev"], pr["nxt"], pr["gap"]
        hashed = np.array([int(r[:8], 16) / 2**32 < P for r in ids])
        both = hashed[prev] & hashed[nxt]
        rows = {}
        for lab, sel in (("after_unblock", act[prev] == 1), ("after_natural_end", act[prev] != 1)):
            m = both & sel & (gap < 5)
            n = k = n_int = k_int = 0
            for a, b in zip(prev[m], nxt[m]):
                fa, fb = fp.get(ids[a]), fp.get(ids[b])
                if fa is None or fb is None:
                    continue
                n += 1
                same = fa[0] == fb[0] and abs(fa[1] - fb[1]) <= TOL
                k += same
                if not (COS_WIN[0] <= fa[2] <= COS_WIN[1]):
                    n_int += 1; k_int += same
            rows[lab] = (n, k, n_int, k_int)
        r = {"sample_id": s}
        for lab, (n, k, n_int, k_int) in rows.items():
            r[f"{lab}_pairs"] = n; r[f"{lab}_same_end_frac"] = k / n if n else np.nan
            r[f"{lab}_internal_pairs"] = n_int
            r[f"{lab}_internal_same_end_frac"] = k_int / n_int if n_int else np.nan
        rc_rows.append(r)

        # 5' clip by decision for long reads with a mapped 5' end
        for lab, lo, hi in (("450-2k", 450, 2000), (">=2k", 2000, 10**9)):
            for code, name in ((1, "unblock"), (2, "sequence"), (0, "none")):
                sel = np.where((act == code) & (L >= lo) & (L < hi))[0]
                clips = [fp[ids[i]][3] for i in sel if ids[i] in fp]
                dep = decb[sel]; dep = dep[dep >= 0]
                if len(clips) < 20:
                    continue
                c = np.array(clips)
                kp_rows.append({"sample_id": s, "length_class": lab, "decision": name,
                                "reads_with_5prime": len(c),
                                "clip_median": float(np.median(c)),
                                "clip_gt_100_frac": float((c > 100).mean()),
                                "clip_gt_200_frac": float((c > 200).mean()),
                                "bases_at_decision_median": float(np.median(dep)) if len(dep) else np.nan})
        print(s, rows, flush=True)
    write(out / "recapture_ends.tsv", rc_rows)
    write(out / "keep_decisions_5prime.tsv", kp_rows)


if __name__ == "__main__":
    main(*sys.argv[1:])
