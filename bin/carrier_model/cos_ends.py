#!/usr/bin/env python3
"""Where carrier reads begin and end on the lambda genome.

A ligation-library read starts at an adapted molecule end. Intact lambda has
ends only at its two cos sites (reference 0 and G), so where carrier reads
begin says how the carrier's adaptable ends divide between cos ends and
breaks. Observed: cos-start reads sit 5-15 bp outside the cos coordinate (the
12-nt overhang filled in by end-prep), so a cos start is an extrapolated 5' end
within [-60, +20] bp of the cos site, measured into the genome.

Input per sample: <workdir>/<s>.ends.paf.gz (minimap2, run_carrier_model.sh)
and <workdir>/<s>.lambda_ids.tsv (extract_reads.py). Reads are <id>_whole, or
<id>_head + <id>_tail (first/last 1.5 kb).

Read sets:
  hash        unbiased P=0.02 hash sample of all lambda-derived reads (assigned
              to lambda, or ambiguous with lambda among the tied organisms)
  hash_lt450, hash_450_2k, hash_ge2k   the same sample split by read length
  escaped     every lambda-derived read >= 2 kb

For escaped reads that begin at a cos site, full_length_frac is the fraction
whose 3' end reaches the opposite cos site: a lower bound on exp(-breaks per
genome), since read truncation also shortens reads.

Output: <outdir>/cos_ends.tsv, <outdir>/cos_5prime_positions.tsv (2-kb bins by
strand).

Usage: cos_ends.py <workdir> <outdir> <sample_id>...
"""
import gzip, math, sys
from pathlib import Path
import numpy as np

G = 48502
COS_WIN = (-60, 20)
MAX_CLIP = 100
P = 0.02


def hashed(rid):
    return int(rid[:8], 16) / 2**32 < P


def ends(rec):
    """Distances, measured into the genome, of the read's 5' end from the cos
    site it starts beside, and of its 3' end from the opposite cos site.
    None where the unaligned clip is too long to extrapolate across."""
    qlen, qs, qe, strand, ts, te = rec
    clip5, clip3 = qs, qlen - qe
    if strand == "+":
        d5, d3 = ts - clip5, G - (te + clip3)
    else:
        d5, d3 = G - (te + clip5), ts - clip3
    return (d5 if clip5 <= MAX_CLIP else None), (d3 if clip3 <= MAX_CLIP else None), strand


def is_cos(d):
    return d is not None and COS_WIN[0] <= d <= COS_WIN[1]


def wilson(k, n, z=1.96):
    if n == 0:
        return math.nan, math.nan
    ph = k / n; d = 1 + z * z / n
    c = ph + z * z / (2 * n); h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def analyse(sample, work):
    prim = {}
    with gzip.open(work / f"{sample}.ends.paf.gz", "rt") as fh:
        for line in fh:
            if "tp:A:P" not in line:
                continue
            f = line.split("\t")
            rec = (int(f[1]), int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8]))
            q = f[0]
            if q not in prim or rec[2] - rec[1] > prim[q][2] - prim[q][1]:
                prim[q] = rec
    sets = {k: [] for k in ("hash", "hash_lt450", "hash_450_2k", "hash_ge2k", "escaped")}
    hist = np.zeros((2, G // 2000 + 1), dtype=np.int64)
    with open(work / f"{sample}.lambda_ids.tsv") as fh:
        next(fh)
        for line in fh:
            rid, L, _by, _org = line.rstrip("\n").split("\t")
            L = int(L)
            five = prim.get(f"{rid}_whole") or prim.get(f"{rid}_head")
            three = prim.get(f"{rid}_whole") or prim.get(f"{rid}_tail")
            d5 = ends(five)[0] if five else None
            d3, s3 = (ends(three)[1], ends(three)[2]) if three else (None, None)
            s5 = five[3] if five else None
            r = (L, d5, d3, s5 == s3)
            if hashed(rid):
                sets["hash"].append(r)
                sets["hash_lt450" if L < 450 else "hash_450_2k" if L < 2000 else "hash_ge2k"].append(r)
                if five:
                    p5 = five[4] - five[1] if five[3] == "+" else five[5] + five[1]
                    hist[0 if five[3] == "+" else 1, min(max(p5, 0), G) // 2000] += 1
            if L >= 2000:
                sets["escaped"].append(r)
    rows = []
    for name, rs in sets.items():
        n = len(rs)
        m = [r for r in rs if r[1] is not None]
        k = sum(is_cos(r[1]) for r in m)
        lo, hi = wilson(k, len(m))
        cos_start = [r for r in m if is_cos(r[1]) and r[2] is not None and r[3]]
        full = sum(is_cos(r[2]) and 0.9 * G <= r[0] <= 1.1 * G for r in cos_start)
        rows.append({"sample_id": sample, "read_set": name, "reads": n,
                     "five_prime_placed": len(m), "placed_frac": len(m) / n if n else math.nan,
                     "cos_start": k, "f_cos": k / len(m) if m else math.nan,
                     "f_cos_lo": lo, "f_cos_hi": hi,
                     "cos_start_both_ends_placed": len(cos_start),
                     "full_length": full,
                     "full_length_frac": full / len(cos_start) if cos_start else math.nan,
                     "reads_gt_1.1G": sum(r[0] > 1.1 * G for r in rs)})
    return rows, hist


def main(work, out, samples):
    work, out = Path(work), Path(out)
    allrows, hists = [], {}
    for s in samples:
        rows, hist = analyse(s, work)
        allrows += rows; hists[s] = hist
        for r in rows:
            print(f"{s:16} {r['read_set']:12} n={r['reads']:>7} placed={r['placed_frac']:.3f} "
                  f"f_cos={r['f_cos']:.3f} [{r['f_cos_lo']:.3f},{r['f_cos_hi']:.3f}] "
                  f"full={r['full_length']}/{r['cos_start_both_ends_placed']} >1.1G={r['reads_gt_1.1G']}")
    cols = list(allrows[0])
    with open(out / "cos_ends.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in allrows:
            fh.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    with open(out / "cos_5prime_positions.tsv", "w") as fh:
        fh.write("bin_start\t" + "\t".join(f"{s}_{st}" for s in samples for st in "+-") + "\n")
        for i in range(G // 2000 + 1):
            fh.write(f"{i*2000}\t" + "\t".join(str(hists[s][j, i]) for s in samples for j in (0, 1)) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
