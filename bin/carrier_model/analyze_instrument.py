#!/usr/bin/env python3
"""What the instrument did to each deposited read: decisions, escapes, timing,
recapture, and the two transient episodes.

Inputs: <workdir>/<s>.reads.npz (extract_reads.py; deposit order, roles, lengths)
        <workdir>/<s>.instrument.npz (join_instrument.py; same order)
        <instrument_root>/<s>/pore_activity_*.csv[.gz]
Outputs (<outdir>/):
  decisions.tsv        per replicate x read class: none / unblock / sequence shares
  escapes.tsv          per replicate: lambda-derived reads >= 2 kb by decision, and
                       what they cost
  instrument_time.tsv  per replicate x 15-min bin (run clock): yield, decision mix,
                       decision depth, basecaller latency, active channels
  recapture.tsv        per replicate: next read in the same channel after an
                       unblocked lambda-derived read vs after a read that ended
                       naturally
  <workdir>/<s>.pairs.npz  consecutive same-channel read pairs (recapture_ends.py)
  pores.tsv            pores available at the start of each run (pore_activity)

Usage: analyze_instrument.py <workdir> <outdir> <instrument_root>
"""
import csv, glob, gzip, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import BIN_S, SAMPLES, classes as read_classes, write_tsv as write  # noqa: E402

ACTION = {0: "none", 1: "unblock", 2: "sequence"}
SAMPLE_RATE = 5000


def classes(d):
    L = d["length"].astype(np.int64)
    R = read_classes(d)
    lam = R["lambda_derived"]
    return L, {
        "lambda-derived <450": lam & (L < 450),
        "lambda-derived 450-2k": lam & (L >= 450) & (L < 2000),
        "lambda-derived >=2k (escaped)": lam & (L >= 2000),
        "sample": R["sample"],
        "contaminant": R["contaminant"],
        "ambiguous, not lambda": R["ambiguous"] & ~R["ambiguous_lambda"],
        "unassigned": R["none"],
    }, lam, R["sample"]


def main(work, out, iroot):
    work, out, iroot = Path(work), Path(out), Path(iroot)
    dec_rows, esc_rows, t_rows, rc_rows, pore_rows = [], [], [], [], []
    for s in SAMPLES:
        d = np.load(work / f"{s}.reads.npz"); ins = np.load(work / f"{s}.instrument.npz")
        L, cls, lam, smp = classes(d)
        act = ins["action"]; ch = ins["channel"]; t0 = ins["start_time"]; dur = ins["duration"]
        decb = ins["dec_bases"]; tbc = ins["t_basecaller"]; dect = ins["dec_time"]

        # --- decision mix by class
        for name, m in cls.items():
            n = int(m.sum())
            row = {"sample_id": s, "class": name, "reads": n}
            for k, v in ACTION.items():
                row[f"{v}_frac"] = float((act[m] == k).mean()) if n else np.nan
            row["median_bases_at_decision"] = float(np.median(decb[m & (decb >= 0)])) if (m & (decb >= 0)).any() else np.nan
            row["median_read_len"] = float(np.median(L[m])) if n else np.nan
            dec_rows.append(row)

        # --- escapes: lambda-derived >= 2 kb, by decision, with the bases they cost
        esc = lam & (L >= 2000)
        Bc = L[lam].sum()
        r = {"sample_id": s, "escaped_reads": int(esc.sum()),
             "escaped_share_of_lambda_bases": L[esc].sum() / Bc}
        for k, v in ACTION.items():
            m = esc & (act == k)
            r[f"{v}_reads"] = int(m.sum()); r[f"{v}_share_of_lambda_bases"] = L[m].sum() / Bc
            r[f"{v}_median_len"] = float(np.median(L[m])) if m.any() else np.nan
        # unblocked but long: how far past the decision did the strand keep going
        m = esc & (act == 1) & (decb >= 0)
        r["unblocked_escape_median_bases_after_decision"] = float(np.median(L[m] - decb[m])) if m.any() else np.nan
        m = lam & (act == 1) & (L >= 450) & (L < 2000) & (decb >= 0)
        r["ejected_median_bases_after_decision"] = float(np.median(L[m] - decb[m]))
        r["sample_unblocked_frac"] = float((act[smp] == 1).mean())
        r["sample_unblocked_bases_frac"] = L[smp & (act == 1)].sum() / L[smp].sum()
        esc_rows.append(r)

        # --- time series on the run clock
        b = (t0 // BIN_S).astype(np.int64)
        for k in range(int(b.max()) + 1):
            m = b == k
            if m.sum() < 500:
                continue
            c = m & lam
            t_rows.append({
                "sample_id": s, "t_start_h": k * BIN_S / 3600, "reads": int(m.sum()),
                "bases": int(L[m].sum()), "active_channels": int(np.unique(ch[m]).size),
                "sample_per_lambda_derived": (m & smp).sum() / max(c.sum(), 1),
                "sample_base_frac": L[m & smp].sum() / L[m].sum(),
                "lambda_unblock_frac": float((act[c] == 1).mean()) if c.any() else np.nan,
                "lambda_none_frac": float((act[c] == 0).mean()) if c.any() else np.nan,
                "lambda_sequence_frac": float((act[c] == 2).mean()) if c.any() else np.nan,
                "lambda_escape_frac": float((L[c] >= 2000).mean()) if c.any() else np.nan,
                "median_bases_at_decision": float(np.median(decb[c & (decb >= 0)])) if (c & (decb >= 0)).any() else np.nan,
                "median_time_in_basecaller": float(np.nanmedian(tbc[c])) if c.any() else np.nan,
                "median_lambda_duration_s": float(np.median(dur[c])) if c.any() else np.nan,
            })

        # --- recapture: the next read in the same channel
        order = np.lexsort((t0, ch))
        chs, ts, ds = ch[order], t0[order], dur[order]
        same = chs[1:] == chs[:-1]
        gap = ts[1:] - (ts[:-1] + ds[:-1])
        prev = order[:-1][same]; nxt = order[1:][same]; g = gap[same]
        after_unblock = lam[prev] & (act[prev] == 1)
        after_natural = (act[prev] == 2) | ((act[prev] == 0) & ~lam[prev])
        row = {"sample_id": s}
        for lab, m in (("after_unblocked_lambda", after_unblock), ("after_natural_end", after_natural)):
            row[f"{lab}_n"] = int(m.sum())
            row[f"{lab}_median_gap_s"] = float(np.median(g[m]))
            row[f"{lab}_gap_lt_0.5s"] = float((g[m] < 0.5).mean())
            row[f"{lab}_gap_lt_2s"] = float((g[m] < 2).mean())
            row[f"{lab}_next_lambda_derived"] = float(lam[nxt[m]].mean())
        rc_rows.append(row)
        np.savez_compressed(work / f"{s}.pairs.npz", prev=prev, nxt=nxt, gap=g.astype(np.float32))

        # --- pores at start
        f = sorted(glob.glob(str(iroot / s / "pore_activity_*.csv*")))
        if f:
            T = {}
            with (gzip.open(f[0], "rt") if f[0].endswith(".gz") else open(f[0])) as fh:
                next(fh)
                for st, minute, samples in csv.reader(fh):
                    T.setdefault(int(minute), {})[st] = int(samples) / (SAMPLE_RATE * 60)
            avg = lambda ks, a, z: float(np.mean([sum(T.get(mm, {}).get(k, 0) for k in ks) for mm in range(a, z)]))
            pore_rows.append({"sample_id": s,
                              "available_min5_15": avg(["strand", "adapter", "pore", "unblocking"], 5, 15),
                              "sequencing_min5_15": avg(["strand"], 5, 15),
                              "available_h24": avg(["strand", "adapter", "pore", "unblocking"], 1440, 1450),
                              "available_h48": avg(["strand", "adapter", "pore", "unblocking"], 2880, 2890),
                              "available_h71": avg(["strand", "adapter", "pore", "unblocking"], 4260, 4270)})
        print(s, "done", flush=True)

    write(out / "decisions.tsv", dec_rows); write(out / "escapes.tsv", esc_rows)
    write(out / "instrument_time.tsv", t_rows); write(out / "recapture.tsv", rc_rows)
    if pore_rows:
        write(out / "pores.tsv", pore_rows)


if __name__ == "__main__":
    main(*sys.argv[1:])
