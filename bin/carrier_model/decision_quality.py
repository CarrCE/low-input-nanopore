#!/usr/bin/env python3
"""Live decision quality and first-chunk speed, by read class, for every run;
and lambda read statistics under both definitions (the paper's carrier class,
lambda-assigned only; and the model's lambda-derived reads).

Output: <outdir>/decision_quality.tsv, <outdir>/carrier_definitions.tsv
Usage: decision_quality.py <results> <workdir> <outdir> <instrument_root>
"""
import glob, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import SAMPLES, classes  # noqa: E402


def main(ror, work, out, iroot):
    ror, work, out, iroot = Path(ror), Path(work), Path(out), Path(iroot)
    q_rows, c_rows = [], []
    for s in SAMPLES:
        ids = pd.read_csv(work / f"{s}.read_ids.txt.gz", header=None, names=["read_id"])
        d = np.load(work / f"{s}.reads.npz"); act = np.load(work / f"{s}.instrument.npz")["action"]
        L = d["length"].astype(np.int64); R = classes(d)
        strict, like = R["carrier"], R["lambda_derived"]
        f = sorted(glob.glob(str(iroot / s / "AS_timings_*.csv*")) +
                   glob.glob(str(iroot / s / "adaptive_sampling/AS_timings_*.csv*")))[0]
        t = pd.read_csv(f, usecols=["read_id", "mean_qscore", "bases", "samples"]
                        ).drop_duplicates("read_id", keep="last")
        m = ids.merge(t, on="read_id", how="left")
        q = m["mean_qscore"].to_numpy(); b = m["bases"].to_numpy(); n = m["samples"].to_numpy()
        for lab, sel in (("lambda-derived ejected 450-2k", like & (act == 1) & (L >= 450) & (L < 2000)),
                         ("lambda-derived kept >=2k", like & (act == 2) & (L >= 2000)),
                         ("sample kept", R["sample"] & (act == 2))):
            q_rows.append({"sample_id": s, "class": lab, "reads": int(sel.sum()),
                           "decision_qscore_median": float(np.nanmedian(q[sel])),
                           "bases_at_decision_median": float(np.nanmedian(b[sel])),
                           "first_chunk_bases_per_s": float(np.nanmedian(b[sel] / n[sel] * 5000))})
        for lab, c in (("carrier class (lambda-assigned)", strict), ("lambda-derived (model)", like)):
            Lc = L[c]; e = Lc >= 2000
            c_rows.append({"sample_id": s, "carrier_definition": lab, "reads": int(c.sum()),
                           "mean_len": float(Lc.mean()), "median_len": float(np.median(Lc)),
                           "escaped_reads_frac": float(e.mean()),
                           "escaped_bases_frac": float(Lc[e].sum() / Lc.sum())})
        print(s, "done", flush=True)
    pd.DataFrame(q_rows).to_csv(out / "decision_quality.tsv", sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(c_rows).to_csv(out / "carrier_definitions.tsv", sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    main(*sys.argv[1:])
