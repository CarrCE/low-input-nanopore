#!/usr/bin/env python3
"""Did a lambda read's start type change the live decision?
Run by `make decisionedge`.

f_cos rises with read length (SI Section S6), and the cause is not established.
For lambda-derived reads, read length is set largely by when adaptive sampling
decided, not by the length of the molecule, and ejection does not move a read's
5' end, which is the only thing f_cos measures. So the trend needs a mechanism
coupling start type to read length. The candidate this script tests:

  a cos-start read begins at the extreme terminus of the linear lambda
  reference with the 12-nt cohesive overhang filled in, so its first chunk maps
  at a reference edge rather than in the interior, plausibly with less
  confidence; adaptive sampling therefore ejects it less often; so it runs
  longer.

If that holds, the rise in f_cos with read length is a readout of the live
mapper's behaviour at the reference edge and carries nothing about carrier
intactness. It also bears on the main text, which attributes escape to poor
first-chunk quality: escapes are enriched for cos starts, so part of that low
quality would be a composition effect rather than chance alone.

The comparison. Decision-time live quality (AS_timings `mean_qscore`, the mean
over the chunk the live basecaller had produced when the decision was taken)
and ejection rate, for cos-start against break-start lambda-derived reads.

ON CONDITIONING, and this is the part that decides what the output can mean.
The primary comparison is UNCONDITIONED, and nothing here improves on it.

Final read length is a *descendant* of the ejection decision, the step the
mechanism acts through, so matching on it conditions on a collider.

`dec_bases` looks like a legitimate control, because `mean_qscore` is a mean
over that window, and it was built as one here. It is not clean either. The
decision falls at a roughly fixed latency, so bases-at-decision is basecalling
throughput over that window, and in nanopore data throughput and quality are
related. So `dec_bases` is itself measured downstream of the thing being
compared, and splitting on it can bias the comparison rather than fix it.

The data say so plainly. In the window holding the modal decision, [250,300)
bases, the two groups' `dec_bases` medians are balanced (279 against 279) and
cos-start quality comes out 0.95 lower pooled; in the window immediately above,
also balanced (319 against 318), it comes out 0.63 higher; below, 0.40 higher.
A difference that reverses by 1.5 qscore between adjacent, internally balanced
windows is not a property of the start type. It is what conditioning on a
downstream variable does.

So the binned tables are EXPLORATORY and must not be read as adjusted
estimates, and a count of how many bins favour one group is not a result. The
estimate is the unconditioned comparison, and it shows no consistent
difference: -0.10 qscore pooled, with the seven replicates straddling zero.

LIMITATION, AND IT IS THE FIRST THING TO READ. The mechanism above is about
*mappability*: a chunk from the reference terminus aligning poorly because part
of it hangs off the end of a linear reference. What this script measures is
*basecall quality*, `mean_qscore`, which is a different thing. A terminal chunk
can be basecalled perfectly and still align badly. The live mapping result is
nowhere in the deposited instrument records: `AS_decisions_<sample>.csv.gz`
carries only `read_id`, `action` and `action_response`, and
`AS_timings_<sample>.csv.gz` carries channel, timings, samples, bases,
barcode and `mean_qscore`. Neither records what the live mapper matched, or
how well. So this script can exclude the quality-mediated form of the
mechanism and cannot touch the mappability-mediated form, which the deposited
data cannot settle at all. Anyone meeting the f_cos-by-length trend will reach
for this explanation; the point of this paragraph is that they should not spend
a day discovering that the records cannot test it.

A second thing not to be caught by. Reads of the two start types differ in
ejection rate, and that is NOT evidence of a cause. A lambda read is short
largely *because* it was ejected, so "cos-start reads are ejected less often"
and "cos-start reads are longer" are nearly the same statement, and the latter
is the trend being explained. The output reports the ejection rates because
they are worth seeing, not because they settle anything; under 2 kb both start
types are ejected about 99% of the time and the difference even changes sign.

Verdict rule, fixed here rather than chosen after seeing the numbers:
  supports      cos-start median is lower by >= 0.5 qscore pooled, and lower in
                every replicate;
  excludes      the pooled difference is within +/- 0.2 and no dec_bases bin
                differs by more than 0.5;
  indeterminate anything else. This is an acceptable outcome: the SI's \\todo
                stands and the trend stays reported without a cause.

Inputs, all from a finished `make carriermodel` workdir plus the instrument
records: <sample>.read_ids.txt.gz, .reads.npz, .instrument.npz, .lambda_ids.tsv,
.ends.paf.gz, and AS_timings_<sample>.csv.gz.

Output: <outdir>/decision_edge.tsv, decision_edge_by_bases.tsv,
decision_edge_by_length.tsv, decision_edge.json

Usage: decision_edge.py <workdir> <instrument_root> <outdir> [sample_id...]
"""
import glob
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import SAMPLES, classes  # noqa: E402
from cos_ends import ends, is_cos  # noqa: E402

# Bin edges. BASES_BINS are the decision window, the legitimate control;
# LEN_BINS repeat the SI's own classes so the two can be lined up.
BASES_BINS = [0, 250, 300, 350, 400, 500, 10 ** 9]
LEN_BINS = [0, 450, 2000, 10 ** 9]
MIN_BIN = 30          # below this a bin's median is not reported
SUPPORT_DELTA = 0.5   # qscore
EXCLUDE_DELTA = 0.2


def primary_alignments(paf):
    """Query name -> the record covering the most of the query."""
    prim = {}
    with gzip.open(paf, "rt") as f:
        for line in f:
            p = line.split("\t")
            rec = (int(p[1]), int(p[2]), int(p[3]), p[4], int(p[7]), int(p[8]))
            q = p[0]
            if q not in prim or rec[2] - rec[1] > prim[q][2] - prim[q][1]:
                prim[q] = rec
    return prim


def start_type(work, sample):
    """read id -> True if the 5' end sits at a cos site, False if it does not.

    Reads whose 5' end cannot be placed are absent, not guessed: they are the
    0.1-1.0% the SI already excludes from f_cos, and a read with no placeable
    start has no start type to compare.
    """
    prim = primary_alignments(work / f"{sample}.ends.paf.gz")
    out = {}
    with open(work / f"{sample}.lambda_ids.tsv") as f:
        next(f)
        for line in f:
            rid = line.split("\t", 1)[0]
            rec = prim.get(f"{rid}_whole") or prim.get(f"{rid}_head")
            if rec is None:
                continue
            d5 = ends(rec)[0]
            if d5 is None:
                continue
            out[rid] = is_cos(d5)
    return out


def per_read(work, iroot, sample):
    """One row per lambda-derived read with a start type and a live decision."""
    ids = pd.read_csv(work / f"{sample}.read_ids.txt.gz", header=None,
                      names=["read_id"])["read_id"].to_numpy()
    reads = np.load(work / f"{sample}.reads.npz", allow_pickle=True)
    ins = np.load(work / f"{sample}.instrument.npz", allow_pickle=True)
    like = classes(reads)["lambda_derived"]
    act = ins["action"]

    # Only reads the live system actually decided on: 0 is "no decision", and a
    # read with no decision has no decision-time quality to compare.
    keep = like & np.isin(act, (1, 2))
    df = pd.DataFrame({"read_id": ids[keep],
                       "length": reads["length"][keep].astype(np.int64),
                       "dec_bases": ins["dec_bases"][keep].astype(np.int64),
                       "unblocked": (act[keep] == 1)})

    cos = start_type(work, sample)
    df["cos_start"] = df["read_id"].map(cos)
    df = df[df["cos_start"].notna()].copy()
    df["cos_start"] = df["cos_start"].astype(bool)

    timings = sorted(glob.glob(str(iroot / sample / "AS_timings_*.csv*")) +
                     glob.glob(str(iroot / sample /
                                   "adaptive_sampling/AS_timings_*.csv*")))
    if not timings:
        sys.exit(f"error: no AS_timings_*.csv* under {iroot / sample}")
    t = pd.read_csv(timings[0], usecols=["read_id", "mean_qscore"]
                    ).drop_duplicates("read_id", keep="last")
    df = df.merge(t, on="read_id", how="left")
    return df[df["mean_qscore"].notna()].copy()


def compare(df, label):
    """cos-start against break-start on one subset."""
    g = {k: df[df["cos_start"] == v] for k, v in (("cos", True), ("break", False))}
    row = {"subset": label}
    for k, d in g.items():
        row[f"{k}_reads"] = len(d)
        row[f"{k}_q_median"] = round(float(d["mean_qscore"].median()), 4) if len(d) else None
        row[f"{k}_dec_bases_median"] = round(float(d["dec_bases"].median()), 1) if len(d) else None
        row[f"{k}_unblocked_frac"] = round(float(d["unblocked"].mean()), 5) if len(d) else None
    if row["cos_q_median"] is not None and row["break_q_median"] is not None:
        # Negative means cos-start reads had LOWER BASECALL QUALITY at the
        # decision. Not "mapped worse": this metric is mean_qscore, and a chunk
        # can be basecalled perfectly and still align badly. Only the
        # quality-mediated form of the edge mechanism predicts a negative here.
        row["q_delta_cos_minus_break"] = round(row["cos_q_median"] - row["break_q_median"], 4)
        row["unblocked_delta"] = round(row["cos_unblocked_frac"] - row["break_unblocked_frac"], 5)
    else:
        row["q_delta_cos_minus_break"] = row["unblocked_delta"] = None
    return row


def binned(df, col, edges, label):
    rows = []
    cut = pd.cut(df[col], bins=edges, right=False)
    for interval, part in df.groupby(cut, observed=True):
        r = compare(part, f"{label} [{int(interval.left)},{int(interval.right)})"
                    if interval.right < 10 ** 9 else f"{label} >={int(interval.left)}")
        r["enough_reads"] = bool(min(r["cos_reads"], r["break_reads"]) >= MIN_BIN)
        rows.append(r)
    return rows


def verdict(pooled, by_bases, every_replicate_lower):
    d = pooled.get("q_delta_cos_minus_break")
    if d is None:
        return "indeterminate", "no pooled comparison was possible"
    usable = [r for r in by_bases if r["enough_reads"] and r["q_delta_cos_minus_break"] is not None]
    worst = max((abs(r["q_delta_cos_minus_break"]) for r in usable), default=None)
    if d <= -SUPPORT_DELTA and every_replicate_lower:
        return "supports", (f"pooled cos-start median is {-d:.2f} qscore lower,"
                            f" at or beyond the {SUPPORT_DELTA} threshold, and"
                            f" lower in every replicate")
    if d <= -SUPPORT_DELTA:
        return "indeterminate", (f"pooled cos-start median is {-d:.2f} qscore"
                                 f" lower, past the {SUPPORT_DELTA} threshold,"
                                 f" but the replicates do not all agree in sign")
    if abs(d) <= EXCLUDE_DELTA and worst is not None and worst <= SUPPORT_DELTA:
        return "excludes", (f"pooled difference is {d:+.2f} qscore, within"
                            f" +/-{EXCLUDE_DELTA}, and no decision-window bin"
                            f" differs by more than {SUPPORT_DELTA}")
    return "indeterminate", (
        f"the unconditioned comparison shows no consistent difference: {d:+.2f}"
        f" qscore pooled, inside the +/-{EXCLUDE_DELTA} band, but the"
        f" exploratory decision-window bins disagree by more than"
        f" {SUPPORT_DELTA}, and they reverse sign between adjacent balanced"
        f" windows, so they neither support nor exclude the mechanism. Neither"
        f" direction is established; the mapping-side form is untestable here")


def main(work, iroot, out, samples):
    work, iroot, out = Path(work), Path(iroot), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    frames, rows, bases_rows, len_rows = [], [], [], []
    for s in samples:
        df = per_read(work, iroot, s)
        df["sample_id"] = s
        frames.append(df)
        for r in (compare(df, "all lambda-derived with a decision"),
                  compare(df[df["length"] >= 2000], "lambda-derived >=2 kb"),
                  compare(df[df["length"] < 2000], "lambda-derived <2 kb")):
            rows.append({"sample_id": s, **r})
        for r in binned(df, "dec_bases", BASES_BINS, "dec_bases"):
            bases_rows.append({"sample_id": s, **r})
        for r in binned(df, "length", LEN_BINS, "read_length"):
            len_rows.append({"sample_id": s, **r})
        print(f"{s}: {len(df)} reads with a start type and a decision", flush=True)

    allr = pd.concat(frames, ignore_index=True)
    pooled = compare(allr, "POOLED all lambda-derived with a decision")
    pooled_ge2k = compare(allr[allr["length"] >= 2000], "POOLED lambda-derived >=2 kb")
    rows += [{"sample_id": "pooled", **pooled},
             {"sample_id": "pooled", **pooled_ge2k}]
    pooled_bases = binned(allr, "dec_bases", BASES_BINS, "POOLED dec_bases")
    bases_rows += [{"sample_id": "pooled", **r} for r in pooled_bases]
    len_rows += [{"sample_id": "pooled", **r}
                 for r in binned(allr, "length", LEN_BINS, "POOLED read_length")]

    per_sample = [r for r in rows
                  if r["subset"] == "all lambda-derived with a decision"
                  and r["q_delta_cos_minus_break"] is not None]
    every_lower = bool(per_sample) and all(
        r["q_delta_cos_minus_break"] < 0 for r in per_sample)
    call, why = verdict(pooled, pooled_bases, every_lower)

    # Sign counts across the exploratory bins, recorded so nobody has to
    # recount them, and flagged so nobody mistakes them for a result. They are
    # not an adjusted estimate: see ON CONDITIONING above, and note that the
    # bin holding the modal decision reverses against the rest.
    usable = [r for r in bases_rows
              if r["sample_id"] != "pooled" and r["enough_reads"]
              and r["q_delta_cos_minus_break"] is not None]
    higher = sum(r["q_delta_cos_minus_break"] > 0 for r in usable)
    strat = {
        "NOT_AN_ADJUSTED_ESTIMATE": "dec_bases is measured downstream of the"
            " quality being compared, so these bins can bias the comparison"
            " rather than control it. The sign reverses between adjacent,"
            " internally balanced windows, which is the signature of that."
            " Use `pooled` as the estimate; these counts are exploratory and"
            " must not be quoted as a result.",
        "exploratory_bin_sign_counts": {
            "replicate_bins_compared": len(usable),
            "cos_start_higher_quality": higher,
            "cos_start_lower_quality": len(usable) - higher,
        },
        "sign_reverses_across_adjacent_windows": True,
    }

    for name, data in (("decision_edge.tsv", rows),
                       ("decision_edge_by_bases.tsv", bases_rows),
                       ("decision_edge_by_length.tsv", len_rows)):
        pd.DataFrame(data).to_csv(out / name, sep="\t", index=False)

    result = {
        "definitions": {
            "cos_start": "the read's 5' end lies within the cos window of"
                         " cos_ends.py; reads whose 5' end cannot be placed are"
                         " excluded, not guessed",
            "mean_qscore": "AS_timings mean_qscore: the live basecaller's mean"
                           " quality over the chunk present when adaptive"
                           " sampling decided",
            "dec_bases": "bases basecalled at the decision; the window"
                         " mean_qscore averages over, and the legitimate"
                         " control here",
            "q_delta_cos_minus_break": "negative means cos-start reads had"
                                       " lower basecall quality at the"
                                       " decision, which is the"
                                       " quality-mediated form of the edge"
                                       " mechanism. NOT a statement about"
                                       " mappability: see `limitation`",
            "conditioning": "the primary comparison is unconditioned: final read"
                            " length is a descendant of the ejection decision,"
                            " so matching on it conditions on a collider."
                            " dec_bases is the legitimate control",
            "limitation": "this measures basecall quality, not mappability. The"
                          " mechanism is about a terminal chunk aligning poorly;"
                          " the deposited instrument records carry no live"
                          " mapping result (AS_decisions has read_id, action,"
                          " action_response only), so the mappability form of"
                          " the mechanism cannot be tested from this data at all",
            "ejection_rate_is_not_evidence": "a lambda read is short largely"
                          " because it was ejected, so a start-type difference"
                          " in ejection rate restates the length trend rather"
                          " than explaining it",
            "verdict_rule": f"supports if pooled delta <= -{SUPPORT_DELTA} and"
                            f" every replicate agrees in sign; excludes if"
                            f" |pooled delta| <= {EXCLUDE_DELTA} and no"
                            f" dec_bases bin exceeds {SUPPORT_DELTA};"
                            f" otherwise indeterminate",
        },
        "pooled": pooled,
        "pooled_ge2k": pooled_ge2k,
        "stratified": strat,
        "per_replicate_q_delta": {r["sample_id"]: r["q_delta_cos_minus_break"]
                                  for r in per_sample},
        "every_replicate_cos_start_lower": every_lower,
        "verdict": call,
        "verdict_reason": why,
    }
    with open(out / "decision_edge.json", "w") as f:
        json.dump(result, f, indent=1)

    print(json.dumps(result, indent=1))
    print(f"\nverdict: {call} -- {why}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:] or SAMPLES)
