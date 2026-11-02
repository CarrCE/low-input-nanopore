#!/usr/bin/env python3
"""Join MinKNOW per-read metadata onto the DEPOSITED reads.

The deposited (human-masked) reads define the read set: every output array is
aligned to the order of <s>.readlengths.tsv.gz, which is the deposited FASTQ
order (extract_reads.py asserts this). Reads are keyed by their ONT read id,
from <workdir>/<s>.read_ids.txt.gz (extract_reads.py), so the join works the
same whether the pipeline named reads by ONT id or by SRA spot. The instrument files only contribute
metadata keyed on read_id; rows for reads that are not in the deposit are
counted and dropped. The instrument files hold no sequence.

Inputs, in <instrument_dir> (one replicate's folder of the instrument-data
record; see docs/carrier-model.md). Plain or gzipped:
  summary*.tsv                       Dorado sequencing summary
  AS_decisions_*.csv                 read_id, action, action_response
  AS_timings_*.csv                   read_id, batch_time, samples, bases, ...
(an adaptive_sampling/ subfolder, as MinKNOW writes it, is also accepted)

Output: <workdir>/<sample>.instrument.npz with, aligned to deposit order:
  in_summary, channel, mux, start_time (s from run start), duration (s),
  action (0 none, 1 unblock, 2 sequence), dec_bases (bases seen at decision, -1 none),
  dec_time (batch_time, epoch s; nan none), t_basecaller (as recorded; nan none)
and <workdir>/<sample>.instrument.json with join statistics.

Usage: join_instrument.py <sample_id> <readlengths.tsv.gz> <read_ids.txt.gz> <instrument_dir> <workdir>
"""
import glob, json, sys
from pathlib import Path
import numpy as np
import pandas as pd


def one(pattern_list):
    hits = sorted(set(h for p in pattern_list for h in glob.glob(p)))
    if len(hits) != 1:
        sys.exit(f"expected exactly one file for {pattern_list}, found {hits}")
    return hits[0]


def pats(idir, stem, ext):
    return [str(idir / sub / f"{stem}.{ext}{gz}") for sub in ("", "adaptive_sampling")
            for gz in ("", ".gz")]


def main(sample, rl_path, ids_path, idir, workdir):
    idir, work = Path(idir), Path(workdir)
    stats = {"sample_id": sample}
    dep = pd.read_csv(rl_path, sep="\t", usecols=["read_length"])
    ids = pd.read_csv(ids_path, header=None, names=["read_id"])["read_id"]
    assert len(ids) == len(dep), f"{sample}: {len(ids)} read ids for {len(dep)} reads"
    dep.insert(0, "read_id", ids.to_numpy())
    stats["deposit_reads"] = len(dep)

    summ = pd.read_csv(one(pats(idir, "summary*", "tsv")), sep="\t",
                       usecols=["read_id", "channel", "mux", "start_time", "duration",
                                "sequence_length_template"],
                       float_precision="round_trip")
    stats["summary_rows"] = len(summ)
    stats["summary_duplicate_ids"] = int(summ["read_id"].duplicated().sum())
    summ = summ.drop_duplicates("read_id", keep="first")
    m = dep.merge(summ, on="read_id", how="left", validate="one_to_one")
    in_s = m["channel"].notna().to_numpy()
    stats["deposit_reads_in_summary"] = int(in_s.sum())
    stats["summary_rows_not_in_deposit"] = int(len(summ) - in_s.sum())
    mism = in_s & (m["sequence_length_template"].fillna(-1).to_numpy() != m["read_length"].to_numpy())
    stats["length_mismatches"] = int(mism.sum())

    dec = pd.read_csv(one(pats(idir, "AS_decisions_*", "csv")))
    dec["action_response"] = dec["action_response"].astype(str).str.strip()
    stats["decision_rows"] = len(dec)
    stats["decision_responses"] = dec["action_response"].value_counts().to_dict()
    stats["decision_duplicate_ids"] = int(dec["read_id"].duplicated().sum())
    dec = dec.drop_duplicates("read_id", keep="last")
    code = {"unblock": 1, "sequence": 2}
    dec["action_code"] = dec["action"].map(code).fillna(3).astype(np.int8)
    m = m.merge(dec[["read_id", "action_code"]], on="read_id", how="left")
    stats["decision_rows_not_in_deposit"] = int(len(dec) - m["action_code"].notna().sum())

    tim = pd.read_csv(one(pats(idir, "AS_timings_*", "csv")),
                      usecols=["read_id", "batch_time", "bases", "time_in_basecaller"],
                      float_precision="round_trip")
    stats["timing_rows"] = len(tim)
    tim = tim.drop_duplicates("read_id", keep="last")
    m = m.merge(tim, on="read_id", how="left")
    assert len(m) == len(dep) and (m["read_id"].to_numpy() == dep["read_id"].to_numpy()).all(), \
        "join changed the deposit order"

    np.savez_compressed(
        work / f"{sample}.instrument.npz",
        in_summary=in_s,
        channel=m["channel"].fillna(-1).astype(np.int16).to_numpy(),
        mux=m["mux"].fillna(-1).astype(np.int8).to_numpy(),
        start_time=m["start_time"].to_numpy(dtype=np.float64),
        duration=m["duration"].to_numpy(dtype=np.float32),
        action=m["action_code"].fillna(0).astype(np.int8).to_numpy(),
        dec_bases=m["bases"].fillna(-1).astype(np.int32).to_numpy(),
        dec_time=m["batch_time"].to_numpy(dtype=np.float64),
        t_basecaller=m["time_in_basecaller"].to_numpy(dtype=np.float64))
    vc = m["action_code"].fillna(0).value_counts().to_dict()
    stats["deposit_action_counts"] = {{0: "none", 1: "unblock", 2: "sequence", 3: "other"}[int(k)]: int(v)
                                      for k, v in vc.items()}
    json.dump(stats, open(work / f"{sample}.instrument.json", "w"), indent=1)
    print(json.dumps(stats))


if __name__ == "__main__":
    main(*sys.argv[1:])
