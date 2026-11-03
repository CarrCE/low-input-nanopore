#!/usr/bin/env python3
"""Build the instrument-data record's files for one replicate, from MinKNOW output.

This is how the published instrument records were made; a reader needs it only
to audit them. It reads the raw MinKNOW run folder, which is not public.

Rows are filtered as TEXT, so every retained field is byte-identical to the
instrument output: only rows for reads in the SRA deposit are kept (the
deposited, human-masked read set is the reference), and the summary loses two
columns that carry no information about the reads (`filename`, a local POD5
path, and `barcode`, "unclassified" throughout). Line endings are normalised to
LF (the decision files were written with CRLF). No file holds sequence.

Outputs in <outdir>/<sample_id>/ (gzip, mtime 0, so checksums are reproducible):
  summary_<sample_id>.tsv.gz        Dorado sequencing summary (sup basecalling)
  AS_decisions_<sample_id>.csv.gz   read_id, action, action_response
  AS_timings_<sample_id>.csv.gz     per-read adaptive-sampling decision timing
  pore_activity_<sample_id>.csv.gz  channel-state time per minute (unfiltered)
and <outdir>/<sample_id>/build.json with row counts in and out.

Usage: build_instrument_record.py <sample_id> <readlengths.tsv.gz> <minknow_dir> <outdir>
"""
import glob, gzip, io, json, sys
from pathlib import Path

DROP_SUMMARY = {"filename", "barcode"}


def one(patterns):
    hits = sorted(set(h for p in patterns for h in glob.glob(p)))
    if len(hits) != 1:
        sys.exit(f"expected one file for {patterns}, found {hits}")
    return hits[0]


def gz_writer(path):
    raw = open(path, "wb")
    return raw, io.TextIOWrapper(gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0),
                                 encoding="utf-8", newline="\n")


def filter_file(src, dst, keep_ids, id_col, sep, drop=()):
    raw, out = gz_writer(dst)
    n_in = n_out = 0
    with open(src, encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n").split(sep)
        idx = header.index(id_col)
        cols = [i for i, h in enumerate(header) if h not in drop]
        out.write(sep.join(header[i] for i in cols) + "\n")
        for line in fh:
            f = line.rstrip("\r\n").split(sep)
            n_in += 1
            if f[idx] in keep_ids:
                out.write(sep.join(f[i] for i in cols) + "\n")
                n_out += 1
    out.close(); raw.close()
    return {"source": Path(src).name, "rows_in": n_in, "rows_out": n_out,
            "columns_dropped": sorted(set(header) & set(drop))}


def main(sample, rl_path, idir, outdir):
    idir = Path(idir); od = Path(outdir) / sample
    od.mkdir(parents=True, exist_ok=True)
    with gzip.open(rl_path, "rt") as fh:
        next(fh)
        keep = {line.split("\t", 1)[0] for line in fh}
    info = {"sample_id": sample, "deposit_reads": len(keep)}
    info["summary"] = filter_file(one([str(idir / "summary*.tsv")]),
                                  od / f"summary_{sample}.tsv.gz", keep, "read_id", "\t", DROP_SUMMARY)
    info["AS_decisions"] = filter_file(
        one([str(idir / "AS_decisions_*.csv"), str(idir / "adaptive_sampling/AS_decisions_*.csv")]),
        od / f"AS_decisions_{sample}.csv.gz", keep, "read_id", ",")
    info["AS_timings"] = filter_file(
        one([str(idir / "AS_timings_*.csv"), str(idir / "adaptive_sampling/AS_timings_*.csv")]),
        od / f"AS_timings_{sample}.csv.gz", keep, "read_id", ",")
    src = one([str(idir / "pore_activity_*.csv")])
    raw, out = gz_writer(od / f"pore_activity_{sample}.csv.gz")
    with open(src, encoding="utf-8", newline="") as fh:
        n = 0
        for line in fh:
            out.write(line.rstrip("\r\n") + "\n"); n += 1
    out.close(); raw.close()
    info["pore_activity"] = {"source": Path(src).name, "rows_in": n - 1, "rows_out": n - 1,
                             "columns_dropped": []}
    json.dump(info, open(od / "build.json", "w"), indent=1)
    print(json.dumps(info))


if __name__ == "__main__":
    main(*sys.argv[1:])
