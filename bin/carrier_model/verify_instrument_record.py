#!/usr/bin/env python3
"""Verify the instrument-data record and write its checksum manifest.

For each replicate, run join_instrument.py on the record's (trimmed) files and
require every output array to equal the one built from MinKNOW's original files
(<workdir>/<sample>.instrument.npz). That proves the record loses nothing the
analysis uses. Then write the manifest: file, bytes, rows, md5, sha256, and the
run's flow cell, MinKNOW version and SRA run accession. The manifest is
committed as assets/instrument_files.tsv.

Usage: verify_instrument_record.py <record_dir> <results> <workdir> <repo> <manifest_out>
"""
import csv
import gzip
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import SAMPLES, read_tsv  # noqa: E402


def digests(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            md5.update(chunk); sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def rows(path):
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) - 1


def main(record, results, work, repo, manifest_out):
    record, results, work, repo = map(Path, (record, results, work, repo))
    meas = {r["sample_id"]: r for r in read_tsv(repo / "assets/measurements.tsv")}
    sra = {}
    with open(repo / "assets/samplesheets/all.csv") as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            sra[r["sample_id"]] = r.get("sra_accession", "")
    manifest, ok = [], True
    for s in SAMPLES:
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "join_instrument.py"), s,
                            str(results / s / "competitive" / f"{s}.readlengths.tsv.gz"),
                            str(work / f"{s}.read_ids.txt.gz"), str(record / s), tmp], check=True, capture_output=True)
            a = np.load(work / f"{s}.instrument.npz")
            b = np.load(Path(tmp) / f"{s}.instrument.npz")
            bad = [k for k in a.files if not np.array_equal(a[k], b[k], equal_nan=a[k].dtype.kind == "f")]
        status = "identical" if not bad else f"DIFFERS: {bad}"
        ok &= not bad
        print(f"{s}: join from the record's files {status}", flush=True)
        for f in sorted((record / s).glob("*.gz")):
            md5, sha = digests(f)
            manifest.append({"sample_id": s, "file": f"{s}/{f.name}", "bytes": f.stat().st_size,
                             "rows": rows(f), "md5": md5, "sha256": sha,
                             "flow_cell_id": meas[s]["flow_cell_id"],
                             "minknow_version": meas[s]["minknow_version"],
                             "sra_run": sra.get(s, ""), "rejoin_check": status})
    with open(manifest_out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(manifest)
    print("ALL IDENTICAL" if ok else "SOME DIFFER")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(*sys.argv[1:])
