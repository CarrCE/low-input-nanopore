#!/usr/bin/env python3
"""Pass 1: one stream over a deposited (human-masked) FASTQ, in lockstep with its
competitive read-length table, producing

  <workdir>/<sample>.reads.npz        per-read arrays in deposit order: length,
                                      start time, role, organism, qscore
  <workdir>/<sample>.lambda_ends.fa   lambda-derived reads for end mapping:
                                      every one >= LONG_BP, plus a deterministic
                                      hash sample (fraction P) of all of them.
                                      Reads longer than 2*END_BP are written as
                                      two records, <id>_head and <id>_tail (first
                                      and last END_BP bases), which is all end
                                      mapping needs.
  <workdir>/<sample>.lambda_ids.tsv   read_id, length, sampled_by (hash|long),
                                      call organism

The FASTQ and readlengths.tsv.gz are in the same order (bin/sequencing_summary.py
relies on the same fact); every record's read_id and length are asserted equal,
so a misalignment fails loudly rather than mislabelling reads. Human masking
replaced bases with N and kept read count and length, so the masked deposit is a
complete input: nothing here needs unmasked reads.

Usage: extract_reads.py <sample_id> <fastq[.gz]> <readlengths.tsv.gz> <workdir>
"""
import gzip
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROLES, is_lambda_derived  # noqa: E402

P = 0.02
LONG_BP = 2000
END_BP = 1500


def hashed(rid):
    return int(rid[:8], 16) / 2**32 < P


def main(sample, fastq, rl_path, workdir):
    work = Path(workdir)
    cmd = ["gzip", "-dc", fastq] if fastq.endswith(".gz") else ["cat", fastq]
    fq = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 24)
    rl = gzip.open(rl_path, "rt")
    header = rl.readline().rstrip("\n").split("\t")
    assert header == ["read_id", "call", "organism", "role", "read_length"], header

    orgs = {}
    L, ST, RO, OR, QS = [], [], [], [], []
    out_fa = open(work / f"{sample}.lambda_ends.fa", "w")
    out_ids = open(work / f"{sample}.lambda_ids.tsv", "w")
    out_ids.write("read_id\tread_length\tsampled_by\tcall_organism\n")
    n = 0
    it = iter(fq.stdout)
    for h in it:
        seq = next(it); next(it); next(it)
        fields = h.decode()[1:].rstrip("\n").split("\t")
        rid = fields[0]
        tags = dict(f.split(":", 2)[0::2] for f in fields[1:] if f.count(":") >= 2)
        row = rl.readline().rstrip("\n").split("\t")
        if row[0] != rid:
            sys.exit(f"{sample}: order mismatch at record {n}: fastq {rid} vs readlengths {row[0]}")
        length = int(row[4])
        if length != len(seq) - 1:
            sys.exit(f"{sample}: length mismatch for {rid}: table {length} vs fastq {len(seq)-1}")
        st = datetime.fromisoformat(tags["st"]).timestamp() if "st" in tags else np.nan
        role, org = row[3], row[2]
        L.append(length); ST.append(st); RO.append(ROLES.index(role))
        OR.append(orgs.setdefault(org, len(orgs)))
        QS.append(float(tags.get("qs", "nan")))
        if is_lambda_derived(role, org):
            by = "long" if length >= LONG_BP else ("hash" if hashed(rid) else None)
            if by:
                sq = seq.decode().rstrip("\n")
                if length > 2 * END_BP:
                    out_fa.write(f">{rid}_head\n{sq[:END_BP]}\n>{rid}_tail\n{sq[-END_BP:]}\n")
                else:
                    out_fa.write(f">{rid}_whole\n{sq}\n")
                out_ids.write(f"{rid}\t{length}\t{by}\t{org}\n")
        n += 1
    assert rl.readline() == "", f"{sample}: readlengths has rows beyond the FASTQ"
    if fq.wait() != 0:
        sys.exit(f"{sample}: decompression failed")
    out_fa.close(); out_ids.close()
    np.savez_compressed(
        work / f"{sample}.reads.npz",
        length=np.array(L, dtype=np.int32), st=np.array(ST, dtype=np.float64),
        role=np.array(RO, dtype=np.int8), org=np.array(OR, dtype=np.int16),
        qs=np.array(QS, dtype=np.float32),
        roles=np.array(ROLES), orgs=np.array(list(orgs)), p_hash=P, long_bp=LONG_BP)
    print(f"[carrier-model] {sample}: {n:,} reads")


if __name__ == "__main__":
    main(*sys.argv[1:])
