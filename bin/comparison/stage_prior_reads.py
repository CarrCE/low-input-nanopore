#!/usr/bin/env python3
"""
Lay out a prior study's reads the way wf-metagenomics expects them for the
Kraken2 reanalysis: one directory per sample, named by the sample alias that
the rows in assets/comparison/prior_studies.tsv use.

This is the scripted form of the staging done by hand before the EPI2ME runs of
29 April 2026, and it reproduces those inputs exactly:

  zorzano      copies each deposited FASTQ unchanged to <out>/<alias>/<alias>.fastq.gz.
               Input: RawData/ of FISABIO-BioinformaticsService/
               CommsEarth_Environment-fDNArocks (bin/comparison/fetch_zorzano.sh).

  raghavendra  rewrites each read name as <alias>__<original id>, keeping the
               rest of the header, so a read's sample survives any pooling.
               Input: the per-replicate files written by
               bin/comparison/fetch_raghavendra.sh (Ec_R1.fastq.gz, ...).

Usage:
    stage_prior_reads.py zorzano     --src DIR --out DIR
    stage_prior_reads.py raghavendra --src DIR --out DIR
"""
import argparse
import gzip
import shutil
import sys
from pathlib import Path

ZORZANO = {
    # alias            deposited file in RawData/
    "Blank":             "ExtractedWaterControl.fastq.gz",
    "Microbialite":      "MicrobialitesNonIrradiated.fastq.gz",
    "Microb_086MGy":     "MicrobialitesIrradiated0.874MGy.fastq.gz",
    "Microb_1045MGy":    "MicrobialitesIrradiated10.5MGy.fastq.gz",
    "Stromat":           "StromatolitesNonIrradiated.fastq.gz",
    "Stromat_1045MGy":   "StromatolitesIrradiated10.5MGy.fastq.gz",
    "OxideIron":         "OxideStromatolitesNonIrradiated.fastq.gz",
    "OxideIron_1045MGy": "OxideStromatolitesIrradiated.fastq.gz",
    "Carbonate":         "Carbonates.fastq.gz",
}

RAGHAVENDRA = ["Ec_R1", "Ec_R2", "Ec_R3", "YSC_R1", "YSC_R2", "YSC_R3",
               "Mix1_R1", "Mix1_R2", "Mix2_R1", "Mix2_R2"]


def stage_zorzano(src: Path, out: Path) -> None:
    for alias, name in ZORZANO.items():
        f = src / name
        if not f.is_file():
            sys.exit(f"error: {f} not found")
        (out / alias).mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f, out / alias / f"{alias}.fastq.gz")
        print(f"{alias:<18} <- {name}")


def stage_raghavendra(src: Path, out: Path) -> None:
    for alias in RAGHAVENDRA:
        f = src / f"{alias}.fastq.gz"
        if not f.is_file():
            sys.exit(f"error: {f} not found")
        (out / alias).mkdir(parents=True, exist_ok=True)
        n = 0
        with gzip.open(f, "rt") as fin, \
                gzip.open(out / alias / f"{alias}.fastq.gz", "wt") as fout:
            while True:
                header = fin.readline()
                if not header:
                    break
                seq, plus, qual = fin.readline(), fin.readline(), fin.readline()
                if not qual:
                    sys.exit(f"error: truncated record at the end of {f}")
                if not header.startswith("@"):
                    sys.exit(f"error: bad FASTQ header in {f}: {header[:80]!r}")
                rid, sep, rest = header[1:].rstrip("\n").partition(" ")
                rest = rest.lstrip()
                fout.write(f"@{alias}__{rid}{(' ' + rest) if rest else ''}\n")
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)
                n += 1
        print(f"{alias:<8} {n:>6} reads")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("study", choices=["zorzano", "raghavendra"])
    ap.add_argument("--src", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    {"zorzano": stage_zorzano, "raghavendra": stage_raghavendra}[a.study](a.src, a.out)


if __name__ == "__main__":
    main()
