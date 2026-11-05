#!/usr/bin/env python3
"""
Count Kraken2 hits per sample from wf-metagenomics per-read assignments
(<out>/reads_assignments/<alias>_lineages.kraken2.assignments.tsv), using the
rule behind each study's `kraken2_*` rows in assets/comparison/prior_studies.tsv.

  zorzano      a hit is any read classified to any taxon (column 1 == "C"); no
               organism filter, matching how Zorzano et al. report hits.
  raghavendra  a hit for a target species is a read whose lineage (column 6)
               contains that species' name. Strain and sub-species assignments
               therefore count, and genus-only or higher assignments do not,
               matching the paper's species-rank reporting. Targets:
               Escherichia coli, Saccharomyces cerevisiae. Reads whose
               "<alias>__" tag does not name the sample are skipped and counted.

Columns of the assignments file (Kraken2 output with EPI2ME's lineage
appended): C/U, read id, taxid, length, k-mer hits, lineage.

This is the scripted form of the counting done for the EPI2ME runs of 29 April
2026 and applies the same rules.

Usage:
    count_kraken2_hits.py zorzano     --run DIR [--label q1]
    count_kraken2_hits.py raghavendra --run DIR [--label q10]

Writes a TSV to stdout:
    study label sample target reads_in bases_in hit_reads hit_bases
"""
import argparse
import sys
from pathlib import Path

ZORZANO = ["Blank", "Microbialite", "Microb_086MGy", "Microb_1045MGy",
           "Stromat", "Stromat_1045MGy", "OxideIron", "OxideIron_1045MGy",
           "Carbonate"]
RAGHAVENDRA = ["Ec_R1", "Ec_R2", "Ec_R3", "YSC_R1", "YSC_R2", "YSC_R3",
               "Mix1_R1", "Mix1_R2", "Mix2_R1", "Mix2_R2"]
TARGETS = ["Escherichia coli", "Saccharomyces cerevisiae"]


def records(path: Path):
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            try:
                length = int(p[3])
            except ValueError:
                continue
            yield p[0], p[1], length, (p[5] if len(p) >= 6 else "")


def assignments(run: Path, sample: str) -> Path:
    for base in (run / "reads_assignments", run / "output" / "reads_assignments"):
        f = base / f"{sample}_lineages.kraken2.assignments.tsv"
        if f.is_file():
            return f
    sys.exit(f"error: no assignments for {sample} under {run}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("study", choices=["zorzano", "raghavendra"])
    ap.add_argument("--run", required=True, type=Path,
                    help="wf-metagenomics output directory")
    ap.add_argument("--label", default="", help="free text, e.g. q1 or q10")
    a = ap.parse_args()

    print("study\tlabel\tsample\ttarget\treads_in\tbases_in\thit_reads\thit_bases")
    if a.study == "zorzano":
        for s in ZORZANO:
            n = b = hn = hb = 0
            for cu, _rid, length, _lin in records(assignments(a.run, s)):
                n += 1; b += length
                if cu == "C":
                    hn += 1; hb += length
            print(f"zorzano\t{a.label}\t{s}\tany\t{n}\t{b}\t{hn}\t{hb}")
        return

    mistagged = 0
    for s in RAGHAVENDRA:
        n = b = 0
        hits = {t: [0, 0] for t in TARGETS}
        for _cu, rid, length, lineage in records(assignments(a.run, s)):
            if rid.split("__", 1)[0] != s:
                mistagged += 1
                continue
            n += 1; b += length
            for t in TARGETS:
                if t in lineage:
                    hits[t][0] += 1; hits[t][1] += length
        for t in TARGETS:
            print(f"raghavendra\t{a.label}\t{s}\t{t}\t{n}\t{b}\t{hits[t][0]}\t{hits[t][1]}")
    if mistagged:
        print(f"warning: {mistagged} reads had a sample tag that did not match "
              f"their file and were skipped", file=sys.stderr)


if __name__ == "__main__":
    main()
