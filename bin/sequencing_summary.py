#!/usr/bin/env python3
"""
Per-replicate sequencing summary: yield, read length and read quality, before
and after the carrier and contaminant are removed.

The pipeline reports per-organism counts and a pooled read-length distribution,
but nothing that answers the first question anyone asks of a sequencing run --
how many reads, how many bases, how long, how good -- and nothing at all about
read quality, because quality never enters the assignment. It lives only in the
FASTQ header.

This joins the two: read length and role come from
`<sample>.readlengths.tsv.gz`, and the mean read qscore comes from the `qs:f:`
tag of the corresponding FASTQ record.

The join is positional. minimap2 emits alignments in input order and
`assign_reads.py` streams the resulting BAM in that order, so the Nth row of
readlengths corresponds to the Nth FASTQ record. That is asserted on every read
rather than assumed; a mismatch aborts rather than silently pairing the wrong
quality with the wrong role. Positional joining is what makes this feasible at
all -- a read_id-keyed dictionary over ~10 million reads costs gigabytes, and
this needs none.

On the qscore: ONT reports the mean read quality in error-probability space,
`-10*log10(mean(10^(-Q/10)))`, not as an arithmetic mean of Phred values. The
`qs:f:` tag already carries that number, so it is used as given. Medians of it
are order statistics and need no such care.

Three aggregates are reported per replicate:

  all               every read in the FASTQ. Note this is ALREADY the output of
                    depletion-mode adaptive sampling -- carrier reads were
                    rejected on the instrument during the run -- so there is no
                    pre-adaptive-sampling read set to compare against.
  carrier_removed   what is left after computationally removing everything
                    attributed to the CARRIER, INCLUDING
                    ambiguous classes tied only among carrier organisms. This
                    is the in-silico analogue of what adaptive sampling did on
                    the instrument, and is the like-for-like comparison.
  carrier_contaminant_removed
                    as above but also removing the contaminant. This is what is
                    left once every organism traceable to the carrier
                    preparation is gone.

  The subsets are deliberately NOT named "depleted": depletion-mode adaptive
  sampling is a real-time instrument behaviour that targeted the carrier alone,
  and every read here is already its output. There is no pre-adaptive-sampling
  read set to compare against.
  contaminant       the carrier-derived contaminant on its own, so its read
                    length and quality can be compared with the community's
  community         reads assigned to a community organism (role=sample)

Usage:
    python3 bin/sequencing_summary.py --out results/summary/sequencing_summary.tsv
    python3 bin/sequencing_summary.py lowinput_s1_r1 --results results
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("replicates", nargs="*", help="default: every finished replicate")
    p.add_argument("--results", default="results")
    p.add_argument("--samplesheets", default="assets/samplesheets")
    p.add_argument("--mode", default="competitive")
    p.add_argument("--out", default=None, help="TSV output path")
    return p.parse_args()


def fastq_path(sample_id, sheet_dir):
    """Resolve the FASTQ from the samplesheets rather than guessing a name."""
    for p in sorted(Path(sheet_dir).glob("*.csv")):
        for line in p.read_text().splitlines():
            if line.startswith(sample_id + ","):
                f = line.split(",")
                if len(f) > 3 and f[3]:
                    return Path(f[3])
    return Path("data") / f"{sample_id}.fastq"


def qs_from_header(line):
    """Mean read qscore from the ONT `qs:f:` tag."""
    for field in line.rstrip("\n").split("\t")[1:]:
        if field.startswith("qs:f:"):
            return float(field[5:])
    return float("nan")


def collect(sample_id, results, sheet_dir, mode):
    rl_path = Path(results) / sample_id / mode / f"{sample_id}.readlengths.tsv.gz"
    fq_path = fastq_path(sample_id, sheet_dir)
    if not rl_path.is_file():
        sys.exit(f"error: {rl_path} not found")
    if not fq_path.is_file():
        sys.exit(f"error: {fq_path} not found")

    lengths, qscores, roles, orgs = [], [], [], []
    opener = gzip.open if str(fq_path).endswith(".gz") else open

    with gzip.open(rl_path, "rt") as rl, opener(fq_path, "rt") as fq:
        rl.readline()                                   # header
        for n, rline in enumerate(rl):
            head = fq.readline()
            if not head:
                sys.exit(f"error: {fq_path} ran out of reads at row {n+1}; the "
                         f"FASTQ and the assignments do not describe the same run")
            fq.readline(); fq.readline(); fq.readline()  # seq, +, qual

            rid, _call, org, role, length = rline.rstrip("\n").split("\t")
            fq_id = head[1:].split(None, 1)[0].split("\t", 1)[0]
            if fq_id != rid:
                sys.exit(f"error: positional join broke at row {n+1}: FASTQ has "
                         f"{fq_id}, assignments have {rid}. The two files are not "
                         f"in the same order; this script cannot be used as-is.")
            lengths.append(int(length)); qscores.append(qs_from_header(head))
            roles.append(role); orgs.append(org)

    return (np.array(lengths, dtype=np.int64), np.array(qscores, dtype=np.float64),
            np.array(roles), np.array(orgs))


def attributable(roles, orgs, drop_roles, drop_orgs):
    """Reads attributable to the given roles, ties included.

    An `ambiguous:A|B` read whose tied organisms all belong to the roles being
    dropped came from one of them whichever way the tie falls, so it is
    attributable regardless. Keeping it would flatter whatever remains.
    """
    mask = np.isin(roles, drop_roles)
    for i in np.flatnonzero(roles == "ambiguous"):
        tied = orgs[i].split(":", 1)[1].split("|") if ":" in orgs[i] else []
        if tied and all(t in drop_orgs for t in tied):
            mask[i] = True
    return mask


def main():
    args = parse_args()
    root = Path(__file__).resolve().parent.parent
    results = (root / args.results) if not Path(args.results).is_absolute() else Path(args.results)
    sheets = (root / args.samplesheets) if not Path(args.samplesheets).is_absolute() else Path(args.samplesheets)

    reps = args.replicates or sorted(
        d.name for d in results.iterdir()
        if (d / args.mode / f"{d.name}.readlengths.tsv.gz").is_file())
    if not reps:
        sys.exit(f"error: no finished replicates under {results}")

    rows = []
    for sid in reps:
        # Carrier/contaminant organism names come from the contig map for this
        # replicate's reference set, so nothing is hardcoded to lambda/E. coli.
        set_name = None
        for p in sorted(sheets.glob("*.csv")):
            for line in p.read_text().splitlines():
                if line.startswith(sid + ","):
                    f = line.split(",")
                    if len(f) > 5 and f[5]:
                        set_name = Path(f[5]).stem
        cmap = results / "references" / (set_name or sid.rsplit("_r", 1)[0]) / "contig_map.tsv"
        carrier_orgs, contam_orgs = set(), set()
        if cmap.is_file():
            for line in cmap.read_text().splitlines()[1:]:
                f = line.split("\t")
                if len(f) >= 3 and f[2] == "carrier":
                    carrier_orgs.add(f[1])
                elif len(f) >= 3 and f[2] == "contaminant":
                    contam_orgs.add(f[1])

        print(f"[sequencing-summary] {sid}: reading", file=sys.stderr)
        lengths, qs, roles, orgs = collect(sid, results, sheets, args.mode)

        drop_carrier = attributable(roles, orgs, ("carrier",), carrier_orgs)
        drop_both = attributable(roles, orgs, ("carrier", "contaminant"),
                                 carrier_orgs | contam_orgs)
        groups = {
            "all":                          np.ones(len(lengths), dtype=bool),
            "carrier_removed":              ~drop_carrier,
            "carrier_contaminant_removed":  ~drop_both,
            "contaminant":                  roles == "contaminant",
            "community":                    roles == "sample",
        }
        for name, m in groups.items():
            if not m.any():
                continue
            L, Q = lengths[m], qs[m]
            good = ~np.isnan(Q)
            rows.append({
                "sample_id": sid, "subset": name,
                "reads": int(m.sum()), "bases": int(L.sum()),
                "median_read_length": int(np.median(L)),
                "mean_read_length": round(float(L.mean()), 1),
                "median_qscore": round(float(np.median(Q[good])), 2) if good.any() else "",
            })
        print(f"  reads {len(lengths):,}  -carrier {int((~drop_carrier).sum()):,}  "
              f"-carrier-contam {int((~drop_both).sum()):,}  "
              f"community {int((roles=='sample').sum()):,}", file=sys.stderr)

    cols = ["sample_id", "subset", "reads", "bases", "median_read_length",
            "mean_read_length", "median_qscore"]
    out = []
    out.append("\t".join(cols))
    for r in rows:
        out.append("\t".join(str(r[c]) for c in cols))
    text = "\n".join(out) + "\n"

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text)
        print(f"[sequencing-summary] -> {args.out}", file=sys.stderr)
    print(text)


if __name__ == "__main__":
    main()
