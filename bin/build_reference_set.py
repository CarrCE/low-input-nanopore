#!/usr/bin/env python3
"""
Assemble one combined reference from a reference-set TSV plus genomes fetched
by `datasets download genome accession`.

Emits the single FASTA that competitive assignment maps against, together with
the contig -> organism -> role map that turns per-contig alignments back into
per-organism counts, and a genome-size table used later for coverage breadth.

Contig names are left untouched so alignments remain traceable to the exact
NCBI record; the mapping to organism is carried alongside rather than by
renaming.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference-tsv", required=True,
                   help="assets/references/<set>.tsv")
    p.add_argument("--genome-dir", required=True,
                   help="directory containing ncbi_dataset/data/<accession>/*.fna")
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-contig-map", required=True)
    p.add_argument("--out-genome-sizes", required=True)
    p.add_argument("--out-provenance", default=None)
    return p.parse_args()


def read_reference_tsv(path):
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    if not rows:
        sys.exit(f"error: no reference rows parsed from {path}")

    required = {"organism", "accession", "role"}
    missing = required - set(header or [])
    if missing:
        sys.exit(f"error: {path} is missing required column(s): {sorted(missing)}")

    # A reference set is only usable if the sample fractions form a distribution.
    frac = sum(float(r.get("theoretical_dna_fraction") or 0)
               for r in rows if r["role"] == "sample")
    if abs(frac - 1.0) > 0.01:
        print(f"warning: sample theoretical_dna_fraction sums to {frac:.6f}, not 1.0",
              file=sys.stderr)
    return rows


def find_genome_files(genome_dir: Path, accession: str):
    """Locate the FASTA(s) datasets wrote for one accession."""
    for base in (genome_dir / "ncbi_dataset" / "data" / accession, genome_dir / accession):
        if base.is_dir():
            hits = sorted(list(base.glob("*.fna")) + list(base.glob("*.fna.gz"))
                          + list(base.glob("*.fa")) + list(base.glob("*.fasta")))
            if hits:
                return hits
    hits = sorted(genome_dir.rglob(f"{accession}*.fna")) or \
           sorted(genome_dir.rglob(f"*{accession}*.fna.gz"))
    if hits:
        return hits
    sys.exit(f"error: no FASTA found for accession {accession} under {genome_dir}")


def open_maybe_gz(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def main():
    args = parse_args()
    genome_dir = Path(args.genome_dir)
    rows = read_reference_tsv(args.reference_tsv)

    seen_contigs = {}
    provenance = []

    with open(args.out_fasta, "w") as out_fa, \
         open(args.out_contig_map, "w") as out_map, \
         open(args.out_genome_sizes, "w") as out_sizes:

        out_map.write("contig\torganism\trole\taccession\n")
        out_sizes.write("organism\trole\taccession\tgenome_size\tn_contigs\ttheoretical_dna_fraction\n")

        for row in rows:
            organism, accession, role = row["organism"], row["accession"], row["role"]
            files = find_genome_files(genome_dir, accession)

            total_len, n_contigs, cur = 0, 0, None
            digest = hashlib.sha256()

            for f in files:
                with open_maybe_gz(f) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            cur = line[1:].split()[0]
                            if cur in seen_contigs:
                                sys.exit(
                                    f"error: duplicate contig id {cur!r} in {organism} "
                                    f"({accession}); already claimed by {seen_contigs[cur]}. "
                                    "Competitive assignment requires unique contig ids.")
                            seen_contigs[cur] = organism
                            n_contigs += 1
                            out_map.write(f"{cur}\t{organism}\t{role}\t{accession}\n")
                            out_fa.write(line)
                        else:
                            seq = line.strip()
                            total_len += len(seq)
                            digest.update(seq.upper().encode())
                            out_fa.write(line)

            if n_contigs == 0:
                sys.exit(f"error: {accession} ({organism}) contributed no sequence")

            out_sizes.write(f"{organism}\t{role}\t{accession}\t{total_len}\t{n_contigs}\t"
                            f"{row.get('theoretical_dna_fraction', '')}\n")
            provenance.append({
                "organism": organism, "accession": accession, "role": role,
                "genome_size": total_len, "n_contigs": n_contigs,
                "sequence_sha256": digest.hexdigest(),
                "source_files": [f.name for f in files],
            })
            print(f"[refs] {organism:38s} {accession:18s} {role:12s} "
                  f"{total_len:>12,} bp  {n_contigs:>4d} contig(s)")

    if args.out_provenance:
        with open(args.out_provenance, "w") as fh:
            json.dump({"reference_set": args.reference_tsv,
                       "genomes": provenance}, fh, indent=2)

    print(f"[refs] combined reference: {len(seen_contigs):,} contigs, "
          f"{sum(g['genome_size'] for g in provenance):,} bp")


if __name__ == "__main__":
    main()
