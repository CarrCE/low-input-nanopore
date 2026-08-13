#!/usr/bin/env python3
"""Build the committed human-masking test set from real reads.

Run once to regenerate assets/testdata/. The output is committed; this script
exists so that the fixture's provenance is reproducible rather than asserted.

Every read in the fixture is real sequence:

  * human            -- GIAB HG002/NA24385, public and consented
  * community        -- this study's own lowinput_s1 reads
  * conserved region -- this study's reads that HRRT flags but the pipeline
                        positively attributes to an organism. These are the
                        false positives the rescue exists to prevent, taken
                        from a real HRRT run rather than imagined.
  * chimera          -- a real community read spliced to a real human read at a
                        recorded offset, so the expected boundary is exact
  * junk             -- the only synthetic component; random sequence that
                        nothing should ever claim

The fixture ships its own alignment records so the test needs no mapping and no
reference downloads. Those records carry REAL query intervals extracted from
real alignments; reference coordinates are not preserved, because the masker
consumes query intervals only. Chimera records are recomputed from the splice
offset, which is what makes their expected outcome exact rather than inferred.
"""

import argparse
import gzip
import os
import random
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
from mask_human import merge, total_len  # noqa: E402

RNG = random.Random(20260812)


def read_fastq(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip("\n")
            p = fh.readline()
            q = fh.readline().rstrip("\n")
            yield h[1:].split()[0], s, q


def paf_intervals(path, ident_min=None):
    """read_id -> merged query intervals; optionally identity-filtered."""
    out, ident = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            m = re.search(r"de:f:([0-9.]+)", line)
            i = 1 - float(m.group(1)) if m else None
            if ident_min is not None and i is not None and i < ident_min:
                continue
            out.setdefault(f[0], []).append((int(f[2]), int(f[3])))
            if i is not None:
                ident[f[0]] = max(ident.get(f[0], 0), i)
    return {k: merge(v) for k, v in out.items()}, ident


def sam_query_intervals(path):
    """read_id -> merged query intervals, parsed from CIGAR."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t")
            if int(f[1]) & 4:
                continue
            cigar = f[5]
            if cigar == "*":
                continue
            qpos, start = 0, None
            for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                n = int(n)
                if op in "SH":
                    if start is None:
                        qpos += n
                    continue
                if op in "MI=X":
                    if start is None:
                        start = qpos
                    qpos += n
            if start is not None and qpos > start:
                out.setdefault(f[0], []).append((start, qpos))
    return {k: merge(v) for k, v in out.items()}


def bam_query_intervals(path, wanted):
    import pysam
    out = {}
    with pysam.AlignmentFile(path, "rb", check_sq=False) as fh:
        for a in fh.fetch(until_eof=True):
            if a.is_unmapped or a.query_name not in wanted:
                continue
            s, e = a.query_alignment_start, a.query_alignment_end
            if e > s:
                out.setdefault(a.query_name, []).append((s, e))
    return {k: merge(v) for k, v in out.items()}


def cigar_for(intervals, length):
    """Minimal CIGAR encoding a set of query intervals."""
    parts, cur = [], 0
    for s, e in intervals:
        if s > cur:
            parts.append(f"{s - cur}S")
        parts.append(f"{e - s}M")
        cur = e
    if cur < length:
        parts.append(f"{length - cur}S")
    return "".join(parts) or f"{length}S"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--giab-fastq", required=True)
    p.add_argument("--giab-human-paf", required=True)
    p.add_argument("--giab-community-sam", required=True)
    p.add_argument("--study-fastq", required=True)
    p.add_argument("--study-assignments", required=True)
    p.add_argument("--study-flagged", required=True)
    p.add_argument("--study-bam", required=True)
    p.add_argument("--study-human-paf", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--n-human-high", type=int, default=500)
    p.add_argument("--n-human-divergent", type=int, default=200)
    p.add_argument("--n-community", type=int, default=500)
    p.add_argument("--n-conserved", type=int, default=200)
    p.add_argument("--n-chimera", type=int, default=20)
    p.add_argument("--n-junk", type=int, default=100)
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)

    giab_h, giab_ident = paf_intervals(a.giab_human_paf)
    giab_t = sam_query_intervals(a.giab_community_sam)
    giab_seq = {rid: (s, q) for rid, s, q in read_fastq(a.giab_fastq)}

    flagged = {l.split()[0] for l in open(a.study_flagged) if l.strip()}
    calls = {}
    with gzip.open(a.study_assignments, "rt") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            calls[f[0]] = (f[2], f[1])           # call, organism
    study_h, _ = paf_intervals(a.study_human_paf)

    # Study reads worth including: positively attributed AND actually covered by
    # that organism. The coverage floor matters -- assign_reads will attribute a
    # read on a very small footprint, and such reads are frequently human reads
    # with an incidental organism hit. Including them would bake a wrong
    # expectation into the fixture.
    cand = {r for r, (c, _o) in calls.items() if c == "assigned"}
    study_t = bam_query_intervals(a.study_bam, cand)
    study_seq = {}
    for rid, s, q in read_fastq(a.study_fastq):
        if rid in cand:
            study_seq[rid] = (s, q)

    def covered(rid):
        L = len(study_seq[rid][0])
        return total_len(study_t.get(rid, [])) / L if L else 0.0

    community = [r for r in cand
                 if r not in flagged and r in study_seq and covered(r) >= 0.5]
    conserved = [r for r in cand
                 if r in flagged and r in study_seq and covered(r) >= 0.5]
    RNG.shuffle(community)
    RNG.shuffle(conserved)

    hi = [r for r, i in giab_ident.items()
          if i >= 0.90 and r in giab_seq
          and total_len(giab_h[r]) / len(giab_seq[r][0]) >= 0.5]
    lo = [r for r, i in giab_ident.items() if i < 0.90 and r in giab_seq]
    # The 9 human reads that graze an organism are the sharpest cases in the
    # set: if the rescue is too permissive, these are what it wrongly releases.
    grazers = [r for r in giab_t if r in giab_seq]
    hi = [r for r in hi if r not in grazers]
    RNG.shuffle(hi)
    RNG.shuffle(lo)

    fq = open(os.path.join(a.outdir, "human_masking.fastq"), "w")
    man = open(os.path.join(a.outdir, "human_masking_manifest.tsv"), "w")
    asg = open(os.path.join(a.outdir, "human_masking.assignments.tsv"), "w")
    paf = open(os.path.join(a.outdir, "human_masking.human.paf"), "w")
    sam = open(os.path.join(a.outdir, "human_masking.target.sam"), "w")
    ids = open(os.path.join(a.outdir, "human_masking.flagged_ids.txt"), "w")

    man.write("read_id\tcategory\texpect\texpect_detail\tread_len\tsource\n")
    asg.write("read_id\torganism\tcall\trole\tas_best\tas_runnerup\tmargin\t"
              "read_length\taligned_bases\n")
    sam.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:fixture\tLN:250000000\n")

    def emit(rid, seq, qual, cat, expect, detail, human, target, call, organism,
             source):
        fq.write(f"@{rid}\n{seq}\n+\n{qual}\n")
        ids.write(rid + "\n")
        man.write(f"{rid}\t{cat}\t{expect}\t{detail}\t{len(seq)}\t{source}\n")
        asg.write(f"{rid}\t{organism}\t{call}\tnone\t0\t0\t0\t{len(seq)}\t0\n")
        for s, e in human:
            n = e - s
            paf.write(f"{rid}\t{len(seq)}\t{s}\t{e}\t+\tchm13\t3117292070\t"
                      f"1\t{n}\t{n}\t{n}\t60\tde:f:0.0500\n")
        if target:
            # SEQ and QUAL are '*'. The masker reads only query intervals, which
            # pysam recovers from the CIGAR alone -- verified on all records --
            # and storing the sequence again would double the fixture's size for
            # nothing. Keeping the file as text rather than BAM means a reviewer
            # can read the expected intervals directly.
            sam.write(f"{rid}\t0\tfixture\t1\t60\t{cigar_for(target, len(seq))}"
                      f"\t*\t0\t0\t*\t*\n")

    n = 0
    for rid in hi[:a.n_human_high]:
        s, q = giab_seq[rid]
        emit(rid, s, q, "human_high_identity", "masked_fully",
             "no organism accounts for any base", giab_h.get(rid, []), [],
             "unmapped", "unassigned", "GIAB HG002 GM24385_1")
        n += 1
    for rid in lo[:a.n_human_divergent]:
        s, q = giab_seq[rid]
        emit(rid, s, q, "human_divergent", "masked_fully",
             "divergent from CHM13 but still unattributed",
             giab_h.get(rid, []), [], "unmapped", "unassigned",
             "GIAB HG002 GM24385_1")
        n += 1
    for rid in grazers:
        s, q = giab_seq[rid]
        emit(rid, s, q, "human_grazing_organism", "masked_partial",
             "human read with an incidental organism hit; only the organism-"
             "accounted bases may survive",
             giab_h.get(rid, []), giab_t.get(rid, []), "assigned",
             "incidental", "GIAB HG002 GM24385_1")
        n += 1
    for rid in community[:a.n_community]:
        s, q = study_seq[rid]
        emit(rid, s, q, "community", "rescued", "positively attributed",
             study_h.get(rid, []), study_t.get(rid, []), "assigned",
             calls[rid][1], "lowinput_s1_r1")
        n += 1
    for rid in conserved[:a.n_conserved]:
        s, q = study_seq[rid]
        emit(rid, s, q, "conserved_region", "rescued",
             "HRRT flagged it, an organism accounts for it: the false-positive "
             "case", study_h.get(rid, []), study_t.get(rid, []), "assigned",
             calls[rid][1], "lowinput_s1_r1")
        n += 1

    # ---- chimeras: exact, recomputed boundaries --------------------------
    pool_c = [r for r in community[a.n_community:] if covered(r) >= 0.8]
    pool_h = [r for r in hi[a.n_human_high:]
              if total_len(giab_h[r]) >= 400 and len(giab_seq[r][0]) >= 400]
    for i in range(min(a.n_chimera, len(pool_c), len(pool_h))):
        cr, hr = pool_c[i], pool_h[i]
        cs, cq = study_seq[cr]
        hs, hq = giab_seq[hr]
        cl = min(len(cs), 2000)
        hl = min(len(hs), 1500)
        seq, qual = cs[:cl] + hs[:hl], cq[:cl] + hq[:hl]
        tgt = merge([(s, min(e, cl)) for s, e in study_t.get(cr, []) if s < cl])
        hum = merge([(cl + s, cl + min(e, hl))
                     for s, e in giab_h.get(hr, []) if s < hl])
        rid = f"chimera_{i:03d}_{cr[:8]}_{hr[:8]}"
        emit(rid, seq, qual, "chimera", "masked_partial",
             f"organism 0-{cl}, human {cl}-{cl + hl}; the human half must be "
             f"masked and the organism half kept",
             hum, tgt, "assigned", calls[cr][1],
             f"spliced lowinput_s1_r1:{cr[:8]} + GIAB:{hr[:8]}")
        n += 1

    for i in range(a.n_junk):
        L = RNG.randint(300, 2000)
        seq = "".join(RNG.choice("ACGT") for _ in range(L))
        rid = f"junk_{i:03d}"
        emit(rid, seq, "5" * L, "unalignable_junk", "masked_fully",
             "nothing claims it, so nothing protects it", [], [],
             "no_survivor", "unassigned", "synthetic random sequence")
        n += 1

    for fh in (fq, man, asg, paf, sam, ids):
        fh.close()
    print(f"[testset] {n} reads written to {a.outdir}")


if __name__ == "__main__":
    main()
