# Partitioning *E. coli*: why the contaminant cannot simply be subtracted

## The problem

The lambda carrier is not pure lambda. Lambda is propagated in an *E. coli* K-12
host, and genomic DNA from that production host is carried over into the prep.
When ~1000 ng of lambda carrier goes into library prep alongside well under a
nanogram of sample DNA, even a small mass fraction of carried-over host DNA is
large compared with the sample. Any honest accounting therefore has to include
*E. coli* K-12 as a third category alongside "community" and "carrier" — which is
why both reference sets list it with `role=contaminant`
(`assets/references/lowinput_s1.tsv`, `assets/references/lowinput_s2.tsv`).

That alone would be routine. The complication specific to `lowinput_s1` is that
its community *also* contains *E. coli*:

| | Community member | Carrier-derived contaminant |
|---|---|---|
| Organism | *Escherichia coli* | *Escherichia coli* K-12 |
| Strain | B-1109 (Zymo D6311) | K-12 substr. MG1655 |
| Accession | `GCF_028743555.1` | `GCF_000005845.2` |
| Genome size | 4,925,141 bp | 4,641,652 bp |
| Source | Zymo Research submission (Zymo's own NCBI BioProject) | NCBI reference genome |
| Role in the reference set | `sample`, `theoretical_dna_fraction = 0.00089` | `contaminant`, fraction 0 |

These are genuinely different strains. They differ by roughly 283 kb of assembled
sequence and by strain-specific accessory content on both sides, while sharing a
highly similar core genome. That difference is the whole point: because the two
are distinguishable over part of their length, reads **can be partitioned**
between them rather than discarded wholesale.

## What the original analysis did, and why it was wrong

The original `lowinput_s1` analysis took the sequential route. It used `breseq`
to build a reference-guided consensus of the *E. coli* actually present in the
carrier prep — a sensible idea in itself, since the production strain is not
byte-identical to MG1655 — and then subtracted every read that matched that
consensus before mapping the survivors to the community.

The subtraction does not know which *E. coli* a read came from. Any read falling
in the shared core genome matches the consensus and is deleted, regardless of
whether it originated from the K-12 background or from the community's own
B-1109. The result:

- the community's *E. coli* is largely erased along with the contaminant;
- *E. coli*'s measured abundance collapses toward whatever fraction of B-1109 is
  strain-specific, with no way to tell that from a genuine dropout;
- because D6311 is a log-distributed standard and *E. coli* sits at
  `theoretical_dna_fraction = 0.00089`, the damage lands squarely in the part of
  the recovery curve the study is trying to characterise;
- the loss is invisible in the output. Subtracted reads leave no record, so
  nothing in the result set says how many reads were deleted or which organisms
  they were plausibly from.

The same failure mode applies to any carrier/community pair sharing sequence; the
*E. coli*/*E. coli* case here is just the sharpest instance.

## The solution: competitive assignment

`bin/assign_reads.py` never subtracts. Every read is mapped once
(`minimap2 -ax map-ont -N 10 --secondary=yes`) against a single combined index
holding the carrier, the contaminant, and every community member, so all
candidate organisms compete for the same read at the same time.

For each read:

1. Keep the best alignment score (`AS`) achieved against **each organism**.
   Supplementary alignments to the same organism are not summed — that would let
   a fragmented reference outscore a contiguous one.
2. Rank organisms by that score. Let `best` be the winner and `runner_up` the
   next organism.
3. Require a margin: `threshold = max(10, 0.01 * best_AS)` — 10 alignment-score
   units in absolute terms, or 1% of the winning score, whichever is larger. The
   fractional term keeps the test meaningful for long reads, where a 10-point
   difference is noise; the absolute floor keeps it meaningful for short ones.
4. If `best_AS - runner_up_AS >= threshold`, assign the read to `best`.
   Otherwise emit it as an explicit ambiguous class naming every tied organism:
   `ambiguous:Escherichia coli|Escherichia coli K-12`.
5. A read that only one organism aligned at all is assigned to that organism.
   Ambiguity is a statement about two organisms being indistinguishable, not
   about a single alignment being weak; calling weak-but-unique alignments
   "ambiguous" would both mislabel them and drain reads out of the per-organism
   counts.

So `lowinput_s1`'s *E. coli* reads land in three honest bins — B-1109-specific,
K-12-specific, and core-genome-ambiguous — and the contaminant contribution can be
*estimated from the strain-specific fraction* rather than assumed. Nothing is
deleted, the ambiguous class is reported in `counts.tsv` and `metrics.tsv` like
any other class, and `compute_metrics.py` fails the run if the classes do not sum
to the FASTQ read count.

Per-read evidence is preserved in `<sample>.assignments.tsv.gz`
(`read_id, organism, call, role, as_best, as_runnerup, margin, read_length,
aligned_bases`), so the margin distribution and the size of the ambiguous class
are auditable rather than a matter of trust. When reporting *E. coli* abundance
from `lowinput_s1`, report the ambiguous class alongside it — the strain-specific
counts are a lower bound on each strain, and the ambiguous count is the
irreducible overlap.

## An observed side effect: lambda/K-12 ambiguity

The ambiguous class is not only an *E. coli* phenomenon. Runs of these data
produce an `ambiguous:Enterobacteria phage lambda|Escherichia coli K-12` class as
well, which is expected: *E. coli* K-12 MG1655 carries lambda-related prophage
sequence in its own genome, so reads from those regions genuinely align well to
both references. Under sequential subtraction those reads are removed at the
first step and never accounted for. Here they are visible, counted, and kept out
of both organisms' unique tallies.

## `lowinput_s2` avoids the problem entirely

The D6321 spike-in community — *Truepera radiovictrix*, *Imtechella
halotolerans*, *Allobacillus halotolerans* — shares no homology with lambda or
with *E. coli*. There is no strain to confuse with the contaminant and no core
genome to tie on, so the partitioning problem cannot arise. That is what makes
`lowinput_s2` the clean dataset for the headline enrichment result, and
`lowinput_s1` the dataset that stresses the method.

*E. coli* K-12 is still included in `assets/references/lowinput_s2.tsv` as a
contaminant reference. It is not there to be subtracted: it is there so that
carrier-derived *E. coli* reads are attributed explicitly instead of silently
landing in the `unassigned` bin, where they would be indistinguishable from
sequencing failure and would quietly inflate the apparent noise floor.

## Caveats

- The margin rule is a heuristic, not a likelihood. `max(10, 1% of AS)` was
  chosen to be strict enough that core-genome reads do not get awarded to one
  strain by chance, and loose enough that a genuinely strain-specific read is not
  called ambiguous. Both thresholds are exposed as
  `--min-margin-abs` / `--min-margin-frac` in `bin/assign_reads.py`; the
  sensitivity of the *E. coli* split to them has not been characterised.
- Competitive assignment recovers a *partition*, not a deconvolution. It bounds
  each strain's contribution; it does not apportion the shared core between them.
- A quantitative comparison of what sequential subtraction costs on these data —
  the `--mode both` delta table — is not implemented yet. See item 2 in
  `TODO.md`.
