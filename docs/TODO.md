# Known gaps

What is not done, not sourced, or not yet trustworthy in this repository.

Items 1–5 are open. Section 6 records questions that were open, were
investigated, and were closed; they are kept because the manuscript's
Supplementary Information refers to several of them, and because "we checked and
it did not matter" is a result a reader may want to audit rather than take on
trust.

---

## 1. Reads are not deposited, and `--fetch_from_sra` is not implemented

The FASTQs for `lowinput_s1` (r1–r3) and `lowinput_s2` (r0–r3) have not been
submitted to SRA/ENA. Until they are, nothing here is clone-and-run for anyone
outside the lab. Consequences today:

- Every samplesheet row must carry a local `fastq` path. `parseSamplesheet` in
  `main.nf` errors if that column is blank, and again if the file does not
  exist.
- The `sra_accession` column in `assets/samplesheets/*.csv` is reserved but
  empty in every row.
- `params.fetch_from_sra` exists in `nextflow.config` but is a placeholder:
  setting it exits immediately with
  `--fetch_from_sra is not implemented yet`.

To close: deposit the reads, populate `sra_accession`, then add a `FETCH_READS`
process (`fasterq-dump`/ENA FTP, pinned in a container) and switch the
samplesheet parser to prefer the accession when `--fetch_from_sra` is set. The
container does not exist yet either — `docker/tools` has no SRA toolkit.

## 2. The prior-study reanalysis is pinned but has not been re-run through the scripted path

`comparison/kraken2_db.manifest.tsv` records exactly what the `kraken2_q1` and
`kraken2_q10` rows classified against — URL, byte count, checksum, and the
provenance of each checksum:

- Kraken2 PlusPF-8, 2024-12-28 build (`k2_pluspf_08gb_20241228.tar.gz`,
  5,925,280,339 B, MD5 `01b8b1eb…`, published upstream)
- NCBI new_taxdump 2025-01-01 (139,761,991 B, MD5 `171470a1…`, SHA-256
  `7ff98c65…`, computed locally — NCBI publishes no sidecar for archived dumps)

Both URLs were verified to be the ones `epi2me-labs/wf-metagenomics` v2.14.1
(commit `a57ff73c…`) hard-codes for `--database_set PlusPF-8`, so the pin records
what was actually used rather than a plausible substitute.

**The gap:** the committed numbers were produced through the EPI2ME desktop
application and have *not* been regenerated through
`comparison/run_kraken2_reanalysis.sh`. Doing so would convert "pinned and
reproducible in principle" into "reproduced", and is the natural response if a
reviewer questions the prior-study values. It needs a ~5.5 GiB download, ~8 GiB
of RAM, and the Zorzano raw reads, which are not staged locally (the Basapathi
Raghavendra reads are, under `data/raghavendra_2023/`).

Nothing here touches the main pipeline, which does not use Kraken2 at all.

## 3. `COVERAGE_PROFILE` measures alignments, not assignments

`samtools depth` runs over the primary alignments to community contigs. It is
not restricted to the reads `assign_reads.py` awarded to each organism, so an
organism whose sequence is shared with an abundant relative accumulates depth
from that relative's reads. Pooled over the replicates:

| organism | alignment depth | attributable | ratio |
|---|---:|---:|---:|
| *E. faecalis* | 2.09× | 0.0016× | 0.1% |
| *E. coli* B-1109 | 158× | 2.08× | 1.3% |
| *S. enterica* | 0.71× | 0.18× | 25% |
| the other ten | — | — | 97–104% |

*E. faecalis* is awarded 4 reads in `lowinput_s1_r1` and shows 0.88× there:
1,942 reads tie between it and *Listeria*, and roughly half place their primary
alignment on its genome. For B-1109 the top 1% of 1 kb bins hold 79% of the
depth, peaking above 10⁵× near 3.06 and 3.09 Mb — λ carrier reads on
λ-related prophage sequence.

**Mostly addressed without re-running.** `bin/coverage_attribution.py`
(`make attribution`) computes both depths per replicate from `counts.tsv` and the
reference genome sizes, and the 1× interpretability threshold is now applied to
the attributable one. That is the correct criterion regardless of how depth is
measured, and applying it settles the question for the manuscript:

- 17 pairs clear 1× on alignment depth; **14** clear it on attributable depth
- the three that leave are exactly the *E. coli* B-1109 rows (0.83×, 0.62×,
  0.62× attributable), which should never have been reported as characterised
- nothing enters, and the 14 that remain agree between the two depths to within
  0.4%, so their statistics may be read as the organism's own

Both coverage figures mark the panels where the depths disagree rather than
dropping them, since the *E. coli* panel is the cautionary example the SI
explains.

The same script fixed a second error: the coverage globs match `test_s2`, the
40,000-read smoke-test subsample, so the denominator had been 45 pairs rather
than 42. `--exclude` now drops it in both the table and the figure.

**What is left.** Only *E. coli* B-1109 in the pooled figure, where attributable
depth is 2.08× and so above threshold, but the profile drawn is still alignment
depth. Showing its true pooled profile needs the BAM filtered by assignment
before `samtools depth`, which needs a full re-map — the qname BAMs are not
retained (only `test_s2.qname.bam`) and `assignments.tsv.gz` carries no alignment
coordinates. At 2.08× pooled, roughly 0.7× per replicate, that is marginal
territory where uniformity statistics are weak anyway. Deferred, not blocking.

## 4. Declared parameters that nothing reads

`params.min_readlen`, `params.community`, `params.carrier_accession` and
`params.contaminant_accession` are defined in `nextflow.config` but no process
uses them. Either implement them or delete them — advertised knobs that do
nothing are a reproducibility hazard.

(`params.min_qscore` and `params.coverage_window` were on this list and are now
both consumed.)

## 5. Genome-set caveats not assessed

Two `lowinput_s1` references are contig-level rather than complete
(*S. cerevisiae* `GCA_030867715.1`, *C. neoformans* `GCA_028975465.1`). Complete
chromosome-level substitutes from different strains are noted in
`data/readme.md` (`GCF_000146045.2`, `GCF_000149245.1`). A fragmented reference
can lose alignment-score competitions at contig boundaries; the size of that
effect has not been measured. Stated as a limitation in the manuscript.

Relatedly, `data/readme.md` and `assets/references/*.tsv` describe the same
genomes twice and can drift. The TSVs are authoritative — the pipeline reads
them; `data/readme.md` should point at them rather than restate them.

---

## 6. Closed

### `--breseq_consensus`: implemented, run, and found unnecessary

The original `lowinput_s1` analysis subtracted against a `breseq`
reference-guided consensus of the *E. coli* actually present in the carrier prep,
on the reasonable grounds that the λ production host need not be identical to
MG1655. That option is implemented here (`EXTRACT_CONTAMINANT_READS` →
`BRESEQ_CONSENSUS` → `MAP_CONSENSUS` → `assign_reads.py --consensus-hits`) so the
earlier work can be reproduced faithfully rather than approximated. Its
accounting is asserted by `make check` (`tests/consensus_accounting.py`, 7
checks).

It turns out not to be needed. `bin/contaminant_divergence.sh` runs breseq over
the competitively-assigned contaminant reads of a finished run:

| replicate | reads | mapped | depth | SNPs | indels | structural |
|---|---:|---:|---:|---:|---:|---:|
| s1_r1 | 159,229 | 92.3% | 31.2× | 20 | 0 | 0 |
| s1_r2 | 91,920 | 91.2% | 19.7× | 16 | 0 | 0 |
| s1_r3 | 76,826 | 91.1% | 19.1× | 10 | 1 | 0 |
| s2_r0 | 577,647 | 86.7% | 40.8× | 100 | 8 | 0 |
| s2_r1 | 435,677 | 86.5% | 31.5× | 84 | 3 | 0 |
| s2_r2 | 579,853 | 87.8% | 44.8× | 100 | 3 | 0 |
| s2_r3 | 446,495 | 89.4% | 39.9× | 69 | 3 | 0 |

86–92% of contaminant reads map to stock MG1655 in every replicate, no replicate
shows a single structural variant, and the worst case is ~100 SNPs in 4,641,652
bp — one per 46 kb. A 1 kb ONT read has a ~2% chance of overlapping even one.
Refining the reference to the strain actually present cannot move a read across
the assignment margin, so the stock reference is adequate. This also disposes of
an objection to the competitive-vs-sequential comparison: at this divergence,
subtracting against a consensus would have removed the same reads, so the
measured cost of subtraction is not an artifact of reference choice.

Depth is the operative constraint on the option: below roughly 10× breseq
predicts missing coverage across the whole reference and returns a deleted
genome instead of a consensus. `--breseq_min_depth` (default 10) catches that up
front, and a whole-reference `DEL` check catches it after. The bundled test
profile reaches ~0.3× and therefore cannot exercise the option.

**One observation left unresolved.** S2 calls roughly 4× the variants of S1, and
it is not a depth artifact: s1_r1 (31.2×) and s2_r1 (31.5×) are depth-matched and
call 20 vs 84. Calls are stable within each session — 9 positions recur in all
three S1 replicates, 52 in all four S2 — so the difference is real and between
sessions. Two candidate explanations this design cannot separate: a different NEB
λ lot (the sessions are 14 months apart), or seed sets that are not comparable in
kind, because S1's community contains *E. coli* B-1109 so ambiguous reads never
reach breseq, while in S2 nothing competes with *E. coli* and every
enterobacterial read is assigned to the contaminant. The second also explains
S2's lower mapping rate. Separating them needs a λ-only control library, which is
a different experiment. Reported in the SI as an open observation rather than
adjudicated.

### Coverage-artifact analysis

Was missing; `COVERAGE_PROFILE` wrote depth files nothing consumed. Now
`bin/coverage_summary.py` and `bin/plot_coverage.py` produce the coverage figure
and table, and `bin/coverage_dropouts.py` locates low-coverage regions and
annotates them against a GFF3. The two *Listeria* dropouts both turned out to be
mobile genetic elements at ~half median depth.

### Mojarro 2019 reads/bases were unsourced literals

Traced to Table 1 ("Low-Input Carrier Sequencing Metrics"), row "*B. subtilis*
reads": 5 reads, 5,270 bases, with the 2 pg input from that paper's abstract. The
row now carries `verified=TRUE` and `classifier=published_table1` in
`comparison/prior_studies.tsv`. The drop-unverified-rows flag is no longer
needed.

### The Zorzano comparison mixed two classifiers across the two axes

Resolved by carrying three variants of each condition. The published figure uses
`kraken2_q1`, where both axes come from the same reanalysis of the published raw
reads. The defective hybrid is retained only so the earlier figure can be
reproduced for audit, and is not published.

### `lowinput_s2_r0`'s input mass

Recorded as 0.223 ng in `assets/measurements.tsv`, with
`sample_dna_basis=raw_fluorescence_extrapolated` — an extrapolation from the raw
fluorescence of the calibration standards and the sample, below the Qubit HS
reporting range. It is excluded from the headline statistics on that basis, with
the reason recorded in `include_reason`, and `bin/check_measurements.py` enforces
that an exclusion carries a stated reason.

The earlier form of this entry described r0 as having *no* mass and being
excluded because the field was blank. That is the anti-pattern the measurements
file exists to prevent: an exclusion resting on absence reverses itself the
moment the value is filled in, which is exactly what happened once — filling in
0.223 silently returned r0 to a figure and moved a published mean.

### No automated test beyond the smoke run

`make check` (`tests/consensus_accounting.py`) asserts seven properties of the
consensus-subtraction accounting over the smoke-test BAM, including that no read
is lost. There is still no CI workflow, and no assertion over the full-run
outputs.
