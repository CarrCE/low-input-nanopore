# Carrier-depletion model

How much depletion-mode adaptive sampling with a lambda carrier should enrich a
sample, and what the seven runs show instead. Backs Supplementary Section S6,
Figure S1 and Table S11 of the manuscript ([`paper-crosswalk.md`](paper-crosswalk.md)).

```bash
make carriermodel        # bin/carrier_model/run_carrier_model.sh
```

## Inputs

| Input | Where it comes from |
|---|---|
| Deposited, human-masked reads | SRA/ENA (`--fetch_from_sra` caches them in `sra/`), a `--mask_human` run's `results/<sample>/human/`, or a folder named by `READS=` |
| Read assignments and per-sample totals | a finished `make all` (`results/`) |
| Lambda genome | `results/references/<set>/combined.fasta`, contig `NC_001416.1` |
| MinKNOW instrument records | the instrument-data record, [doi:10.5281/zenodo.23907136](https://doi.org/10.5281/zenodo.23907136), fetched to `data/instrument/<sample>/` by `make instrumentdata` (or named by `INSTRUMENT=`); checksums in [`assets/instrument_files.tsv`](../assets/instrument_files.tsv) |

Every read-level quantity comes from the deposited reads, the ones a reader
downloads. Human masking replaced bases with N but kept read count and length,
so nothing here depends on unmasked data. The samplesheet's local `fastq`
column is not used for this reason: in the lab's own checkout it points at the
reads before masking.

Copies fetched from SRA/ENA work too. The archive renames each read by run and
spot (`SRR40180147.1`), keeps the ONT read id as the next header field, and
drops the ONT header tags. `extract_reads.py` recovers the ONT read id, which
keys the instrument records, and run time comes from the sequencing summary
rather than from the header tags.

The instrument records hold no sequence. For each run they are the Dorado
sequencing summary (sup basecalling), the adaptive-sampling decision and timing
logs, and pore activity, with rows restricted to the deposited reads.
`bin/carrier_model/build_instrument_record.py` made them from the MinKNOW output,
and `bin/carrier_model/verify_instrument_record.py` checks that the analysis
gives identical arrays from them as from the originals. Every deposited read is
in its run's summary, and the summary's read length equals the deposited length
for all 60.4 million reads.

## Lambda-derived reads

The model counts **lambda-derived reads**: reads that competitive assignment
gave to lambda, plus reads it called ambiguous with lambda among the tied
organisms. In D6311 part of lambda's left arm is shared with a prophage in the
community's *E. coli* B-1109, so reads from it tie and are called ambiguous.
Leaving them out removes reads from one end of lambda only, which biases where
reads appear to start and the bases per capture, and it does so in D6311 but
hardly in D6321 (12-13% of lambda-derived reads in D6311, about 1% in D6321).
With about 1,000 ng of lambda against well under 1 ng of sample, a read that
fits lambda as well as anything else is lambda.

The paper's **carrier** read class (lambda-assigned only) is unchanged and is
used everywhere outside the model: the read-class accounting, Table S9 and
Figure 4. `results/carrier_model/carrier_definitions.tsv` gives both side by side.

## The model

The simple model assumes (A1) capture per molecule, independent of length;
(A2) sample molecules read whole; (A3) every carrier capture truncated at
`b_eject` = 465 bp, the pooled median carrier read. Then base enrichment is
`L_carrier / b_eject`, about 100x for intact lambda.

The generalisation (`bin/carrier_model/model.py`) replaces A3 with the measured
bases per lambda-derived read, which include the reads that escaped ejection,
and restates A1 per adapted **end**, since a ligation read starts at an adapter
and lambda's ends are of two kinds:

    N_s / N_c = (m_s / m_c) * (L_eff / Ln_s)                      (1)
    L_eff     = G / (phi + beta)                                  (2)
    f_cos     = phi / (phi + beta)                                (3)

`N` = reads, `m` = mass into library prep, `Ln_s` = mean sample read length,
`G` = 48,502 bp, `phi` = fraction of cos ends that take an adapter, `beta` =
breaks per genome, `f_cos` = fraction of lambda-derived reads whose 5' end lies
within [-60, +20] bp of a cos site (`cos_ends.py`; cos starts sit 5-15 bp
outside the cos coordinate, where end-prep filled in the 12-nt overhang).
Intact lambda with both ends adaptable has `L_eff = G`. Every other read class
is taken as observed, and enrichment is computed as `bin/compute_metrics.py`
does. The decomposition into `phi` and `beta` is exactly identified, so it is
tested only by `0 <= phi <= 1`, by agreement between replicates sharing a lambda
lot, and by the bound `exp(-beta) >= F`, where `F` is the fraction of escaped
cos-start reads that reach the opposite cos site.

## Findings

D6321 r0 is excluded from every range below, as from the paper's headline
statistics (its sample mass is extrapolated).

1. **The simple model over-predicts D6311 by 1.4-2.1x**, because 465 bp is a
   median. Ejection works as designed (modal ejected length 455 bp in every
   run), but the 1.5-3.8% of lambda-derived reads that escape ejection carry
   36-60% of lambda-derived bases, so the mean per capture is 695-1,112 bp.
2. **Intact lambda with the measured bases per read reproduces D6311 within
   11% and D6321 within 24%**, with no fitted parameter.
3. **Lambda is not intact, and most cos ends are not adapted.** Only 36-47% of
   lambda-derived reads start at a cos site. Equations (1)-(3) give 0.50-0.73
   breaks per genome and 34-58% of cos ends adapted (D6311 replicates agree:
   `phi` 0.34-0.37, `beta` 0.54-0.67). The two effects roughly cancel, which is
   why "intact lambda" predicts well. Hypothesis, untested: the cohesive ends
   anneal and escape end-prep. It predicts that heating lambda (65 C, then
   rapid cooling) before end-prep raises `f_cos` and lowers enrichment.
4. **D6321 r0 is physically inconsistent** (`phi` = 1.76): its reads imply a
   sequenceable sample mass of 0.064 ng against 0.223 ng extrapolated. For
   r1-r3 the implied masses are 0.21, 0.38 and 0.34 ng against 0.26, 0.32 and
   0.40 ng measured.
5. **Length independence holds within D6311.** Across organisms (*E. coli*
   excluded, see [`ecoli-partitioning.md`](ecoli-partitioning.md)), the
   within-replicate slope of log bases per fg against log mean read length is
   -0.02 (model: 0), and of log reads per fg -1.02 (model: -1), over a
   threefold range of read length.
6. **Escapes are keep decisions on a poor first chunk.** 81-96% of escaped
   lambda-derived reads carry an adaptive-sampling `sequence` (keep)
   decision; across the seven runs, 6 of 1.49 million are failed unblocks. Kept reads alone hold
   30-58% of lambda-derived bases. Their first chunk had a live mean qscore of
   5.0-6.0 against 10.3-11.5 for ejected reads and 10.3-11.6 for kept sample
   reads, and moved at 207-284 against 363-409 bases per second, so the live
   mapper could not recognise lambda and, in depletion mode, kept the read.
7. **Ejection degrades over each run.** Between the first 6 h and the last
   12 h, escapes rise from 1.1-2.6% to 3.3-8.0% of lambda-derived reads, the
   share that adaptive sampling kept from 1.3-3.6% to 3.6-10.0%, and the median
   unblocked read from 442-469 to 480-555 bp. Community reads per
   lambda-derived read change by -28% to +2%, so enrichment is front-loaded.
8. **Two episodes where ejection largely stopped** (D6321 r2 at 25.0-26.5 h;
   D6321 r1 at 29.75-31.25 h): the unblocked share of lambda-derived reads fell
   from about 0.86 to 0.30 for about 1.5 h, with no change in active channels,
   decision depth or basecaller latency, and no mux scan or refuel. Enrichment
   rose inside them. The cause is unknown.
9. **Recapture of the same molecule end is ruled out.** Across the seven runs, in 8,854 consecutive
   same-channel pairs whose first read was an ejected read starting at an
   internal lambda position, no second read starts at the same end (within
   20 bp, gap under 5 s); 95% upper bound about 0.03%.
10. **Decisions.** Ejection is decided at a median 312-327 bases, and the strand
    leaves a median 197-202 bases later. Sample reads wrongly ejected: 3.1-4.0%
    of reads and 0.6-0.7% of bases in D6311; under 0.1% in D6321.
11. **Pores at the start of each run** (minutes 5-15): 393-453 per flow cell.

## Outputs

| File | Contents |
|---|---|
| `results/summary/carrier-model.{pdf,png,csv,json}` | Figure S1 |
| `results/summary/carrier-model-table.{csv,json}` | Table S11, unrounded |
| `results/carrier_model/decomposition.tsv` | per run: read classes, length structure, every prediction |
| `results/carrier_model/cos_ends.tsv`, `cos_5prime_positions.tsv` | where lambda-derived reads start and end |
| `results/carrier_model/time_resolved.tsv` | 15-min bins from the first deposited read |
| `results/carrier_model/decisions.tsv`, `escapes.tsv`, `decision_quality.tsv` | adaptive-sampling decisions by read class |
| `results/carrier_model/instrument_time.tsv`, `pores.tsv` | run-clock time course and pore counts |
| `results/carrier_model/recapture.tsv`, `recapture_ends.tsv`, `keep_decisions_5prime.tsv` | recapture tests |
| `results/carrier_model/per_organism_s1.tsv` | finding 5 |
| `results/carrier_model/carrier_definitions.tsv` | carrier class vs lambda-derived reads |

## Open questions

- What stopped ejection during the two episodes, and why community captures
  rose then.
- Whether more decision chunks for reads with an unmappable, low-quality start
  would recover the 30-58% of lambda-derived bases now spent on kept reads.
