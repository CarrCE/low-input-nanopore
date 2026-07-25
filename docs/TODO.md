# Known gaps

An honest list of what is not done, not sourced, or not yet trustworthy in this
repository. Items 1–7 are the ones that block publication or that a reader of the
results needs to know about; 8 onward are smaller loose ends.

## 1. Reads are not deposited, and `--fetch_from_sra` is not implemented

The FASTQs for `lowinput_s1` (r1–r3) and `lowinput_s2` (r0–r3) have not been
submitted to SRA/ENA. Consequences today:

- Every samplesheet row must carry a local `fastq` path. `parseSamplesheet` in
  `main.nf` exits with an error if that column is blank, and again if the file
  does not exist.
- The `sra_accession` column in `assets/samplesheets/*.csv` is reserved but
  empty in every row.
- `params.fetch_from_sra` exists in `nextflow.config` but is a placeholder:
  setting it exits immediately with
  `error: --fetch_from_sra is not implemented yet`.

To close: deposit the reads, populate `sra_accession`, then add a
`FETCH_READS` process (`fasterq-dump`/ENA FTP, pinned in a container) and switch
the samplesheet parser to prefer the accession when `--fetch_from_sra` is set.
The container for it does not exist yet either — `docker/tools` has no SRA
toolkit.

## 2. Only the competitive mode is wired

`nextflow.config` and `--help` advertise `--mode competitive | sequential | both`,
and `main.nf` validates all three, but the workflow body only builds the
competitive path. `--mode sequential` and `--mode both` currently pass validation
and then silently run the competitive path anyway, which is worse than failing.

Not implemented:

- the sequential subtraction chain (carrier → contaminant → community);
- the `--breseq_consensus` path, i.e. building a reference-guided consensus of
  the *E. coli* actually present in the carrier prep and subtracting reads that
  match it. The `withLabel: breseq` container
  (`quay.io/biocontainers/breseq:0.40.1--h3be2455_0`) is pinned in
  `conf/base.config` and its platform/root quirks are documented, but no process
  carries that label;
- the per-organism delta table that `--mode both` is supposed to emit — which is
  the artifact that would actually quantify how much the sequential method costs
  (see `ecoli-partitioning.md`).

Interim mitigation: reject `sequential` and `both` at parse time instead of
accepting them.

## 3. No aggregation across replicates, and no figures

`COMPUTE_METRICS` emits one `metrics.tsv` and one `summary.json` per sample. There
is no step that:

- pools replicates by the `experiment` column into mean/SD/CI per organism;
- honours `include_in_headline` when computing the headline statistics (the
  column is parsed into `meta` and then never used downstream);
- produces any manuscript figure — the analysis container ships matplotlib and
  nothing calls it;
- uses the per-read `readlengths.tsv.gz` for the adaptive-sampling
  ejection-signature analysis it was written for.

`modules/local/` and `comparison/` are empty placeholders for this work.

## 4. Mojarro 2019 reads/bases values are unsourced literals

The prior-study comparison carries reads-per-femtogram and bases-per-femtogram
figures attributed to Mojarro et al. 2019 that were transcribed as bare numbers.
They need a citation with the exact table/figure they came from, and a
restatement of that paper's denominator (DNA into library prep vs DNA into
extraction vs cell count) so it can be confirmed that it matches the convention
in `bin/compute_metrics.py`. Until sourced, do not report them.

## 5. The Zorzano comparison mixes two classifiers across the two axes

The comparison against Zorzano et al. currently puts values derived from one
classifier on one axis and values derived from a different classifier on the
other. Classifier choice changes both the numerator (what counts as assigned)
and the effective denominator, so the two axes are not measured on a common
basis and the comparison is not yet defensible. Resolve by reanalysing both
studies through a single, pinned classifier (see item 6) or by stating the
comparison as classifier-conditional with both values shown.

## 6. Prior-study reanalysis needs pinned databases

Any reanalysis of prior studies must pin exactly what it classified against, or
it is not reproducible and not comparable:

- Kraken2 PlusPF-8 database: `k2_pluspf_08gb_20241228.tar.gz`
- NCBI taxonomy: `new_taxdump_2025-01-01.zip`

Neither is currently referenced anywhere in the pipeline: there is no Kraken2
container, no database download step, and no record of these versions in
`conf/base.config`. Pin both (URL + checksum) alongside the container tags before
running the comparison.

## 7. `lowinput_s2_r0` has an unquantified input mass

For r0, fewer D6321 preps went into the extraction and the Qubit HS reading came
back off-scale low, so the sample DNA mass into library prep is not quantified.
Because every per-femtogram metric divides by that mass, r0:

- has a blank `library_dna_ng` in `assets/samplesheets/lowinput_s2.csv`;
- is flagged `include_in_headline=0`;
- yields blank `dna_fg`, `reads_per_fg`, `bases_per_fg` and a null
  `input_sample_fraction`/`enrichment` in its outputs (by design — the code emits
  blanks rather than a fabricated denominator).

It should be reported separately as a detection-at-unquantified-input result, not
folded into the headline per-fg statistics. That separation is currently a
convention in the samplesheet only; see item 3 — nothing enforces it downstream.

---

## 8. Declared parameters that nothing reads

`params.min_qscore`, `params.min_readlen`, `params.coverage_window`,
`params.community`, `params.carrier_accession` and `params.contaminant_accession`
are defined in `nextflow.config` but no process uses them. In particular
`COVERAGE_PROFILE` emits per-base `samtools depth` output and never bins by
`coverage_window`. Either implement them or delete them — advertised knobs that
do nothing are a reproducibility hazard.

## 9. No automated test beyond the smoke run

`-profile docker,test` runs 40,000 reads end to end, but nothing asserts anything
about the result: there is no expected-output check, and no CI workflow. At
minimum, assert that the read accounting reconciles and that the counts file
contains the expected organism classes.

## 10. `data/readme.md` and the reference TSVs describe the same genomes twice

The accession tables in `data/readme.md` and in `assets/references/*.tsv` are
maintained separately and can drift. The TSVs are authoritative (the pipeline
reads them); `data/readme.md` should point at them rather than restate them.

## 11. Genome-set caveats not yet assessed

Two `lowinput_s1` references are contig-level rather than complete
(*S. cerevisiae* `GCA_030867715.1`, *C. neoformans* `GCA_028975465.1`). Complete
chromosome-level substitutes from different strains are noted in `data/readme.md`
(`GCF_000146045.2`, `GCF_000149245.1`). The effect of the fragmented references on
competitive assignment — a fragmented reference can lose alignment-score
competitions at contig boundaries — has not been measured.
