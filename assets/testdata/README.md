# Human-masking test set

1,520 reads that `tests/human_masking.py` runs `bin/mask_human.py` over,
asserting every read against `human_masking_manifest.tsv`. No network, no
mapping, no reference download: the fixture ships the alignment records the
masker consumes, so `make check` runs it in well under a second.

Regenerate with `tests/make_masking_testset.py` (see its docstring for the
inputs). The output is committed so the fixture's provenance is reproducible
rather than asserted.

## Composition

| category | n | expected | what it tests |
|---|---:|---|---|
| `human_high_identity` | 500 | fully masked | sensitivity: unattributed human must not survive |
| `human_divergent` | 197 | fully masked | sensitivity where identity to CHM13 is poor — which is **not** the same as poor read quality |
| `human_grazing_organism` | 9 | not released intact | real human reads that incidentally hit an organism; if the rescue is too permissive, these are what it wrongly keeps |
| `community` | 500 | untouched | the masker must not touch attributed reads |
| `conserved_region` | 194 | not destroyed; ≥90% untouched | **the false-positive test** |
| `chimera` | 20 | human half masked, organism-aligned bases kept | exact, recomputed boundaries |
| `unalignable_junk` | 100 | fully masked | nothing claims it, so nothing protects it |

## The fixture may not hold what the deposit withholds

The study reads here come from the file **before** masking, because that is what
the masker has to be fed. So a fixture read whose released form is masked would
publish, inside this repository, sequence that SRA deliberately does not — and
no test of the masker can notice, because the masker is asserted against the
fixture's own frozen records. It agrees with itself while disagreeing with what
was actually released.

That is not hypothetical, and the way it happened is worth understanding,
because nothing misbehaved.

Six `conserved_region` reads were selected here by the two rules below — an
`assigned` call and ≥50% organism coverage — and emitted with `expect=rescued`.
Both rules are about **attribution**. Neither asks the question that actually
decides release: does the read carry at least `--chimera_min_bp` (150) of
sequence that no organism accounts for? The chimera rule is an independent path
to masking, and the selection logic never modelled it. All six carried 160–2,211
human-exclusive bp, so the real run masked 25,900 of their 96,117 bases — and so
did `tests/human_masking.py`, correctly, every time it ran. It did not fail
because the ≥90% allowance below exists precisely to tolerate a minority of
genuine chimeras among these reads. That allowance is right. What was missing is
that a read the masker masks is a read whose unmasked form must not be committed.

So this was never a behaviour bug or a drift: the fixture's expectation for those
six was wrong from the day it was written, and the leak is structural. The
fixture has to hold the **pre-masking** read in order to test masking, and
nothing asked whether that read was safe to publish. The generator had a floor
to keep its *expectations* honest and nothing at all to keep its *contents*
releasable. The six were removed on 15 Aug 2026, which is why
`conserved_region` is 194 rather than 200.

Two things now prevent a recurrence. `tests/make_masking_testset.py` requires
`--deposited-masked-ids` and excludes every listed read from every pool, so a
regeneration cannot reintroduce them. `tests/fixture_deposit_agreement.py`,
wired into `make check`, asserts the result: base for base against the deposited
reads when they are present, and against `deposited_masked_ids.txt` when they
are not.

Three GIAB reads were collapsed at the same time. `human_grazing_organism` was
filtered out of the high-identity pool but not the divergent one, so those three
carried two manifest rows with contradicting expectations — and because every
consumer builds a dict keyed on read id, the conflict resolved silently to
whichever row came last. `human_divergent` is 197 rather than 200 for that
reason, and `test_each_read_appears_once` now asserts it cannot recur. This one
was a stale test rather than a leak: GIAB reads are public and consented.

`deposited_masked_ids.txt` holds the 1,147 reads the deposit masks within the
leading 500,000 records of `lowinput_s1_r1` — the window the generator draws
from. It is derived from the submitted file itself (a read is masked iff its
released sequence contains an `N`), not from the masking run's own bookkeeping,
and the deposit is immutable once submitted. Regenerating it would mean a new
submission, at which point the fixture has to be rebuilt anyway.

## Provenance

Every read is real sequence except the junk.

- **Human** — Giab HG002/NA24385 (Ashkenazim son), UCSC ultralong ONT
  PromethION, `GM24385_1.fastq.gz`, from the leading 120 MB of
  `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/`.
  Public and consented for exactly this purpose. Reads were kept at
  300–12,000 bp to bound the fixture's size, then split at 0.90 identity to
  T2T-CHM13v2.0 (`GCF_009914755.1`), measured with `minimap2 -c` and the
  `de:f` tag.
- **Community and conserved-region** — this study's own `lowinput_s1_r1`
  reads, from the first 500,000 of that file, run through the real pipeline
  with `--mask_human`. The conserved-region reads are the actual false
  positives a real HRRT run produced, not an imagining of what one looks like.
- **Chimeras** — a real community read spliced to a real GIAB read at a
  recorded offset, so the expected boundary is exact rather than inferred.
- **Junk** — the only synthetic component: random sequence, seeded, so nothing
  should ever claim it.

## Two selection rules that matter

**Attributed reads need a coverage floor.** `bin/assign_reads.py` will attribute
a read to an organism on a very small footprint — a 8,118 bp read on a 55 bp
yeast hit — and such reads are frequently *human* reads with an incidental
organism match. Community and conserved-region reads therefore require ≥50% of
the read to be covered by its organism. Without that floor the fixture would
encode the wrong expected outcome for exactly the reads that discriminate.

Note what this floor does **not** cover: it reasons about attribution, and
attribution is not the only thing that decides an outcome. A well-covered,
confidently attributed read is still masked if ≥`--chimera_min_bp` of it is
human-exclusive. That gap is what put six reads in here with the wrong
expectation, and it is why `--deposited-masked-ids` is a required argument
rather than another rule of this kind — it is stated in terms of the released
outcome, so it holds no matter which path to masking a read takes.

**`conserved_region` is not asserted to be untouched.** A minority of these are
genuine chimeras — attributed over half their length *and* carrying human
sequence no organism accounts for — and are correctly masked in part. The
assertions are that none is ever destroyed, and that ≥90% survive intact. If
most were being partially blanked, the human intervals would be too permissive.

## A note on the alignment records

`human_masking.target.sam` carries **real query intervals** extracted from real
alignments, with `SEQ`/`QUAL` set to `*` and reference coordinates not
preserved. The masker consumes query intervals only, pysam recovers them from
the CIGAR alone (verified on every record), and storing the sequence twice would
double the fixture for nothing. Keeping it as text rather than BAM means the
expected intervals can be read directly. Chimera records are recomputed from the
splice offset, which is what makes their expected outcome exact.

Identity in `human_masking.human.paf` is recorded in the `de:f` tag. Do not
compute identity from PAF columns 10 and 11: they tie when minimap2 runs without
`-c`, and on this study's own data they scored alignments at 0.15–0.57 that
`-c` scores at 0.67–0.98.
