# The 15 Aug 2026 history rewrite

This repository's git history was rewritten once, before it was ever public, to
remove read data that the SRA deposit withholds. This note records what was
removed and why, because a reader who notices that early commit SHAs do not
match anything published elsewhere deserves an explanation, and because a
silently rewritten history is worth less than a documented one.

## What was in there

`assets/testdata/human_masking.fastq.gz` is a fixture of real reads that
`bin/mask_human.py` is tested against. It has to hold reads as they were
**before** masking — that is what the masker takes as input — and 694 of its
reads come from this study's own `lowinput_s1_r1`.

Six of those six hundred and ninety-four should not have been there. Their
released form in the deposit is masked: 25,900 of their 96,117 bases are served
as `N` by SRA, of which 3,989 bp are intervals that align to T2T-CHM13v2.0 and
that no organism in the reference panel accounts for. The repository was
therefore distributing sequence the deposit deliberately does not.

## How it happened

Nothing malfunctioned, which is the part worth understanding.

The fixture generator selected `conserved_region` reads on two criteria: a
positive `assigned` call, and at least 50% of the read covered by its organism.
All six pass comfortably — 53.7% to 93.6% lambda coverage. But both criteria ask
about **attribution**, and attribution is not the only thing that decides
release. A confidently attributed read is still masked if at least
`--chimera_min_bp` (150) of it is human-exclusive. All six carry 160–2,211 such
bases.

`tests/human_masking.py` masked those six correctly on every run it ever made.
It passed because it deliberately tolerates a minority of genuine chimeras among
conserved-region reads, which is the right allowance — those reads really are
chimeras. The gap was never in the masker or in the test. A test of the masker
asks *were the right bases masked*; whether the repository may distribute a read
is a different question, and correctness on the first says nothing about the
second.

## What was done

1. The six reads were removed from the fixture, which is why `conserved_region`
   holds 194 reads rather than 200. See `assets/testdata/README.md`.
2. `tests/fixture_deposit_agreement.py` was added and wired into `make check`.
   It compares every study-derived fixture read against its **released** form —
   base for base when the deposited reads are present, and against a committed
   list of the ids the deposit masks when they are not. It is phrased in terms
   of the released outcome rather than any masking rule, so it holds whichever
   path to masking a read takes.
3. `tests/make_masking_testset.py` now requires `--deposited-masked-ids` and
   filters the candidate set once, so a regeneration cannot reintroduce them.
4. The history was rewritten with `git filter-repo`, replacing the fixture blobs
   at the single commit that introduced them. Every commit, message, author and
   date was preserved.

## What the SHAs do and do not tell you

Of the 107 commits, **75 kept their original hashes** — everything up to 29 July
2026. From `Merge pull request #27` onward they changed, and not only because of
the fixture: thirteen merge commits had been **GPG-signed by GitHub's web
interface**, and a signature cannot survive the content it signs being rewritten,
so `filter-repo` strips them. Removing a signature changes the commit object,
which changes its hash, which changes every descendant's parent pointer. The
cascade reaches commits whose own content never changed at all.

The consequence worth stating plainly: **those merges no longer carry a verified
signature**, and any hash cited elsewhere for a commit after 29 July 2026 will
not resolve here. This is unavoidable in any history rewrite; the alternative was
to leave the withheld reads in place.

Verified afterwards by scanning every blob of every commit in the rewritten
history for the six read ids and for 60-base probes taken from inside the
withheld regions: zero sequence matches. The six ids do still appear, in
`assets/testdata/deposited_masked_ids.txt`, which names them on purpose — an id
is not sequence, and anyone can derive that same list from the public deposit.

## What this does not claim

A rewrite reaches this repository, not copies of it. The pre-rewrite history
exists in a private repository retained as the development record, and in the
first author's local working copy. Both are held by people already trusted with
the unmasked reads, which is the only reason that is acceptable rather than a
gap.

Nothing about the published results changed. The rewrite touched only
`assets/testdata/`, which no part of the pipeline reads.
