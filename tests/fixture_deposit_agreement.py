#!/usr/bin/env python3
"""Assert that no committed fixture read leaks sequence the deposit withholds.

`assets/testdata/` ships real reads from this study, and this study's reads are
released to SRA in MASKED form. The fixture has to hold the read as it was
BEFORE masking -- that is what the masker takes as input -- so a fixture read
whose released form is masked publishes, inside the repository, sequence the
deposit deliberately withholds.

No test of the masker can catch that, and not because such a test would be
wrong: a test of the masker asks whether the right bases were masked, which is a
question about behaviour. This is a question about what the repository
DISTRIBUTES, and the two are independent. The masker can be perfectly correct on
a read that should never have been committed in the first place.

That is exactly what happened, and nothing misbehaved when it did. Six
`conserved_region` reads were picked by rules about ATTRIBUTION -- an `assigned`
call, >=50% organism coverage -- while what decides release is a different
question: whether >=`--chimera_min_bp` of the read is human-exclusive. All six
carried 160-2,211 such bases, so the deposit masks 25,900 of their 96,117, and
`tests/human_masking.py` masked them too, correctly, on every run. It passed
because it deliberately tolerates a minority of genuine chimeras there. The
fixture still shipped the unmasked input. The leak was found by hand, comparing
against the deposit; this test is what makes that comparison automatic, and it
is phrased in terms of the RELEASED OUTCOME so that it holds whichever path to
masking a read takes.

Two checks, because the deposit is 41 GB and is not present in a fresh clone:

  1. ALWAYS -- no fixture read may appear in `deposited_masked_ids.txt`, the ids
     the deposit masks over the pool the fixture draws from. Cheap, offline, and
     enough to catch a regeneration that reintroduces the same reads.
  2. WHEN THE READS ARE PRESENT -- every study-derived fixture read is compared
     base for base against its deposited counterpart. This is the real
     invariant; check 1 is its offline shadow. Point it at the reads with
     --deposit, or leave them in the `sra/` download cache.

`deposited_masked_ids.txt` is not a frozen intermediate of the kind that caused
this bug. It is derived from the submitted files themselves -- a read is masked
iff its released sequence contains an N -- and the deposit is immutable once
submitted. If it ever needs regenerating, that means a NEW submission, and the
fixture has to be rebuilt against it regardless.

Run:  python3 tests/fixture_deposit_agreement.py [--deposit reads.fastq.gz]
"""

import argparse
import glob
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "assets", "testdata")

# The pool the fixture draws from: the leading records of this replicate. The
# generator reads a head-truncated copy, so a read outside this window cannot be
# in the fixture and the committed id list need not cover the whole file.
SOURCE_SAMPLE = "lowinput_s1_r1"
POOL_RECORDS = 500_000


def read_fastq(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip("\n")
            fh.readline()
            q = fh.readline().rstrip("\n")
            yield h[1:].split()[0], s, q


def study_reads():
    """Fixture reads derived from this study, as {read_id: category}.

    Chimeras are spliced from a study read and a GIAB read; the manifest records
    only the study half's 8-character prefix, so they are matched by prefix and
    reported separately.
    """
    full, prefixes = {}, {}
    with open(os.path.join(DATA, "human_masking_manifest.tsv")) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rid, category, source = f[0], f[1], f[5]
            if source == SOURCE_SAMPLE:
                full[rid] = category
            elif category == "chimera":
                prefixes[source.split(":")[1].split()[0]] = rid
    return full, prefixes


def find_deposit(explicit):
    if explicit:
        return explicit
    for pattern in (f"sra/{SOURCE_SAMPLE}.fastq.gz",
                    f"sra/{SOURCE_SAMPLE}*.fastq.gz",
                    f"data/{SOURCE_SAMPLE}.masked.fastq.gz"):
        hits = sorted(glob.glob(os.path.join(ROOT, pattern)))
        if hits:
            return hits[0]
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--deposit", default=None,
                   help=f"the deposited {SOURCE_SAMPLE} reads; default: the "
                        f"sra/ download cache")
    p.add_argument("--require-deposit", action="store_true",
                   help="fail rather than skip when the reads are absent. Use "
                        "this before tagging a release.")
    args = p.parse_args()

    full, prefixes = study_reads()
    if not full:
        sys.exit(f"error: no fixture reads sourced from {SOURCE_SAMPLE}; the "
                 f"manifest's source column has changed and this test is no "
                 f"longer checking anything")

    fail = 0

    # ---- 1. the offline check -------------------------------------------
    id_path = os.path.join(DATA, "deposited_masked_ids.txt")
    masked = {l.strip() for l in open(id_path) if l.strip()}
    if not masked:
        sys.exit(f"error: {id_path} is empty")

    leaked = sorted(r for r in full if r in masked)
    masked_prefixes = {r[:8] for r in masked}
    leaked += sorted(cid for pref, cid in prefixes.items()
                     if pref in masked_prefixes)
    if leaked:
        print(f"FAIL {len(leaked)} fixture read(s) are masked in the deposit; "
              f"the repository would publish sequence SRA withholds:")
        for r in leaked[:10]:
            print(f"       {r}\t{full.get(r, 'chimera')}")
        fail += 1
    else:
        print(f"  ok   none of the {len(full)} study-derived fixture reads "
              f"({len(prefixes)} chimera halves included) is masked in the "
              f"deposit")

    # ---- 2. the real invariant, when the reads are here ------------------
    deposit = find_deposit(args.deposit)
    if deposit is None or not os.path.isfile(deposit):
        msg = (f"the deposited {SOURCE_SAMPLE} reads are not present, so the "
               f"base-for-base comparison was NOT run. Fetch them with "
               f"`--fetch_from_sra`, or pass --deposit.")
        if args.require_deposit:
            print(f"FAIL {msg}")
            fail += 1
        else:
            print(f"  SKIP {msg}")
    else:
        want = set(full)
        want_pref = dict(prefixes)
        fixture = {}
        for rid, s, _q in read_fastq(os.path.join(DATA, "human_masking.fastq.gz")):
            if rid in full or rid in prefixes.values():
                fixture[rid] = s

        differs, n = [], 0
        for rid, s, _q in read_fastq(deposit):
            n += 1
            if rid in want:
                want.discard(rid)
                if s != fixture[rid]:
                    differs.append((rid, full[rid], s.count("N")))
            pref = rid[:8]
            if pref in want_pref:
                cid = want_pref.pop(pref)
                # The chimera's study half is a prefix of this read. Whatever
                # the deposit masks inside that prefix must be masked in the
                # fixture too, so any N the deposit carries and the fixture does
                # not is a leak.
                fx = fixture[cid]
                span = min(len(fx), len(s))
                bad = sum(1 for i in range(span)
                          if s[i] == "N" and fx[i] != "N")
                if bad:
                    differs.append((cid, "chimera", bad))
            if not want and not want_pref:
                break
            if n >= POOL_RECORDS:
                break

        if want or want_pref:
            print(f"FAIL {len(want) + len(want_pref)} fixture read(s) were not "
                  f"found in the first {n:,} records of {deposit}. The fixture "
                  f"and the deposit do not describe the same reads.")
            for r in sorted(want)[:10]:
                print(f"       {r}")
            fail += 1
        elif differs:
            print(f"FAIL {len(differs)} fixture read(s) differ from their "
                  f"deposited form:")
            for rid, cat, extra in differs[:10]:
                print(f"       {rid}\t{cat}\t{extra} base(s) the deposit masks "
                      f"and the fixture does not")
            fail += 1
        else:
            print(f"  ok   all {len(full)} study-derived reads and "
                  f"{len(prefixes)} chimera halves are byte-identical to "
                  f"{os.path.basename(deposit)}")

    if fail:
        print(f"\nFAILED: {fail} check(s)")
        return 1
    print("\nfixture/deposit agreement: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
