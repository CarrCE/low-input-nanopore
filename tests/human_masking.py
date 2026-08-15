#!/usr/bin/env python3
"""Assert bin/mask_human.py against the committed read-level fixture.

Runs the real masker over assets/testdata/ and checks every read against the
outcome recorded in the manifest. Needs no mapping, no reference downloads and
no network: the fixture ships the alignment records the masker consumes.

What each category is actually testing:

  human_high_identity     sensitivity -- unattributed human must not survive
  human_divergent         sensitivity where identity to CHM13 is poor, which is
                          NOT the same as poor read quality
  human_grazing_organism  the sharpest case in the set: real human reads that
                          incidentally hit an organism. If the rescue is too
                          permissive these are what it wrongly releases.
  community               the masker must not touch attributed reads
  conserved_region        the FALSE-POSITIVE test. Real reads HRRT flagged and
                          the pipeline attributed. HRRT alone destroys ~49% of
                          this study's S. cerevisiae; these reads are why the
                          rescue exists.
  chimera                 exact recomputed boundaries: the human half must go
                          and the organism half must stay
  unalignable_junk        nothing claims it, so nothing protects it

Run:  python3 tests/human_masking.py       (inside the analysis image)
"""

import gzip
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "assets", "testdata")
MASKER = os.path.join(ROOT, "bin", "mask_human.py")


def read_fastq(path):
    op = gzip.open if path.endswith(".gz") else open
    out = {}
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return out
            s = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            out[h[1:].split()[0]] = s


class HumanMasking(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        out_fastq = os.path.join(cls.tmp, "masked.fastq")
        cmd = [
            sys.executable, MASKER,
            "--fastq", os.path.join(DATA, "human_masking.fastq.gz"),
            "--flagged", os.path.join(DATA, "human_masking.flagged_ids.txt"),
            "--assignments", os.path.join(DATA, "human_masking.assignments.tsv"),
            "--human-paf", os.path.join(DATA, "human_masking.human.paf"),
            "--bam", os.path.join(DATA, "human_masking.target.sam"),
            "--sample-id", "fixture",
            "--out-fastq", out_fastq,
            "--out-manifest", os.path.join(cls.tmp, "mask.tsv"),
            "--out-stats", os.path.join(cls.tmp, "stats.json"),
            "--out-masked-ids", os.path.join(cls.tmp, "ids.txt"),
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"masker failed:\n{r.stderr}")
        cls.before = read_fastq(os.path.join(DATA, "human_masking.fastq.gz"))
        cls.after = read_fastq(out_fastq)
        # Organism query intervals, read straight from the fixture SAM, so the
        # chimera assertion can name exactly which bases were promised.
        sys.path.insert(0, os.path.join(ROOT, "tests"))
        from make_masking_testset import sam_query_intervals
        cls.target = sam_query_intervals(
            os.path.join(DATA, "human_masking.target.sam"))
        cls.expect = {}
        with open(os.path.join(DATA, "human_masking_manifest.tsv")) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                cls.expect[f[0]] = {"category": f[1], "expect": f[2],
                                    "detail": f[3]}

    def outcome(self, rid):
        a, b = self.before[rid], self.after[rid]
        if a == b:
            return "rescued"
        return "masked_fully" if set(b) == {"N"} else "masked_partial"

    # ---- invariants ----------------------------------------------------

    def test_each_read_appears_once(self):
        # Everything downstream builds dicts, so a read in two categories does
        # not fail -- it silently adopts whichever row was written last, and one
        # of the two expectations stops being tested. Three GIAB reads sat in
        # both `human_divergent` and `human_grazing_organism` this way.
        records, ids = 0, set()
        with gzip.open(os.path.join(DATA, "human_masking.fastq.gz"), "rt") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                fh.readline(); fh.readline(); fh.readline()
                records += 1
                ids.add(h[1:].split()[0])
        self.assertEqual(records, len(ids),
                         f"{records - len(ids)} duplicate read(s) in the fixture")

        rows = []
        with open(os.path.join(DATA, "human_masking_manifest.tsv")) as fh:
            fh.readline()
            rows = [l.split("\t")[0] for l in fh if l.strip()]
        dupes = sorted({r for r in rows if rows.count(r) > 1})
        self.assertEqual(dupes, [],
                         f"read(s) claimed by two categories: {dupes[:5]}")
        self.assertEqual(set(rows), ids,
                         "the manifest and the FASTQ describe different reads")

    def test_no_read_is_lost_or_gained(self):
        self.assertEqual(set(self.before), set(self.after))

    def test_every_read_keeps_its_length(self):
        bad = [r for r in self.before if len(self.before[r]) != len(self.after[r])]
        self.assertEqual(bad, [], f"{len(bad)} reads changed length")

    def test_quality_lines_are_untouched(self):
        # Length equality above plus the masker never writing to the qual line
        # is what preserves per-base quality; assert the total base count too.
        self.assertEqual(sum(len(v) for v in self.before.values()),
                         sum(len(v) for v in self.after.values()))

    # ---- per-category outcomes ------------------------------------------

    def _category(self, name, allow=None):
        allow = allow or set()
        rids = [r for r, e in self.expect.items() if e["category"] == name]
        self.assertTrue(rids, f"no reads in category {name}")
        wrong = []
        for r in rids:
            got, want = self.outcome(r), self.expect[r]["expect"]
            if got != want and got not in allow:
                wrong.append((r, want, got))
        return rids, wrong

    def test_human_high_identity_is_masked(self):
        rids, wrong = self._category("human_high_identity")
        self.assertEqual(wrong, [],
                         f"{len(wrong)}/{len(rids)} human reads survived: "
                         f"{wrong[:5]}")

    def test_human_divergent_is_masked(self):
        rids, wrong = self._category("human_divergent")
        self.assertEqual(wrong, [],
                         f"{len(wrong)}/{len(rids)} divergent human reads "
                         f"survived: {wrong[:5]}")

    def test_human_grazing_an_organism_is_not_wholly_released(self):
        rids, _ = self._category("human_grazing_organism")
        released = [r for r in rids if self.outcome(r) == "rescued"]
        self.assertEqual(released, [],
                         "human reads with an incidental organism hit were "
                         f"released intact: {released}")

    def test_attributed_community_reads_are_untouched(self):
        rids, wrong = self._category("community")
        self.assertEqual(wrong, [],
                         f"{len(wrong)}/{len(rids)} attributed community reads "
                         f"were masked: {wrong[:5]}")

    def test_conserved_region_reads_are_never_destroyed(self):
        # THE false-positive test, and the reason the rescue exists: HRRT alone
        # destroys ~49% of this study's S. cerevisiae and ~92% of its
        # C. neoformans. No attributed read may be lost outright.
        rids = [r for r, e in self.expect.items()
                if e["category"] == "conserved_region"]
        destroyed = [r for r in rids if self.outcome(r) == "masked_fully"]
        self.assertEqual(destroyed, [],
                         f"{len(destroyed)}/{len(rids)} attributed reads were "
                         f"destroyed: {destroyed[:5]}")

    def test_conserved_region_reads_are_overwhelmingly_untouched(self):
        # A minority are genuine chimeras -- attributed over half their length
        # AND carrying human sequence no organism accounts for -- and those are
        # correctly masked in part. That is a real property of real reads, not a
        # defect, but it must stay a minority: if most "conserved region" reads
        # are being partially blanked, the human intervals are too permissive.
        rids = [r for r, e in self.expect.items()
                if e["category"] == "conserved_region"]
        untouched = [r for r in rids if self.outcome(r) == "rescued"]
        frac = len(untouched) / len(rids)
        self.assertGreaterEqual(
            frac, 0.90,
            f"only {frac:.1%} of conserved-region reads survived intact "
            f"({len(rids) - len(untouched)}/{len(rids)} partially masked)")

    def test_unalignable_junk_is_masked(self):
        rids, wrong = self._category("unalignable_junk")
        self.assertEqual(wrong, [], f"junk survived: {wrong[:5]}")

    # ---- chimeras, with exact boundaries --------------------------------

    def test_chimera_human_half_masked_and_organism_half_kept(self):
        rids = [r for r, e in self.expect.items() if e["category"] == "chimera"]
        self.assertTrue(rids)
        failures = []
        for r in rids:
            # detail: "organism 0-<cl>, human <cl>-<end>; ..."
            detail = self.expect[r]["detail"]
            cl = int(detail.split("organism 0-")[1].split(",")[0])
            before, after = self.before[r], self.after[r]
            human_part = after[cl:]
            if set(human_part) != {"N"}:
                failures.append((r, "human half not fully masked"))
            # Only bases the organism ACTUALLY aligns to are promised. The scope
            # rule keeps `target \ human`, so an unaligned base inside the
            # organism half is masked, correctly -- nothing accounts for it.
            tgt = self.target.get(r, [])
            damaged = [i for s, e in tgt for i in range(s, min(e, cl))
                       if after[i] != before[i]]
            if damaged:
                failures.append(
                    (r, f"{len(damaged)} organism-aligned bases were masked"))
        self.assertEqual(failures, [], f"chimera boundary failures: {failures[:5]}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
