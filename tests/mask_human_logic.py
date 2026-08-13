#!/usr/bin/env python3
"""Unit tests for the human-masking decision logic in bin/mask_human.py.

These cover the rule itself -- interval algebra and the keep/mask decision --
with no I/O and no containers, so they run anywhere in well under a second and
are wired into `make check`.

The end-to-end behaviour (FASTQ in, masked FASTQ out) is covered separately by
the committed read-level fixture; see assets/testdata/.

The cases that actually discriminate are the last two of the decision block:
a conserved region shared with human must NOT cost an attributed read, and a
chimera must lose its human half while keeping its microbial one. A test set
of obviously-human reads proves almost nothing.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

from mask_human import (  # noqa: E402
    apply_mask, complement, decide, merge, subtract, total_len,
)

CHIMERA_MIN_BP = 150


class IntervalAlgebra(unittest.TestCase):
    def test_merge_coalesces_overlapping_and_adjacent(self):
        self.assertEqual(merge([(0, 10), (5, 20), (30, 40)]), [(0, 20), (30, 40)])
        self.assertEqual(merge([(0, 10), (10, 20)]), [(0, 20)])
        self.assertEqual(merge([]), [])

    def test_merge_handles_nested(self):
        self.assertEqual(merge([(0, 100), (20, 30)]), [(0, 100)])

    def test_subtract(self):
        self.assertEqual(subtract([(0, 100)], [(40, 60)]), [(0, 40), (60, 100)])
        self.assertEqual(subtract([(0, 100)], [(0, 100)]), [])
        self.assertEqual(subtract([(0, 100)], []), [(0, 100)])

    def test_subtract_multiple_holes(self):
        self.assertEqual(subtract([(0, 100)], [(10, 20), (50, 60)]),
                         [(0, 10), (20, 50), (60, 100)])

    def test_complement_and_total_len(self):
        self.assertEqual(complement([(20, 30)], 100), [(0, 20), (30, 100)])
        self.assertEqual(complement([], 50), [(0, 50)])
        self.assertEqual(total_len([(0, 10), (20, 25)]), 15)


class MaskDecision(unittest.TestCase):
    def decide(self, length, human, target, call):
        return decide(length, merge(human), merge(target), call, CHIMERA_MIN_BP)

    def test_wholly_human_is_entirely_masked(self):
        mask, _reason, chimeric = self.decide(1000, [(0, 1000)], [], "unmapped")
        self.assertEqual(mask, [(0, 1000)])
        # Not a chimera: nothing microbial to be chimeric WITH. Counting these
        # as chimeras would report every human read as one.
        self.assertFalse(chimeric)

    def test_unalignable_junk_is_masked(self):
        # No positive attribution, so no reason to release it.
        mask, _reason, _c = self.decide(500, [], [], "no_survivor")
        self.assertEqual(mask, [(0, 500)])

    def test_attributed_read_is_untouched(self):
        mask, _reason, _c = self.decide(1000, [], [(0, 1000)], "assigned")
        self.assertEqual(mask, [])

    def test_conserved_region_does_not_cost_an_attributed_read(self):
        # The false-positive case: HRRT flags an organism read because a
        # conserved region matches human. The organism accounts for those
        # bases, so human-exclusive is 0 and the read is released intact.
        mask, _reason, chimeric = self.decide(1000, [(300, 420)], [(0, 1000)], "assigned")
        self.assertEqual(mask, [])
        self.assertFalse(chimeric)

    def test_chimera_loses_its_human_half_and_keeps_the_rest(self):
        mask, _reason, chimeric = self.decide(1000, [(0, 400)], [(400, 1000)], "assigned")
        self.assertEqual(mask, [(0, 400)])
        self.assertTrue(chimeric)

    def test_chimera_rule_overrides_positive_attribution(self):
        # The privacy-critical case: whole-read best-hit assignment cannot see
        # this, because the read wins its comparison on the microbial half.
        mask, _reason, _c = self.decide(1000, [(0, 600)], [(600, 1000)], "assigned")
        self.assertTrue(total_len(mask) >= 600)

    def test_sub_threshold_human_exclusive_does_not_trip_the_rule(self):
        mask, _reason, chimeric = self.decide(1000, [(0, 100)], [(100, 1000)], "assigned")
        self.assertEqual(mask, [])
        self.assertFalse(chimeric)

    def test_ambiguous_keeps_only_what_an_organism_accounts_for(self):
        # Ambiguous is not positive attribution, so the read is masked -- but
        # organism-accounted bases still survive.
        mask, _reason, _c = self.decide(1000, [(0, 200)], [(200, 900)], "ambiguous")
        self.assertEqual(mask, [(0, 200), (900, 1000)])

    def test_rescue_requires_positive_attribution_not_absence_of_human(self):
        # A read with NO human alignment and NO organism alignment must still
        # be masked. Writing the rule as "mask unless it failed to look human"
        # would release this, which is backwards for a privacy filter.
        mask, _reason, _c = self.decide(800, [], [], "unmapped")
        self.assertEqual(mask, [(0, 800)])


class Masking(unittest.TestCase):
    def test_length_is_preserved(self):
        self.assertEqual(len(apply_mask("ACGT" * 250, [(0, 400)])), 1000)

    def test_blanks_exactly_the_named_interval(self):
        self.assertEqual(apply_mask("ACGTACGTAC", [(2, 5)]), "ACNNNCGTAC")

    def test_empty_mask_returns_the_sequence_unchanged(self):
        self.assertEqual(apply_mask("ACGTACGTAC", []), "ACGTACGTAC")

    def test_multiple_intervals(self):
        self.assertEqual(apply_mask("AAAAAAAAAA", [(0, 2), (8, 10)]), "NNAAAAAANN")


if __name__ == "__main__":
    unittest.main(verbosity=2)
