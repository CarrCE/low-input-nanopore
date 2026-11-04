"""Display names for experiments and replicates in figures.

Figures say "D6311" and "D6321 r2", not the internal ids "lowinput_s1" and
"lowinput_s2_r2". The name is the community standard's catalogue number, read
from the standard_catalogue column of assets/measurements.tsv, so it lives with
the other facts about each replicate rather than in a mapping in plotting code.
A pack-size suffix is dropped: lowinput_s2 used D6321 and D6321-10, "the same
product in a ten-pack" (measurements.tsv header), and the product is what a
figure names. An experiment the file does not name keeps its id.

aggregate_results.py carries its own copy of read_experiment_names(), because
the AGGREGATE process stages that script alone; keep the two in step.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

DEFAULT = Path(__file__).resolve().parent.parent / "assets" / "measurements.tsv"


def read_experiment_names(path=DEFAULT) -> dict[str, str]:
    path = Path(path) if path else None
    if not path or not path.exists():
        return {}
    lines = [l for l in path.read_text().splitlines()
             if l.strip() and not l.lstrip().startswith("#")]
    names: dict[str, str] = {}
    for r in csv.DictReader(lines, delimiter="\t"):
        exp = (r.get("experiment") or "").strip()
        cat = (r.get("standard_catalogue") or "").strip().split("-")[0]
        if exp and cat:
            if names.get(exp, cat) != cat:
                sys.exit(f"error: {path}: experiment {exp} lists two standards "
                         f"({names[exp]} and {cat}); its display name is ambiguous.")
            names[exp] = cat
    return names


def sample_label(sample_id: str, names: dict[str, str]) -> str:
    """lowinput_s2_r3 -> "D6321 r3"; an unrecognised id is returned unchanged."""
    m = re.fullmatch(r"(.+)_(r\d+)", sample_id)
    if not m or m.group(1) not in names:
        return sample_id
    return f"{names[m.group(1)]} {m.group(2)}"
