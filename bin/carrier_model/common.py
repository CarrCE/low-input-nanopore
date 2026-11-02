"""Shared definitions for the carrier-model analysis (bin/carrier_model/)."""
from pathlib import Path

import numpy as np

# The seven replicates of this study, in run order.
SAMPLES = ["lowinput_s1_r1", "lowinput_s1_r2", "lowinput_s1_r3",
           "lowinput_s2_r0", "lowinput_s2_r1", "lowinput_s2_r2", "lowinput_s2_r3"]

# Display names: experiments by the standard they used, as in the manuscript.
EXPERIMENT_NAME = {"lowinput_s1": "D6311", "lowinput_s2": "D6321"}

LAMBDA = "Enterobacteria phage lambda"
ROLES = ["sample", "carrier", "contaminant", "ambiguous", "none"]

ESCAPE_BP = 2000           # lambda-derived reads at or above this escaped ejection
BIN_S = 900                # time bins, seconds


def display(sample_id):
    """lowinput_s1_r2 -> 'D6311 r2'."""
    exp, rep = sample_id.rsplit("_", 1)
    return f"{EXPERIMENT_NAME[exp]} {rep}"


def is_lambda_derived(role, organism):
    """A read the model counts as carrier: assigned to lambda, or an ambiguous
    call whose tied organisms include lambda. See model.py for why."""
    return role == "carrier" or (role == "ambiguous" and LAMBDA in organism)


def classes(d):
    """Per-read masks from a <sample>.reads.npz (extract_reads.py)."""
    role = d["role"]
    roles = list(d["roles"])
    orgs = list(d["orgs"])
    amb_lambda = np.array([o.startswith("ambiguous") and LAMBDA in o for o in orgs],
                          dtype=bool)[d["org"]]
    R = {r: role == roles.index(r) for r in roles}
    R["lambda_derived"] = R["carrier"] | (R["ambiguous"] & amb_lambda)
    R["ambiguous_lambda"] = R["ambiguous"] & amb_lambda
    return R


def write_tsv(path, rows):
    """Write a list of dicts as TSV, floats to 6 significant figures."""
    path = Path(path)
    cols = list(rows[0])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], (float, np.floating)) else str(r[c])
                               for c in cols) + "\n")


def read_tsv(path, comment="#"):
    import csv
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith(comment)), delimiter="\t"))
