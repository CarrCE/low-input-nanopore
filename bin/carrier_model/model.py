#!/usr/bin/env python3
"""Carrier + depletion-mode adaptive sampling: the simple model, generalised.

The simple model (Supplementary Section S10) assumes

  (A1) capture is per molecule, independent of length;
  (A2) a sample molecule is read whole;
  (A3) a carrier molecule is truncated at b_eject = 465 bp,

giving enrichment_bases = L_carrier / b_eject. 465 bp is the pooled median of
the carrier reads, so the model treats every carrier capture as an ejection.

Generalisation.

A3 is replaced by measurement: the bases per carrier capture are observed,
including reads that escaped ejection. Every other read class (contaminant,
ambiguous not involving lambda, unassigned) is taken as observed too.

A1 is restated per adapted END, because a ligation read begins at an adapter
and the carrier's ends are of two kinds -- cos ends and breaks -- that need not
be equally adaptable. With the same per-end adaptation and capture efficiency
for every non-cos end:

  N_s / N_c = (m_s / m_c) * (L_eff / Ln_s)                                 (1)

  L_eff = 2 / (adaptable carrier ends per bp) = G / (phi + beta)          (2)

where N_c counts lambda-derived reads (see below), Ln_s is the sample's
number-average molecule length (estimated by its mean read length, A2), phi the
fraction of cos ends that can take an adapter, and beta the breaks per genome.
Intact lambda with both cos ends adaptable has L_eff = G = 48,502 bp, which
recovers the simple model with the measured mean in place of 465 bp. Where
lambda-derived reads begin (cos_ends.py) gives

  f_cos = phi / (phi + beta)                                               (3)

so (1)-(3) decompose an observed replicate into phi and beta. Physical
constraints: 0 <= phi <= 1, beta >= 0. The decomposition is exactly identified,
so it is a test only through those constraints, through agreement between
replicates sharing a lambda lot, and through the full-length bound
exp(-beta) >= (fraction of escaped cos-start reads that reach the other cos).

Lambda-derived reads. Every model quantity counts reads that competitive
assignment gave to lambda PLUS reads it called ambiguous with lambda among the
tied organisms. In D6311 a stretch of lambda's left arm is shared with a
prophage in E. coli B-1109, so reads from it tie and are called ambiguous;
leaving them out would remove reads from one end of lambda only, biasing f_cos
and the bases per capture. The paper's "carrier" read class (lambda-assigned
only) is unchanged and is used everywhere outside the model.
"""
from __future__ import annotations

AVOGADRO = 6.02214076e23
DA_PER_BP = 650.0
LAMBDA_BP = 48502
B_EJECT_SIMPLE = 465.0


def simple_enrichment(m_s, L_s, m_c, L_c=LAMBDA_BP, b_ej=B_EJECT_SIMPLE):
    """(bases, reads) enrichment under the simple model (A1-A3)."""
    n = lambda m, L: (m * 1e-9) / (L * DA_PER_BP) * AVOGADRO
    ns, nc = n(m_s, L_s), n(m_c, L_c)
    bc = min(L_c, b_ej)
    base_frac = ns * L_s / (ns * L_s + nc * bc)
    read_frac = ns / (ns + nc)
    mass_frac = m_s / (m_s + m_c)
    return base_frac / mass_frac, read_frac / mass_frac


def predict(obs: dict, L_eff: float, Ln_s: float | None = None) -> dict:
    """Predict sample reads, bases and enrichment from observed carrier captures.

    obs: m_s, m_c (ng); N_s, B_s (sample reads, bases); N_c (lambda-derived
    reads); N_tot, B_tot (all reads, bases). Ln_s defaults to mean sample
    read length.
    """
    b_s = obs["B_s"] / obs["N_s"]
    Ln_s = Ln_s or b_s
    N_s = obs["N_c"] * (obs["m_s"] / obs["m_c"]) * (L_eff / Ln_s)             # (1)
    B_s = N_s * b_s
    N_tot = obs["N_tot"] - obs["N_s"] + N_s
    B_tot = obs["B_tot"] - obs["B_s"] + B_s
    in_frac = obs["m_s"] / (obs["m_s"] + obs["m_c"])
    return {"N_s": N_s, "B_s": B_s,
            "E_reads": (N_s / N_tot) / in_frac, "E_bases": (B_s / B_tot) / in_frac}


def implied_L_eff(obs: dict, Ln_s: float | None = None) -> float:
    """Invert (1): the carrier mass per pair of adaptable ends the reads imply."""
    Ln_s = Ln_s or obs["B_s"] / obs["N_s"]
    return (obs["N_s"] / obs["N_c"]) * (obs["m_c"] / obs["m_s"]) * Ln_s


def decompose(L_eff: float, f_cos: float) -> tuple[float, float]:
    """(phi, beta) from (2) and (3)."""
    k = LAMBDA_BP / L_eff
    return f_cos * k, (1 - f_cos) * k


def implied_m_s(obs: dict, L_eff: float, Ln_s: float | None = None) -> float:
    """Invert (1) for the sample mass: the sequenceable sample mass the reads imply."""
    Ln_s = Ln_s or obs["B_s"] / obs["N_s"]
    return (obs["N_s"] / obs["N_c"]) * obs["m_c"] * (Ln_s / L_eff)


if __name__ == "__main__":
    for L_c in (48502, 23250, 20000, 18600):
        print(f"L_c={L_c:>6}  simple model = {simple_enrichment(0.1, 2000, 400.0, L_c)[0]:.0f}x")
