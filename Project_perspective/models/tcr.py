"""TCR kinetic proofreading cascade — amplifying / non-normal.

Source
------
McKeithan (1995) PNAS 92:5042; canonical non-normal biological cascade.
Ganguli, Huh & Sompolinsky (2008) PNAS 105:18970 use the same operator
structure to argue for functional non-normality in biology.

State
-----
C_1, ..., C_N are occupations of the N kinetic-proofreading intermediates
(sequential phosphorylation states of the TCR-pMHC complex). Ligand
concentration L drives entry into C_1 and is treated as fixed input for
the linearisation.

Dynamics
--------
    dC_1/dt = k_on * L - (k_off + k_p) * C_1
    dC_i/dt = k_p * C_{i-1} - (k_off + k_p) * C_i     (2 <= i < N)
    dC_N/dt = k_p * C_{N-1} - k_off * C_N

Operating point
---------------
Linearisation around the L-driven equilibrium; only the intra-cascade
Jacobian in C is used (the L input drops out under linearisation). The
matrix does not depend on L, so nondimensionalisation is not needed --
the cascade is intrinsically dimensionless (rates in 1/time).

Region
------
The subdiagonal-only structure is a classic non-normal cascade: transient
signal amplification through the chain, exponential decay of the
spectrum. Same skeleton as a Jordan block, which is maximally non-normal
at a given spectrum.
"""
from __future__ import annotations

import numpy as np

from ._base import Model


def build(N: int = 5, k_off: float = 0.05, k_p: float = 0.5) -> Model:
    """Cascade Jacobian for the N-stage proofreading chain."""
    A = np.zeros((N, N))
    for i in range(N):
        A[i, i] = -k_off - k_p if i < N - 1 else -k_off
    for i in range(1, N):
        A[i, i - 1] = k_p
    return Model(
        key="tcr",
        name=f"TCR kinetic proofreading ({N}-stage)",
        domain="immune signalling",
        region="amplifying",
        A=A,
        x_eq=np.ones(N),                      # intrinsically dimensionless
        params={"N": N, "k_off": k_off, "k_p": k_p},
        citation="McKeithan 1995 PNAS 92:5042; Ganguli et al. 2008 PNAS 105:18970",
        operating_point="linearisation at the L-driven cascade equilibrium",
        scaling="rates in 1/s; states dimensionless (occupation)",
        notes=("Subdiagonal-only structure; classic Jordan-like non-normal "
               "cascade. Spectrum is single eigenvalue with algebraic "
               "multiplicity effectively N."),
    )
