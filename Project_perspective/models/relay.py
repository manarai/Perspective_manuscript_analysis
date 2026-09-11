"""Contractive symmetric signalling relay -- contractive.

Design intent
-------------
The `contractive` region requires a model whose norm is monotonically
non-increasing under the flow -- i.e., G_max = 1 at t = 0 and never
returns above it. That happens exactly when A is normal (or dissipative
in some Lyapunov metric) with negative-real spectrum. A generic linear
cascade (bidiagonal decay with subdiagonal coupling) is non-normal and
can transiently amplify, so is NOT the right choice for this slot.

We use a symmetric diffusion-plus-decay relay -- a discretised
reaction -- diffusion chain with detailed balance. Physical
justification: reversible enzyme-cascade signalling where each stage is
in fast equilibrium with its neighbours (e.g., cytoplasmic diffusion
between compartments, or nearest-neighbour phosphotransfer with equal
forward/reverse rates). Detailed balance guarantees symmetric A and
therefore contractive dynamics.

Dynamics (chain of N nodes with symmetric nearest-neighbour coupling
and uniform decay)
    dx_i/dt = D (x_{i-1} - 2 x_i + x_{i+1}) - gamma * x_i,   i = 1..N
    Neumann boundaries: x_0 = x_1, x_{N+1} = x_N.

Region
------
Contractive: symmetric A, purely real negative spectrum, ||exp(At)||
monotonically decreasing. Reference point in the descriptor space
where Henrici nu = 0 exactly, G_max = 1, numerical abscissa equals
spectral abscissa.
"""
from __future__ import annotations

import numpy as np

from ._base import Model


def build(N: int = 4, D: float = 0.5, gamma: float = 0.2,
          decay_spread: float = 0.0) -> Model:
    """Relay Jacobian (tridiagonal, Neumann BCs, minus per-node decay).

    decay_spread : float in [0, 1)
        0.0 gives uniform per-node decay = gamma; the Jacobian is
        symmetric and hence normal (nu = 0 exactly). Non-zero values
        set gamma_i = gamma * (1 + decay_spread * (i - (N-1)/2) / ((N-1)/2))
        so decays span [gamma*(1 - decay_spread), gamma*(1 + decay_spread)]
        linearly. Diffusion off-diagonals remain symmetric; only the
        diagonal is perturbed, which does not break the symmetry of A
        (a diagonal shift preserves symmetry). This variant tests
        whether contractivity survives biologically plausible per-node
        parameter heterogeneity.
    """
    L = np.zeros((N, N))
    for i in range(N):
        if i > 0:
            L[i, i - 1] = D
            L[i, i] -= D
        if i < N - 1:
            L[i, i + 1] = D
            L[i, i] -= D
    if decay_spread == 0.0:
        decays = np.full(N, gamma)
    else:
        mid = (N - 1) / 2.0
        decays = gamma * (1.0 + decay_spread * (np.arange(N) - mid) / mid)
    A = L - np.diag(decays)
    key = "relay" if decay_spread == 0.0 else "relay_het"
    name_suffix = "" if decay_spread == 0.0 else f" (heterogeneous, spread={decay_spread})"
    return Model(
        key=key,
        name=f"Diffusion-decay relay{name_suffix}",
        domain="signalling",
        region="contractive",
        A=A,
        x_eq=np.ones(N),
        params={"N": N, "D": D, "gamma": gamma, "decay_spread": decay_spread},
        citation="Textbook (Murray, Mathematical Biology Vol I, 2002)",
        operating_point="linear model; the operator IS the Jacobian",
        scaling="dimensionless states",
        notes=("Symmetric part unchanged by per-node decay heterogeneity "
               "(diagonal perturbation preserves matrix symmetry), so nu = 0 "
               "holds under this perturbation by construction. To break "
               "normality one would need to make the off-diagonals directional "
               "(k_forward != k_backward), which converts the relay into a "
               "cascade -- a different mechanism (see TCR)."),
    )
