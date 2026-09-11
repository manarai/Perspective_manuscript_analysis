"""Repressilator -- oscillatory (damped-focus regime).

Source
------
Elowitz & Leibler (2000) Nature 403:335.

Dynamics (cyclic, three genes, monotonic Hill repression)
    dx_i/dt = alpha / (1 + x_{i-1}^n) - x_i,   i = 1..3 (cyclic)

Operating point -- explicit, damped-focus, NOT the limit cycle
--------------------------------------------------------------
The symmetric fixed point x_1* = x_2* = x_3* = x* solves x* + x*^{n+1} =
alpha. Its Jacobian is a 3x3 circulant

    A = [[-1, 0, -beta], [-beta, -1, 0], [0, -beta, -1]]

with beta = alpha * n * x*^{n-1} / (1 + x*^n)^2. The eigenvalues are

    lambda_0 = -1 - beta
    lambda_{1,2} = -1 + beta/2 +/- i * beta * sqrt(3)/2

so the fixed point is a stable focus for beta < 2 and Hopf-unstable
for beta > 2. In the oscillating regime (beta >> 2) the fixed point is
UNSTABLE and linearising there gives Re lambda > 0, describing escape
from the fixed point rather than the limit cycle -- a Floquet analysis
on the orbit would be needed and produces a different object.

We choose parameters deliberately in the damped-focus regime -- the
fixed point is stable, oscillations are transient, and the operator
descriptors of A are well-defined at the fixed point. This is the
correct linearisation to compare to the other models on the same
footing. n = 2, alpha = 5.0 gives beta ~ 1.39 (safely below the Hopf).

Region
------
Oscillatory: leading eigenvalue is a complex-conjugate pair with
|Im|/|Re| > 1. Non-normal in the eigenvector geometry.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from ._base import Model


def _fixed_point(alpha: float, n: float) -> float:
    """Solve x + x^{n+1} = alpha for x > 0."""
    return brentq(lambda x: x + x ** (n + 1) - alpha, 1e-6, alpha + 1.0)


def build(alpha: float = 5.0, n: float = 2.0,
          decay_spread: float = 0.0) -> Model:
    """Repressilator Jacobian at the symmetric fixed point.

    decay_spread : float in [0, 1)
        0.0 gives the standard symmetric model (all decay rates = 1); the
        Jacobian is circulant and therefore exactly normal (nu = 0). Non-
        zero values break the symmetry with per-gene decay rates
        1 - s, 1, 1 + s. The fixed point is no longer symmetric-scaled
        and the Jacobian is no longer circulant, so nu > 0 while
        oscillation is retained. This variant is the "heterogeneous
        repressilator" cited in the manuscript.
    """
    x_star = _fixed_point(alpha, n)
    beta = alpha * n * x_star ** (n - 1) / (1.0 + x_star ** n) ** 2
    if beta >= 2.0:
        raise ValueError(
            f"beta = {beta:.3f} >= 2 => Hopf-unstable fixed point; "
            f"choose alpha / n to give beta < 2 for damped-focus regime"
        )
    d1, d2, d3 = 1.0 - decay_spread, 1.0, 1.0 + decay_spread
    A = np.array(
        [
            [-d1, 0.0, -beta],
            [-beta, -d2, 0.0],
            [0.0, -beta, -d3],
        ]
    )
    S = np.diag([x_star, x_star, x_star])
    A = np.linalg.inv(S) @ A @ S
    key = "repressilator" if decay_spread == 0.0 else "repressilator_het"
    name_suffix = "" if decay_spread == 0.0 else f" (heterogeneous, spread={decay_spread})"
    return Model(
        key=key,
        name=f"Repressilator (damped-focus){name_suffix}",
        domain="synthetic circuit",
        region="oscillatory",
        A=A,
        x_eq=np.array([x_star, x_star, x_star]),
        params={"alpha": alpha, "n": n, "beta": beta, "x_star": x_star,
                "decay_spread": decay_spread},
        citation="Elowitz & Leibler 2000 Nature 403:335",
        operating_point=(f"symmetric fixed point x*={x_star:.3f}; "
                         f"beta={beta:.3f} < 2 (stable focus, damped oscillation)"),
        scaling="states scaled by fixed-point magnitude",
        notes=("EXPLICITLY the damped-focus regime, not the limit cycle. "
               "Symmetric variant (decay_spread=0) is circulant and normal by "
               "construction; heterogeneous decay rates break the symmetry "
               "and produce non-normal oscillator geometry."),
    )
