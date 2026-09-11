"""PU.1 - GATA1 mutual repression toggle -- low-dimensional decision.

Source
------
Huang, Guo, May & Enver (2007) Dev. Biol. 305:695 (haematopoietic
lineage-priming model); Chickarmane, Enver & Peterson (2009) PLoS
Comput. Biol. 5:e1000268.

Dynamics
--------
    dx/dt = alpha * x^n / (K^n + x^n) + beta / (K^n + y^n) - gamma * x
    dy/dt = alpha * y^n / (K^n + y^n) + beta / (K^n + x^n) - gamma * y

Both factors self-activate and cross-repress. In the parameter range
used by Huang et al. the system is tristable: a symmetric unstable
fixed point on the x = y diagonal, and two stable asymmetric fixed
points where one factor dominates (PU.1-hi / GATA1-lo, and vice versa).

Operating point
---------------
The stable PU.1-hi / GATA1-lo fixed point. Linearisation gives a 2x2
Jacobian with two real negative eigenvalues (stable node) but strongly
asymmetric off-diagonals -- the natural signature of a decision-region
operator.

Region
------
Low-dimensional decision (labelled `toggle`). Should be more
contractive than the amplifying models but with mild non-normality
from the self-activation asymmetry.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import fsolve

from ._base import Model


def _hill(u: float, n: float, K: float) -> float:
    return u ** n / (K ** n + u ** n)


def _dhill(u: float, n: float, K: float) -> float:
    Kn = K ** n
    un = u ** n
    return n * u ** (n - 1) * Kn / (Kn + un) ** 2


def _repression_derivative(u: float, n: float, K: float) -> float:
    Kn = K ** n
    un = u ** n
    return -n * u ** (n - 1) * Kn / (Kn + un) ** 2


def _rhs(state, alpha, beta, gamma, n, K):
    x, y = state
    dx = alpha * _hill(x, n, K) + beta / (K ** n + y ** n) - gamma * x
    dy = alpha * _hill(y, n, K) + beta / (K ** n + x ** n) - gamma * y
    return np.array([dx, dy])


def build(alpha: float = 1.0, beta: float = 1.0, gamma: float = 1.0,
          n: float = 4.0, K: float = 0.5,
          asymmetry: float = 0.0) -> Model:
    """PU.1 / GATA1 toggle Jacobian at the PU.1-hi stable attractor.

    asymmetry : float in [0, 1)
        0.0 gives fully symmetric parameters (Huang et al. 2007 default);
        the Jacobian happens to be symmetric under our nondimensionalisation
        and so nu = 0 by construction. Non-zero values scale alpha, beta,
        gamma for the GATA1 side by (1 - asymmetry), breaking the exchange
        symmetry — probes whether nu = 0 at this attractor is a symmetry
        artifact or a robust feature.
    """
    alpha_y = alpha * (1.0 - asymmetry)
    beta_y = beta * (1.0 - asymmetry)
    gamma_y = gamma * (1.0 - asymmetry)
    def rhs(state, *args):
        x, y = state
        dx = alpha * _hill(x, n, K) + beta / (K ** n + y ** n) - gamma * x
        dy = alpha_y * _hill(y, n, K) + beta_y / (K ** n + x ** n) - gamma_y * y
        return np.array([dx, dy])
    guess = np.array([1.2, 0.1])
    x_star, y_star = fsolve(rhs, guess)
    # Analytic Jacobian at (x*, y*), with (potentially) asymmetric params.
    J11 = alpha * _dhill(x_star, n, K) - gamma
    J12 = beta * _repression_derivative(y_star, n, K)
    J21 = beta_y * _repression_derivative(x_star, n, K)
    J22 = alpha_y * _dhill(y_star, n, K) - gamma_y
    J = np.array([[J11, J12], [J21, J22]])
    eigs = np.linalg.eigvals(J)
    if np.any(eigs.real >= 0):
        raise ValueError(
            f"fsolve returned an unstable fixed point (eigs={eigs}); "
            f"perturb the initial guess or check parameters"
        )
    S = np.diag([x_star, y_star])
    A = np.linalg.inv(S) @ J @ S
    key = "toggle" if asymmetry == 0.0 else "toggle_het"
    name_suffix = "" if asymmetry == 0.0 else f" (heterogeneous, asymmetry={asymmetry})"
    return Model(
        key=key,
        name=f"PU.1 - GATA1 toggle{name_suffix}",
        domain="haematopoiesis",
        region="toggle",
        A=A,
        x_eq=np.array([x_star, y_star]),
        params={"alpha": alpha, "beta": beta, "gamma": gamma, "n": n, "K": K,
                "asymmetry": asymmetry},
        citation="Huang et al. 2007 Dev. Biol. 305:695",
        operating_point=(f"PU.1-hi / GATA1-lo stable fixed point "
                         f"(x*={x_star:.3f}, y*={y_star:.3f})"),
        scaling="each state scaled by fixed-point magnitude",
        notes=("Stable node. Symmetric variant has nu = 0 exactly (Jacobian "
               "coincidentally symmetric under this scaling). Heterogeneous "
               "variant probes whether contractive assignment survives "
               "parameter perturbation."),
    )
