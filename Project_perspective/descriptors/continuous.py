"""Continuous-time operator descriptors for the analytic Jacobian A.

Definitions follow Perspective_Route_Revision_Spec.md §1. All descriptors
take a real square Jacobian A obtained by analytic linearisation of a
mechanistic ODE at a stated operating point, expressed in nondimensional
coordinates that scale each state variable by its equilibrium magnitude
(or by the characteristic concentration used in the source paper).

Dimensionful descriptors (alpha, tau, G_max, t_star, omega) are metric-
relative and must be reported with their scaling. Dimensionless
descriptors (Q, g, nu, r_PR) are time-scaling invariant and are the
primary axes for Fig 2.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional

import numpy as np
from scipy.linalg import expm, svdvals


@dataclass
class Descriptors:
    alpha: float          # spectral abscissa, max Re lambda_i
    tau: float            # memory time, -1/alpha (inf if alpha >= 0)
    gap: float            # (Re l1 - Re l2) / |Re l1|
    Q: float              # |Im l1| / |Re l1|, oscillation quality of leader
    r_PR: float           # participation-ratio rank of singular spectrum
    omega: float          # numerical abscissa, lambda_max((A+A^T)/2)
    G_max: float          # max_{t in [0,T]} ||exp(At)||_2
    t_star: float         # argmax of ||exp(At)||_2
    nu: float             # Henrici departure-from-normality, in [0,1]
    horizon: float        # T actually used for G_max sweep

    def as_dict(self) -> dict:
        return asdict(self)


def _sort_eigs_by_real(A: np.ndarray) -> np.ndarray:
    """Eigenvalues of A sorted by descending real part (ties broken by |Im|)."""
    lam = np.linalg.eigvals(A)
    order = np.lexsort((-np.abs(lam.imag), -lam.real))
    return lam[order]


def spectral_abscissa(A: np.ndarray) -> float:
    return float(np.max(np.linalg.eigvals(A).real))


def numerical_abscissa(A: np.ndarray) -> float:
    """omega = lambda_max((A + A^T)/2). Governs initial transient growth rate."""
    Sym = 0.5 * (A + A.T)
    return float(np.max(np.linalg.eigvalsh(Sym)))


def spectral_gap(A: np.ndarray) -> float:
    """(Re l1 - Re l2) / |Re l1|. NaN if fewer than 2 eigenvalues or |Re l1|==0."""
    lam = _sort_eigs_by_real(A)
    if lam.size < 2:
        return float("nan")
    l1r = lam[0].real
    if abs(l1r) < 1e-15:
        return float("nan")
    return float((l1r - lam[1].real) / abs(l1r))


def oscillation_quality(A: np.ndarray) -> float:
    """Q = |Im l1| / |Re l1| for the leading eigenvalue. 0 if l1 real."""
    lam = _sort_eigs_by_real(A)
    l1 = lam[0]
    if abs(l1.real) < 1e-15:
        return float("inf") if abs(l1.imag) > 1e-15 else 0.0
    return float(abs(l1.imag) / abs(l1.real))


def participation_ratio_rank(A: np.ndarray) -> float:
    """r_PR = (sum sigma_i)^2 / sum sigma_i^2 on the singular spectrum of A.

    Effective rank in [1, n]; small when a few singular values dominate,
    close to n when the spectrum is flat.
    """
    s = svdvals(A)
    num = float(np.sum(s)) ** 2
    den = float(np.sum(s * s))
    return num / den if den > 0 else float("nan")


def henrici_index(A: np.ndarray) -> float:
    """nu = sqrt(||A||_F^2 - sum |lambda_i|^2) / ||A||_F, in [0, 1].

    0 for normal A; approaches 1 as A becomes maximally non-normal at
    fixed spectrum.
    """
    lam = np.linalg.eigvals(A)
    frob2 = float(np.sum(A * A))
    sum_lam2 = float(np.sum((lam.conj() * lam).real))
    num2 = max(0.0, frob2 - sum_lam2)
    denom = float(np.sqrt(frob2)) if frob2 > 0 else 1.0
    return float(np.sqrt(num2) / (denom + 1e-15))


def transient_gain(A: np.ndarray,
                   horizon: Optional[float] = None,
                   n_grid: int = 800,
                   horizon_multiplier: float = 20.0,
                   min_horizon: float = 5.0) -> tuple[float, float, float]:
    """G_max, t_star, horizon actually used.

    If horizon is None, set T = max(min_horizon, horizon_multiplier / |alpha|)
    when alpha < 0, else T = min_horizon (the operator is unstable and
    ||exp(At)|| diverges — G_max here reports the peak inside [0, T] only).
    Spec §1: horizon must be reported alongside G_max.
    """
    alpha = spectral_abscissa(A)
    if horizon is None:
        if alpha < -1e-12:
            horizon = max(min_horizon, horizon_multiplier / abs(alpha))
        else:
            horizon = min_horizon
    ts = np.linspace(0.0, horizon, n_grid)
    gains = np.empty_like(ts)
    for i, t in enumerate(ts):
        gains[i] = np.linalg.norm(expm(A * t), 2)
    idx = int(np.argmax(gains))
    return float(gains[idx]), float(ts[idx]), float(horizon)


def compute(A: np.ndarray,
            horizon: Optional[float] = None,
            n_grid: int = 800) -> Descriptors:
    """All descriptors for a real Jacobian A."""
    A = np.asarray(A, dtype=np.float64)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"A must be square 2-D; got shape {A.shape}")

    lam = _sort_eigs_by_real(A)
    alpha = float(lam[0].real)
    tau = float(-1.0 / alpha) if alpha < 0 else float("inf")
    G_max, t_star, hz = transient_gain(A, horizon=horizon, n_grid=n_grid)
    return Descriptors(
        alpha=alpha,
        tau=tau,
        gap=spectral_gap(A),
        Q=oscillation_quality(A),
        r_PR=participation_ratio_rank(A),
        omega=numerical_abscissa(A),
        G_max=G_max,
        t_star=t_star,
        nu=henrici_index(A),
        horizon=hz,
    )
