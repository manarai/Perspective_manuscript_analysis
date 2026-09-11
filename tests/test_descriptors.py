"""Sanity tests for the descriptor module against hand-solvable A."""

import math
import numpy as np
import pytest

from descriptors.continuous import compute, henrici_index, numerical_abscissa
from nulls.spectrum_matched import spectrum_matched_normal


def test_normal_diagonal_has_zero_nonnormality():
    A = np.diag([-1.0, -2.0, -3.0])
    d = compute(A)
    assert d.nu == pytest.approx(0.0, abs=1e-12)
    assert d.alpha == pytest.approx(-1.0)
    assert d.tau == pytest.approx(1.0)
    assert d.omega == pytest.approx(-1.0)         # sym part = A
    assert d.G_max == pytest.approx(1.0, rel=1e-6)
    assert d.t_star == pytest.approx(0.0, abs=1e-6)


def test_pure_rotation_is_normal_and_oscillatory():
    A = np.array([[0.0, -1.0], [1.0, 0.0]])       # skew-symmetric, spectrum +/- i
    d = compute(A)
    assert d.nu == pytest.approx(0.0, abs=1e-12)
    assert d.alpha == pytest.approx(0.0, abs=1e-12)
    assert d.omega == pytest.approx(0.0, abs=1e-12)
    assert math.isinf(d.Q) or d.Q > 1e6           # Re -> 0
    assert d.G_max == pytest.approx(1.0, rel=1e-6)


def test_shear_block_is_maximally_nonnormal():
    """A = [[-1, k], [0, -1]] shares eigenvalues with -I but has k-dependent
    non-normality and transient amplification. G_max > 1 for large k."""
    k = 10.0
    A = np.array([[-1.0, k], [0.0, -1.0]])
    d = compute(A, horizon=10.0, n_grid=2000)
    assert d.alpha == pytest.approx(-1.0)
    assert d.nu > 0.5
    assert d.G_max > 2.0
    assert d.t_star > 0.0
    assert d.omega > d.alpha                       # numerical > spectral
    assert numerical_abscissa(A) == pytest.approx((k / 2) - 1.0, rel=1e-6)


def test_surrogate_spectrum_and_normality_2d_oscillatory():
    """Non-normal A with a complex-conjugate pair — the case the naive
    real-Schur surrogate would silently corrupt."""
    A = np.array([[-0.2, 3.0], [-1.0, -0.2]])      # eigenvalues -0.2 +/- i*sqrt(3)
    lam_src = np.sort_complex(np.linalg.eigvals(A))
    A_null = spectrum_matched_normal(A)
    lam_null = np.sort_complex(np.linalg.eigvals(A_null))
    np.testing.assert_allclose(lam_null, lam_src, atol=1e-10)
    err = np.max(np.abs(A_null @ A_null.T - A_null.T @ A_null))
    assert err < 1e-10
    # Henrici uses sqrt of a subtraction; near-normal matrices leak
    # machine-eps into the numerator, so tolerate ~sqrt(eps).
    assert henrici_index(A_null) == pytest.approx(0.0, abs=1e-6)
    assert henrici_index(A) > 0.3


def test_surrogate_higher_dim_mixed_real_and_complex():
    rng = np.random.default_rng(0)
    A = rng.standard_normal((6, 6))
    A_null = spectrum_matched_normal(A)
    lam_src = np.sort_complex(np.linalg.eigvals(A))
    lam_null = np.sort_complex(np.linalg.eigvals(A_null))
    np.testing.assert_allclose(lam_null, lam_src, atol=1e-8)
    assert henrici_index(A_null) < 1e-6


def test_participation_ratio_extremes():
    from descriptors.continuous import participation_ratio_rank
    # Flat singular spectrum: identity has r_PR = n
    assert participation_ratio_rank(np.eye(5)) == pytest.approx(5.0)
    # Rank-1 has r_PR = 1
    v = np.array([[1.0, 2.0, 3.0]])
    A = v.T @ v
    assert participation_ratio_rank(A) == pytest.approx(1.0, rel=1e-10)
