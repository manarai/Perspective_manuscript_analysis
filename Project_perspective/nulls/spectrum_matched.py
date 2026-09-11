"""Spectrum-matched normal surrogate.

Given a real Jacobian A, construct a real normal matrix A_null with
identical eigenvalues but orthogonal eigenvectors. Purpose: isolate the
contribution of eigenvector geometry (non-normality) from the spectrum
in Fig 2b of the Perspective revision.

Construction — real block-normal form
--------------------------------------
The naive "complex Schur, drop off-diagonals" recipe from the spec is
subtly wrong for a general non-normal real A: the columns of the
complex Schur Q do not come in conjugate pairs, so `Q diag(T) Q^H` has
a residual imaginary part. We use a construction that is unconditionally
correct:

1. Compute eigenvalues of A.
2. Pair complex-conjugate eigenvalues; keep real eigenvalues as
   singletons.
3. Build a real block-diagonal matrix N with a 1x1 block `[[r]]` for
   each real eigenvalue r, and a 2x2 block `[[a, b], [-b, a]]` for each
   conjugate pair a +/- ib.

Each 2x2 block has eigenvalues a +/- ib and satisfies block @ block.T ==
block.T @ block, so it is normal. A block-diagonal matrix built from
normal blocks is normal. The full spectrum is a exact match by
construction.

The eigenvector geometry is completely destroyed — which is the point.
The block-diagonal basis is arbitrary; only the spectrum is preserved.
"""

from __future__ import annotations

import numpy as np


def _pair_conjugates(eigs: np.ndarray, tol: float = 1e-8):
    """Return (reals, pairs) where pairs is a list of (a, b) with b > 0."""
    reals = []
    pairs = []
    remaining = list(eigs)
    while remaining:
        lam = remaining.pop(0)
        if abs(lam.imag) < tol:
            reals.append(lam.real)
            continue
        target = lam.conjugate()
        best = None
        best_err = np.inf
        for j, other in enumerate(remaining):
            err = abs(other - target)
            if err < best_err:
                best_err = err
                best = j
        if best is None or best_err > tol * max(1.0, abs(lam)):
            raise ValueError(
                f"eigenvalue {lam} has no conjugate partner within tol={tol}; "
                f"is A really real?"
            )
        remaining.pop(best)
        a = lam.real
        b = abs(lam.imag)
        pairs.append((a, b))
    return reals, pairs


def spectrum_matched_normal(A: np.ndarray, tol: float = 1e-8) -> np.ndarray:
    """Real normal matrix with the same eigenvalues as A."""
    A = np.asarray(A, dtype=np.float64)
    n = A.shape[0]
    eigs = np.linalg.eigvals(A)
    reals, pairs = _pair_conjugates(eigs, tol=tol)

    N = np.zeros((n, n), dtype=np.float64)
    k = 0
    for r in reals:
        N[k, k] = r
        k += 1
    for a, b in pairs:
        N[k, k] = a
        N[k, k + 1] = b
        N[k + 1, k] = -b
        N[k + 1, k + 1] = a
        k += 2
    assert k == n, f"block layout filled {k} of {n} rows"

    lam_src = np.sort_complex(np.linalg.eigvals(A))
    lam_null = np.sort_complex(np.linalg.eigvals(N))
    spec_err = float(np.max(np.abs(lam_src - lam_null)))
    if spec_err > tol:
        raise AssertionError(
            f"spectrum-match failed: max |dlambda| = {spec_err:.2e} > {tol:.1e}"
        )
    normality_err = float(np.max(np.abs(N @ N.T - N.T @ N)))
    if normality_err > tol:
        raise AssertionError(
            f"surrogate is not normal: max |[N, N^T]| = "
            f"{normality_err:.2e} > {tol:.1e}"
        )
    return N
