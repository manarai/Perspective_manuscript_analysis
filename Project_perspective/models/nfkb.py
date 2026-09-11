"""NF-kB - IkB negative feedback module -- amplifying / non-normal.

Source
------
Hoffmann et al. (2002) Science 298:1241; minimal-loop version follows
Nelson et al. (2004) Science 306:704 and Ashall et al. (2009) Science
324:242.

State (linearised, 3D minimal loop)
-----------------------------------
x = free nuclear NF-kB (relative to equilibrium)
m = IkB mRNA
y = IkB protein

Dynamics (linear around the equilibrium; feedback closes through y)
    dx/dt = -gamma_x * x - k * y
    dm/dt = alpha * x - gamma_m * m
    dy/dt = rho * m - gamma_y * y

Operating point
---------------
Around the tonic-signalling steady state, with all variables scaled by
their equilibrium magnitudes (so x, m, y are order 1 deviations). Rates
in 1/min, taken from the Hoffmann-family fits: NF-kB nuclear turnover
gamma_x ~ 0.1/min; IkB mRNA turnover gamma_m ~ 0.03/min; IkB protein
turnover gamma_y ~ 0.02/min. Loop coupling alpha = rho = 0.05 and
sequestration k = 0.10 place the fixed point in the damped-focus
regime (loop gain ~ 2.5e-4, well below the Hopf threshold ~ 7.8e-4 for
these decay rates).

Region
------
Cyclic negative feedback with three well-separated timescales -- the
Jacobian is not symmetric, and its eigenvector geometry differs from
its spectrum. Reactive: excursions in x can be transiently amplified
before the loop closes.
"""
from __future__ import annotations

import numpy as np

from ._base import Model


def build(gamma_x: float = 0.10, gamma_m: float = 0.03, gamma_y: float = 0.02,
          alpha: float = 0.05, rho: float = 0.05, k: float = 0.10) -> Model:
    A = np.array(
        [
            [-gamma_x, 0.0, -k],
            [alpha, -gamma_m, 0.0],
            [0.0, rho, -gamma_y],
        ]
    )
    return Model(
        key="nfkb",
        name="NF-kB -- IkB negative feedback",
        domain="transcriptional regulation",
        region="amplifying",
        A=A,
        x_eq=np.ones(3),
        params={"gamma_x": gamma_x, "gamma_m": gamma_m, "gamma_y": gamma_y,
                "alpha": alpha, "rho": rho, "k": k},
        citation="Hoffmann 2002; Nelson 2004; Ashall 2009",
        operating_point="linearisation at the tonic-signalling steady state",
        scaling="states scaled by equilibrium magnitude; rates in 1/min",
        notes=("Cyclic feedback with three timescales; non-symmetric Jacobian "
               "with eigenvector geometry distinct from spectrum."),
    )
