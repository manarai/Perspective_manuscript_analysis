"""Goodwin oscillator -- oscillatory (damped-focus regime).

Source
------
Goodwin (1965) Adv. Enz. Regul. 3:425; canonical 3-variable
transcription -- translation -- repression loop. Same damped-focus
caveat as the repressilator: we linearise at the fixed point in a
regime where the Hopf has not yet fired.

Dynamics
--------
    dx/dt = V / (K^n + z^n) - a * x       (mRNA, repressed by z)
    dy/dt = b * x - c * y                  (protein)
    dz/dt = d * y - e * z                  (active/nuclear form)

Operating point -- explicit, damped-focus
-----------------------------------------
Choose all rates = 1 and V = K = 1. Symmetric steady state
x* = y* = z* satisfies z* + z*^{n+1} = 1. The Jacobian is

    A = [[-1, 0, -beta], [1, -1, 0], [0, 1, -1]]

with beta = V n z*^{n-1} / (K^n + z*^n)^2. Characteristic polynomial
(1 + lambda)^3 = -beta gives

    lambda_0 = -1 - beta^{1/3}
    lambda_{1,2} = -1 + beta^{1/3} * (1/2 +/- i sqrt(3)/2)

Hopf at beta = 8 (real part of lambda_{1,2} crosses zero). For n = 4 we
get beta ~ 0.98, safely below the Hopf and giving damped oscillations.

Region
------
Oscillatory: leading eigenvalue is complex with |Im|/|Re| ~ 1.7.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from ._base import Model


def build(n: float = 7.0, V: float = 1.0, K: float = 1.0) -> Model:
    # Solve z + z^{n+1} = V for z > 0 (with a=b=c=d=e=1, K=1).
    z_star = brentq(lambda z: z + z ** (n + 1) - V, 1e-6, V + 1.0)
    beta = V * n * z_star ** (n - 1) / (K ** n + z_star ** n) ** 2
    if beta >= 8.0:
        raise ValueError(
            f"beta = {beta:.3f} >= 8 => Hopf-unstable; "
            f"choose smaller n for damped-focus"
        )
    A = np.array(
        [
            [-1.0, 0.0, -beta],
            [1.0, -1.0, 0.0],
            [0.0, 1.0, -1.0],
        ]
    )
    S = np.diag([z_star, z_star, z_star])
    A = np.linalg.inv(S) @ A @ S
    return Model(
        key="goodwin",
        name="Goodwin oscillator (damped-focus)",
        domain="transcriptional",
        region="oscillatory",
        A=A,
        x_eq=np.array([z_star, z_star, z_star]),
        params={"n": n, "V": V, "K": K, "beta": beta, "z_star": z_star},
        citation="Goodwin 1965 Adv. Enz. Regul. 3:425",
        operating_point=(f"symmetric fixed point z*={z_star:.3f}; "
                         f"beta={beta:.3f} < 8 (stable focus)"),
        scaling="states scaled by fixed-point magnitude",
        notes=("Damped-focus regime; linearisation inside a limit cycle "
               "at higher n would put Re lambda > 0. Not done here."),
    )
