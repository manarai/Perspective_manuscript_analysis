"""Reactive predator-prey system -- amplifying / non-normal.

Source
------
Neubert & Caswell (1997) Ecology 78:653. Introduced the reactivity /
numerical-abscissa concept in ecology using this exact system.

Dynamics
--------
    dx/dt = r * x * (1 - x/K) - a * x * y     (prey)
    dy/dt = b * a * x * y - m * y              (predator)

Operating point
---------------
Coexistence equilibrium (x*, y*) with x* = m/(b*a) and
y* = (r/a)(1 - x*/K). Linearisation gives

    J = [[ r - 2 r x*/K - a y* ,   -a x* ],
         [       b a y*        , b a x* - m ]]

Parameters (Neubert & Caswell Fig 3 example, unitised): r = 1, K = 10,
a = 0.1, b = 0.5, m = 0.4. That places (x*, y*) at (8, 2). After scaling
each state by its equilibrium magnitude the Jacobian is unchanged in
sign structure and gives a stable node with strictly positive numerical
abscissa -- Neubert & Caswell's headline result.

Region
------
Amplifying: eigenvalues both real-negative, but omega > 0 so a small
perturbation grows before decaying. The system loses contractivity
without losing spectral stability.
"""
from __future__ import annotations

import numpy as np

from ._base import Model


def build(r: float = 1.0, K: float = 10.0, a: float = 0.1,
          b: float = 0.5, m: float = 0.4) -> Model:
    x_eq = m / (b * a)
    y_eq = (r / a) * (1.0 - x_eq / K)
    if x_eq <= 0 or y_eq <= 0:
        raise ValueError(
            f"parameters produce nonpositive coexistence "
            f"equilibrium: x*={x_eq}, y*={y_eq}"
        )
    J = np.array(
        [
            [r - 2 * r * x_eq / K - a * y_eq, -a * x_eq],
            [b * a * y_eq, b * a * x_eq - m],
        ]
    )
    S = np.diag([x_eq, y_eq])
    A = np.linalg.inv(S) @ J @ S               # scale each state by its equilibrium
    return Model(
        key="predator_prey",
        name="Neubert-Caswell predator-prey",
        domain="ecology",
        region="amplifying",
        A=A,
        x_eq=np.array([x_eq, y_eq]),
        params={"r": r, "K": K, "a": a, "b": b, "m": m},
        citation="Neubert & Caswell 1997 Ecology 78:653",
        operating_point="stable coexistence equilibrium",
        scaling="each state scaled by equilibrium magnitude",
        notes=("Original example for the reactivity concept; expect "
               "numerical abscissa > 0 with negative spectral abscissa."),
    )
