"""Balanced excitatory-inhibitory circuit -- amplifying / non-normal.

Source
------
Murphy & Miller (2009) Neuron 61:635. Two-population rate model in the
balanced regime; the canonical demonstration that non-normality
generates large transient responses in cortical circuits even at modest
spectral abscissa.

State
-----
E = excitatory population rate, I = inhibitory population rate
(dimensionless, in units of the linearised deviation from steady state).

Dynamics (already linear at threshold-linear expansion around the fixed
point) -- Murphy & Miller Eq. 1 with the "balanced" convention that the
E->I and E->E weights are similarly large and offset by strong I->E:

    dE/dt = (-1 + w_EE) * E - w_EI * I
    dI/dt = w_IE * E + (-1 - w_II) * I     (self-inhibition subtracts)

Operating point
---------------
Around the balanced steady state at unit rate. All weights positive; the
"balance" is between w_EE and w_IE and between w_EI and w_II. Values
from Murphy & Miller Fig 2 example: w_EE = 4.7, w_EI = 4.5, w_IE = 5.0,
w_II = 4.0. This regime is stable (both eigenvalues negative) but
strongly non-normal (E and I nearly cancel at steady state; small
perturbations expand along the (E, -I) direction before contracting).

Region
------
The paradigmatic amplifying operator: modest, stable spectrum with
G_max often > 10.
"""
from __future__ import annotations

import numpy as np

from ._base import Model


def build(w_EE: float = 4.7, w_EI: float = 4.5,
          w_IE: float = 5.0, w_II: float = 4.0) -> Model:
    A = np.array(
        [
            [-1.0 + w_EE, -w_EI],
            [w_IE, -1.0 - w_II],
        ]
    )
    return Model(
        key="murphy_miller",
        name="Balanced E-I circuit",
        domain="neural",
        region="amplifying",
        A=A,
        x_eq=np.ones(2),
        params={"w_EE": w_EE, "w_EI": w_EI, "w_IE": w_IE, "w_II": w_II},
        citation="Murphy & Miller 2009 Neuron 61:635",
        operating_point="balanced-regime steady state, unit-rate scaling",
        scaling="rates in 1/tau_m; weights dimensionless",
        notes=("Reactive under any reasonable metric; small changes in "
               "w_IE-w_EE move numerical abscissa and G_max sharply."),
    )
