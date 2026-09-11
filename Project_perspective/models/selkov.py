"""Selkov glycolytic oscillator -- oscillatory (metabolic, non-transcriptional).

Source
------
Selkov (1968) Eur. J. Biochem. 4:79. Canonical minimal model of
autocatalytic glycolytic oscillations, entirely metabolic. Included in
the zoo specifically to give the oscillatory region a representative
that is NOT a delayed transcriptional feedback loop; the other three
oscillators (NF-kB, repressilator, Goodwin) all rely on the same
mechanism family (transcription factor represses its own promoter with
a translation-and-nuclear-transport delay).

Dynamics
--------
    dx/dt = -x + a * y + x^2 * y      (x = ADP)
    dy/dt = b - a * y - x^2 * y       (y = fructose-6-phosphate)

The nonlinear term x^2 y is autocatalytic product activation of
phosphofructokinase -- a biochemically distinct mechanism from
transcriptional feedback delays.

Operating point -- explicit, damped-focus, NOT the limit cycle
--------------------------------------------------------------
Fixed point: x* = b, y* = b / (a + b^2). Jacobian at the fixed point:

    J11 = (b^2 - a) / (a + b^2)     J12 = a + b^2
    J21 = -2 b^2 / (a + b^2)        J22 = -(a + b^2)

Hopf when Tr(J) = 0. For a = 0.1, b = 0.3 the fixed point is a damped
focus (Tr < 0, discriminant negative) with clear Q > 1. Above the Hopf
the fixed point is unstable inside a limit cycle -- same caveat as the
repressilator; we linearise deliberately in the sub-Hopf regime.

Region
------
Oscillatory: complex-conjugate leading eigenpair with |Im| / |Re| ~ 3.5.
"""
from __future__ import annotations

import numpy as np

from ._base import Model


def build(a: float = 0.1, b: float = 0.3) -> Model:
    x_star = b
    y_star = b / (a + b ** 2)
    Kab = a + b ** 2
    J11 = (b ** 2 - a) / Kab
    J12 = Kab
    J21 = -2 * b ** 2 / Kab
    J22 = -Kab
    J = np.array([[J11, J12], [J21, J22]])
    if np.max(np.linalg.eigvals(J).real) >= 0:
        raise ValueError(
            f"parameters a={a}, b={b} give unstable fixed point "
            f"(above Hopf); pick a, b in the damped-focus regime"
        )
    S = np.diag([x_star, y_star])
    A = np.linalg.inv(S) @ J @ S
    return Model(
        key="selkov",
        name="Selkov glycolytic oscillator (damped-focus)",
        domain="metabolism",
        region="oscillatory",
        A=A,
        x_eq=np.array([x_star, y_star]),
        params={"a": a, "b": b, "x_star": x_star, "y_star": y_star},
        citation="Selkov 1968 Eur. J. Biochem. 4:79",
        operating_point=(f"stable focus at (x*={x_star:.3f}, y*={y_star:.3f}); "
                         f"sub-Hopf, autocatalytic damped oscillation"),
        scaling="each state scaled by fixed-point magnitude",
        notes=("Metabolic autocatalysis, NOT delayed transcriptional feedback. "
               "Included specifically so the oscillatory region has a "
               "representative from a mechanism family other than "
               "transcription-factor negative-feedback loops."),
    )
