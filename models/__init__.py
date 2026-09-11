"""Model zoo -- 8 published mechanistic systems.

Each build() returns a Model with an analytic Jacobian at a stated
operating point, in a stated nondimensionalisation. Region labels are
the paper's target regions; the point of the illustration is to
demonstrate that unrelated systems in the same region co-locate in the
descriptor plane despite sharing no molecular components, network
topology, or timescale.
"""
from __future__ import annotations

from ._base import Model
from . import (
    tcr,
    nfkb,
    murphy_miller,
    predator_prey,
    repressilator,
    goodwin,
    selkov,
    toggle,
    relay,
)


def all_models() -> list[Model]:
    return [
        tcr.build(),
        nfkb.build(),
        murphy_miller.build(),
        predator_prey.build(),
        repressilator.build(),                     # symmetric, nu = 0 by construction
        repressilator.build(decay_spread=0.6),     # heterogeneous, nu > 0
        goodwin.build(),
        selkov.build(),                            # metabolic oscillator -- distinct mechanism
        toggle.build(),
        toggle.build(asymmetry=0.2),               # heterogeneous, symmetry-artifact probe
        relay.build(),
        relay.build(decay_spread=0.3),             # heterogeneous, symmetry-artifact probe
    ]
