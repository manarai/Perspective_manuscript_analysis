"""Common types for the zoo."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np


@dataclass
class Model:
    key: str                  # short id, e.g. "tcr"
    name: str                 # display name
    domain: str               # immune, transcriptional, neural, ...
    region: str               # amplifying | oscillatory | toggle | contractive
    A: np.ndarray             # Jacobian at operating point, nondim
    x_eq: np.ndarray          # equilibrium magnitudes used for scaling
    params: dict              # cited parameter values
    citation: str             # source paper
    operating_point: str      # description of the point A is taken at
    scaling: str              # nondimensionalisation rule used
    notes: str = ""           # e.g. Hopf choice, floquet caveat


ModelBuilder = Callable[[], Model]
