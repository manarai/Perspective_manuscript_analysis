"""Compute descriptors, surrogates, and scaling sensitivity for the 8 models.

Follows Perspective_Route_Revision_Spec.md §1 order-of-work step 5:
"Run the surrogate null. Look at the result before writing anything --
if the amplifying models do *not* separate from their surrogates,
that changes what the paper can claim."

Outputs
-------
results/descriptors.csv         model | region | A-descriptors | A_null-descriptors
results/sensitivity.csv         model | descriptor | mean | std across scalings
Prints a colour-labelled summary table to stdout so the go / no-go call
can be made immediately.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from descriptors.continuous import compute
from models import all_models
from nulls.spectrum_matched import spectrum_matched_normal


DESCRIPTORS = ["alpha", "tau", "gap", "Q", "r_PR", "omega", "G_max", "t_star",
               "nu", "horizon"]
RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def _rescaled_A(A: np.ndarray, x_eq: np.ndarray, scale: np.ndarray) -> np.ndarray:
    """Apply a per-state rescaling s_i to the Jacobian.

    If the original states are already scaled to unit magnitude at the
    equilibrium, multiplying s_i corresponds to using x_eq_i * s_i as
    the characteristic magnitude of state i instead. A -> S^{-1} A S.
    """
    S = np.diag(scale)
    return np.linalg.inv(S) @ A @ S


def descriptor_row(m, A: np.ndarray) -> dict:
    d = compute(A).as_dict()
    d["model"] = m.key
    d["region"] = m.region
    return d


def run() -> None:
    models = all_models()

    rows = []
    for m in models:
        d_A = descriptor_row(m, m.A); d_A["variant"] = "A"
        A_null = spectrum_matched_normal(m.A)
        d_null = descriptor_row(m, A_null); d_null["variant"] = "A_null"
        rows.append(d_A); rows.append(d_null)

    csv_path = RESULTS_DIR / "descriptors.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["model", "region", "variant"] + DESCRIPTORS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in w.fieldnames})

    # Sensitivity sweep: 3 alternative per-state scalings per model.
    rng = np.random.default_rng(0)
    sens_rows = []
    for m in models:
        n = m.A.shape[0]
        scalings = [
            np.full(n, 0.5),
            np.full(n, 2.0),
            rng.uniform(0.5, 2.0, size=n),
        ]
        vals = {k: [] for k in DESCRIPTORS}
        for s in scalings:
            d = compute(_rescaled_A(m.A, m.x_eq, s))
            for k in DESCRIPTORS:
                vals[k].append(getattr(d, k))
        baseline = compute(m.A)
        for k in DESCRIPTORS:
            arr = np.array(vals[k], dtype=float)
            sens_rows.append({
                "model": m.key,
                "descriptor": k,
                "baseline": getattr(baseline, k),
                "min_over_scalings": float(np.nanmin(arr)),
                "max_over_scalings": float(np.nanmax(arr)),
                "std_over_scalings": float(np.nanstd(arr)),
            })

    sens_path = RESULTS_DIR / "sensitivity.csv"
    with sens_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["model", "descriptor", "baseline",
                                          "min_over_scalings",
                                          "max_over_scalings",
                                          "std_over_scalings"])
        w.writeheader()
        w.writerows(sens_rows)

    _print_summary(rows, sens_rows)


def _print_summary(rows: list[dict], sens_rows: list[dict]) -> None:
    key_cols = ["alpha", "omega", "nu", "G_max", "Q", "r_PR"]
    print()
    print("=" * 95)
    print("DESCRIPTORS OF A AND OF THE SPECTRUM-MATCHED NORMAL SURROGATE")
    print("=" * 95)
    header = f"{'model':16s} {'region':12s} {'variant':7s} " + \
             " ".join(f"{c:>9s}" for c in key_cols)
    print(header)
    print("-" * len(header))
    by_model = {}
    for r in rows:
        by_model.setdefault(r["model"], {})[r["variant"]] = r
    for mkey, variants in by_model.items():
        for vname in ("A", "A_null"):
            r = variants[vname]
            row = f"{mkey:16s} {r['region']:12s} {vname:7s} " + \
                  " ".join(f"{r[c]:+9.3f}" if isinstance(r[c], float) else f"{r[c]:>9}" for c in key_cols)
            print(row)
        # Separation metric: how much do the geometry descriptors move
        # between A and its spectrum-matched surrogate?
        dnu = variants["A"]["nu"] - variants["A_null"]["nu"]
        dG  = variants["A"]["G_max"] / max(1e-12, variants["A_null"]["G_max"])
        print(f"{'':16s} {'':12s} {'delta':7s} "
              f"{'':>9s} {'':>9s} {dnu:+9.3f} {dG:>9.2f}x")
        print()

    print()
    print("=" * 95)
    print("SCALING SENSITIVITY (per-state rescaling by {0.5, 2.0, U(0.5,2.0)})")
    print("=" * 95)
    show = ["nu", "omega", "G_max"]
    header = f"{'model':16s} {'descriptor':10s} {'baseline':>10s} " \
             f"{'min':>10s} {'max':>10s} {'span/|base|':>12s}"
    print(header)
    print("-" * len(header))
    for r in sens_rows:
        if r["descriptor"] not in show:
            continue
        base = r["baseline"]
        span = r["max_over_scalings"] - r["min_over_scalings"]
        rel = span / (abs(base) + 1e-12)
        print(f"{r['model']:16s} {r['descriptor']:10s} {base:+10.4f} "
              f"{r['min_over_scalings']:+10.4f} {r['max_over_scalings']:+10.4f} "
              f"{rel:>12.2f}")


if __name__ == "__main__":
    run()
