"""Compute descriptors, surrogates, and metric-sensitivity distributions.

Protocol (v2 — expanded after the referee-side observation that uniform
scalings DAD^-1 with D = cI leave A unchanged and therefore do not probe
metric dependence at all).

The sweep now does three things per model:

1. Baseline descriptors at the stated operating point.
2. Two INVARIANCE CONTROLS — uniform 0.5 x I and 2.0 x I. These recover
   the baseline to machine precision. They are kept in the sweep and
   asserted as sanity checks on the implementation; they are not counted
   as sensitivity probes.
3. A distribution of GENUINE per-state random probes. Diagonal entries
   are drawn independently from Uniform(0.5, 2.0) and applied as A ->
   diag(s)^{-1} A diag(s). N_PROBES = 500 with a fixed seed.

Outputs
-------
results/descriptors.csv           model | region | A / A_null descriptors
                                  (unchanged from v1)
results/sensitivity_dist.csv      model | descriptor | baseline | mean |
                                  median | p05 | p95 | frac_sign_flip
                                  (frac_sign_flip is meaningful only for
                                   descriptors that can change sign,
                                   populated for omega and alpha).
results/sensitivity.csv           kept as a compatibility shim: baseline
                                  and full min / max over ALL probes (500
                                  random + 2 invariance controls). This is
                                  the file existing tooling reads.
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
N_PROBES = 500
SEED = 0
RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def _rescaled_A(A: np.ndarray, scale: np.ndarray) -> np.ndarray:
    S = np.diag(scale)
    return np.linalg.inv(S) @ A @ S


def descriptor_row(m, A: np.ndarray) -> dict:
    d = compute(A).as_dict()
    d["model"] = m.key
    d["region"] = m.region
    return d


def run() -> None:
    models = all_models()
    rng = np.random.default_rng(SEED)

    # ---- (1) baseline + spectrum-matched surrogate ------------------------
    rows = []
    for m in models:
        d_A = descriptor_row(m, m.A); d_A["variant"] = "A"
        A_null = spectrum_matched_normal(m.A)
        d_null = descriptor_row(m, A_null); d_null["variant"] = "A_null"
        rows.append(d_A); rows.append(d_null)

    (RESULTS_DIR / "descriptors.csv").open("w", newline="").close()
    with (RESULTS_DIR / "descriptors.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["model", "region", "variant"] + DESCRIPTORS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in w.fieldnames})

    # ---- (2) invariance controls (must recover baseline) ------------------
    for m in models:
        n = m.A.shape[0]
        base = compute(m.A)
        for c in (0.5, 2.0):
            d = compute(_rescaled_A(m.A, np.full(n, c)))
            for k in ("nu", "omega", "G_max"):
                assert np.isclose(getattr(d, k), getattr(base, k), atol=1e-8), \
                    (f"invariance control failed: {m.key} at uniform {c}x, "
                     f"{k}: {getattr(d, k)} vs baseline {getattr(base, k)}")

    # ---- (3) genuine per-state random probes ------------------------------
    dist_rows = []          # long format: model, descriptor, N stats
    legacy_rows = []        # min/max over the full envelope (for back-compat)

    for m in models:
        n = m.A.shape[0]
        base = compute(m.A)

        # Random diagonal draws from Uniform(0.5, 2.0), independent per state.
        vals = {k: [] for k in DESCRIPTORS}
        sign_flip_counts = {k: 0 for k in ("alpha", "omega")}
        for _ in range(N_PROBES):
            s = rng.uniform(0.5, 2.0, size=n)
            d = compute(_rescaled_A(m.A, s))
            for k in DESCRIPTORS:
                vals[k].append(getattr(d, k))
            for k in ("alpha", "omega"):
                base_v = getattr(base, k)
                probe_v = getattr(d, k)
                if np.sign(base_v) != np.sign(probe_v) and \
                        abs(base_v) > 1e-12 and abs(probe_v) > 1e-12:
                    sign_flip_counts[k] += 1

        for k in DESCRIPTORS:
            arr = np.array(vals[k], dtype=float)
            baseline_v = getattr(base, k)
            dist_rows.append({
                "model": m.key,
                "descriptor": k,
                "baseline": baseline_v,
                "n_probes": N_PROBES,
                "mean": float(np.nanmean(arr)),
                "median": float(np.nanmedian(arr)),
                "p05": float(np.nanpercentile(arr, 5)),
                "p95": float(np.nanpercentile(arr, 95)),
                "min": float(np.nanmin(arr)),
                "max": float(np.nanmax(arr)),
                "std": float(np.nanstd(arr)),
                "frac_sign_flip": (sign_flip_counts[k] / N_PROBES
                                   if k in sign_flip_counts else ""),
            })
            # Legacy envelope also includes the invariance controls, which
            # under a proper diagonal probe equal baseline exactly; so min /
            # max over the union is min / max of the random probes bracketed
            # by baseline. Report that so anything still reading
            # sensitivity.csv gets the widest observed range.
            legacy_rows.append({
                "model": m.key,
                "descriptor": k,
                "baseline": baseline_v,
                "min_over_probes": float(min(np.nanmin(arr), baseline_v)),
                "max_over_probes": float(max(np.nanmax(arr), baseline_v)),
                "std_over_probes": float(np.nanstd(arr)),
            })

    with (RESULTS_DIR / "sensitivity_dist.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "model", "descriptor", "baseline", "n_probes",
            "mean", "median", "p05", "p95",
            "min", "max", "std", "frac_sign_flip",
        ])
        w.writeheader()
        w.writerows(dist_rows)

    with (RESULTS_DIR / "sensitivity.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "model", "descriptor", "baseline",
            "min_over_scalings", "max_over_scalings", "std_over_scalings",
        ])
        w.writeheader()
        for r in legacy_rows:
            w.writerow({
                "model": r["model"],
                "descriptor": r["descriptor"],
                "baseline": r["baseline"],
                "min_over_scalings": r["min_over_probes"],
                "max_over_scalings": r["max_over_probes"],
                "std_over_scalings": r["std_over_probes"],
            })

    _print_summary(rows, dist_rows)


def _print_summary(rows: list[dict], dist_rows: list[dict]) -> None:
    print()
    print("=" * 100)
    print(f"METRIC-SENSITIVITY DISTRIBUTION ({N_PROBES} random diagonal probes, "
          f"seed={SEED})")
    print("=" * 100)
    header = (f"{'model':16s} {'descriptor':8s} {'baseline':>10s} "
              f"{'p05':>10s} {'median':>10s} {'p95':>10s} "
              f"{'sign_flip':>10s}")
    print(header)
    print("-" * len(header))
    show = ["nu", "omega", "G_max"]
    for r in dist_rows:
        if r["descriptor"] not in show:
            continue
        sf = ""
        if r["frac_sign_flip"] != "":
            sf = f"{r['frac_sign_flip']*100:>7.1f} %"
        print(f"{r['model']:16s} {r['descriptor']:8s} "
              f"{r['baseline']:+10.4f} "
              f"{r['p05']:+10.4f} {r['median']:+10.4f} {r['p95']:+10.4f} "
              f"{sf:>10s}")


if __name__ == "__main__":
    run()
