"""Export S2: metric-sensitivity distribution over 500 random per-state
diagonal probes.

Reads results/sensitivity_dist.csv (long-format: one row per
model x descriptor, carrying baseline, mean/median/p05/p95/min/max/std
across the 500 diagonal draws, plus the sign-flip fraction for
descriptors whose sign carries meaning) and writes:

    supp/S2_sensitivity.csv   twelve rows, per-descriptor columns for
                              baseline, p05, p95, min, max, sign_flip
    supp/S2_sensitivity.md    the same table rendered as markdown

Protocol note (in the md header): uniform scalings (0.5x I and 2x I) do
NOT probe metric dependence, since DAD^-1 = A for D = cI. They are
included in sweep.py as invariance controls only, and are asserted to
recover the baseline to machine precision. Everything below is over the
500 genuinely diagonal draws from Uniform(0.5, 2.0).
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC  = ROOT / "results" / "sensitivity_dist.csv"
OUT_CSV = ROOT / "supp" / "S2_sensitivity.csv"
OUT_MD  = ROOT / "supp" / "S2_sensitivity.md"

DISPLAY = {
    "tcr":               "TCR proofreading",
    "predator_prey":     "Neubert–Caswell predator–prey",
    "murphy_miller":     "Murphy–Miller balanced E–I",
    "selkov":            "Selkov glycolysis",
    "nfkb":              "NF-κB / IκB",
    "repressilator_het": "Repressilator (heterogeneous)",
    "repressilator":    "Repressilator (symmetric)",
    "goodwin":           "Goodwin",
    "toggle":            "Toggle (symmetric)",
    "toggle_het":        "Toggle (heterogeneous)",
    "relay":             "Relay (symmetric)",
    "relay_het":         "Relay (heterogeneous)",
}
ORDER = list(DISPLAY)
DESCRIPTORS = ["nu", "omega", "G_max"]


def load():
    by_model = {}
    with SRC.open() as f:
        for r in csv.DictReader(f):
            by_model.setdefault(r["model"], {})[r["descriptor"]] = {
                "baseline": float(r["baseline"]),
                "p05":  float(r["p05"]),
                "p95":  float(r["p95"]),
                "min":  float(r["min"]),
                "max":  float(r["max"]),
                "median": float(r["median"]),
                "frac_sign_flip":
                    float(r["frac_sign_flip"]) if r["frac_sign_flip"] != "" else None,
            }
    return by_model


def write_csv(rows):
    header = ["model", "display_name"]
    for d in DESCRIPTORS:
        for suffix in ("baseline", "p05", "p95", "min", "max", "frac_sign_flip"):
            header.append(f"{d}_{suffix}")
    with OUT_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for key in ORDER:
            r = rows[key]
            row = [key, DISPLAY[key]]
            for d in DESCRIPTORS:
                s = r[d]
                sf = "" if s["frac_sign_flip"] is None else f"{s['frac_sign_flip']:.3f}"
                row.extend([f"{s['baseline']:.4f}",
                            f"{s['p05']:.4f}", f"{s['p95']:.4f}",
                            f"{s['min']:.4f}", f"{s['max']:.4f}",
                            sf])
            w.writerow(row)


def fmt(v):
    return f"{v:+.3f}"


def write_md(rows):
    lines = []
    lines.append("# Supplementary Table S2 — Metric-sensitivity distribution")
    lines.append("")
    lines.append(
        "**Protocol.** For each of the twelve computations we draw 500 "
        "independent per-state rescalings s ∈ Uniform(0.5, 2.0)^n and "
        "recompute every descriptor at A' = diag(s)⁻¹ A diag(s). "
        "Uniform scalings D = c·I are known to leave DAD⁻¹ = A unchanged; "
        "they are kept in `sweep.py` as invariance controls only and "
        "verified numerically as an implementation check. They are NOT "
        "counted as sensitivity probes. Random seed 0."
    )
    lines.append("")
    lines.append(
        "**Reading.** Q is invariant under change of representation and "
        "is not tabulated. ν, ω and G_max are metric-relative. For ω, the "
        "*sign-flip fraction* is the proportion of the 500 draws in which "
        "the sign of ω flips relative to baseline — a direct measure of "
        "how close the assignment sits to the ω = 0 boundary. The two "
        "amber points in Figure 2a correspond to the two computations with "
        "sign-flip fraction above 50%."
    )
    lines.append("")

    lines.append(
        "| Model | ν baseline | ν p05–p95 | ω baseline | ω p05–p95 | "
        "ω sign-flip | G_max baseline | G_max p05–p95 |"
    )
    lines.append(
        "|---|---:|:---:|---:|:---:|---:|---:|:---:|"
    )

    for key in ORDER:
        r = rows[key]
        n = r["nu"]; om = r["omega"]; g = r["G_max"]
        sf = om["frac_sign_flip"]
        sf_txt = "—" if sf is None else f"{sf*100:.1f} %"
        lines.append(
            f"| {DISPLAY[key]} "
            f"| {n['baseline']:.3f} | {n['p05']:.3f} – {n['p95']:.3f} "
            f"| {fmt(om['baseline'])} | {fmt(om['p05'])} – {fmt(om['p95'])} "
            f"| {sf_txt} "
            f"| {g['baseline']:.3f} | {g['p05']:.3f} – {g['p95']:.3f} |"
        )

    lines.append("")
    lines.append(
        "**How the sign-flip column reads.** A value of 0.0 % means the "
        "amplifier assignment (or the non-amplifier assignment) holds on "
        "every one of the 500 draws — a robust sign statement about that "
        "system under diagonal metric warp. A value at or above 50 % means "
        "the assignment is decided as often by the metric as by the "
        "operator, and any placement in Figure 2a with that mark is "
        "reported as a marginal case rather than a claim. Values in "
        "between are quantified fragility, not anecdote."
    )
    lines.append("")
    lines.append("Source: `results/sensitivity_dist.csv`, generated by `sweep.py` "
                 "(N = 500 diagonal probes, seed = 0). The legacy min/max "
                 "envelope is retained in `results/sensitivity.csv` for "
                 "back-compatibility only.")
    lines.append("")
    OUT_MD.write_text("\n".join(lines))


def main():
    rows = load()
    write_csv(rows)
    write_md(rows)
    print(f"wrote {OUT_CSV.relative_to(ROOT)}")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
