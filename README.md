# Operator Descriptors — code and results

Companion repository for *Operator Descriptors: A Comparative Language for
Biological Dynamics* (Caylor & Terooatea, Perspective, in revision).

All figures and tables in the manuscript are reproducible from the code
and CSVs here. No experimental data — every computation runs on analytic
Jacobians defined in `models/`.

## Layout

```
descriptors/continuous.py   Descriptor definitions: Q, ν, r_PR/n, ω, G_max.
models/                     Nine mechanistic models with analytic Jacobians
                            at stated operating points (TCR, NF-κB,
                            Murphy–Miller, Neubert–Caswell predator–prey,
                            repressilator, Goodwin, Selkov, toggle, relay).
                            Each build() returns a Model dataclass; see
                            models/_base.py.
nulls/spectrum_matched.py   Real block-normal surrogate: same spectrum as
                            A, orthogonal eigenvectors, ν = 0 by construction.
sweep.py                    Runs every model, computes descriptors for A and
                            for its spectrum-matched surrogate, and sweeps
                            per-state rescalings (0.5, 2.0, uniform(0.5, 2.0)).
                            Writes results/descriptors.csv and
                            results/sensitivity.csv, and prints a summary.
make_fig2.py                Builds Figure 2 (three panels) from the CSVs
                            and the analytic Jacobians. Writes figures/fig2.pdf
                            and figures/fig2.png.
results/                    descriptors.csv, sensitivity.csv (source of every
                            number reported in the manuscript).
figures/                    fig2.pdf, fig2.png.
tests/                      pytest suite for descriptor definitions.
response_to_reviewers.md    Response letter to Perspective-route reviewers.
```

## Reproduce

```
python sweep.py              # (re)computes results/*.csv
python make_fig2.py          # rebuilds figures/fig2.{pdf,png}
pytest -q                    # descriptor unit tests
```

Requires Python ≥ 3.10, `numpy`, `scipy`, `matplotlib`.

## The twelve computations

The panel in Figure 2 covers nine source systems in eight domains, with
three symmetry-artifact probes (repressilator, toggle, relay in symmetric
and heterogeneous variants):

| Model                        | n | Region                      | Source paper                              |
|------------------------------|---|-----------------------------|-------------------------------------------|
| TCR proofreading             | 5 | amplifying, non-oscillatory | McKeithan 1995                            |
| Neubert–Caswell predator–prey| 2 | amplifying, non-oscillatory | Neubert & Caswell 1997                    |
| Murphy–Miller balanced E–I   | 2 | amplifying & oscillating    | Murphy & Miller 2009                      |
| Selkov glycolysis            | 2 | amplifying & oscillating    | Selkov 1968                               |
| NF-κB / IκB                  | 3 | amplifying & oscillating    | Hoffmann et al. 2002                      |
| Repressilator (sym / het)    | 3 | oscillating                 | Elowitz & Leibler 2000                    |
| Goodwin oscillator           | 3 | oscillating, not amplifying | Goodwin 1965                              |
| Toggle switch (sym / het)    | 2 | neither                     | Gardner, Cantor & Collins 2000            |
| Relay chain (sym / het)      | 4 | neither                     | (diffusion–decay reference model)         |

## Descriptors

- **Q** = |Im λ₁| / |Re λ₁| — phase advance of the leading mode per unit
  decay time. Invariant under change of representation.
- **ν** = √(‖A‖_F² − Σ|λᵢ|²) / ‖A‖_F — Henrici departure from normality.
  Necessary precondition for transient amplification; metric-relative.
- **r_PR / n** — participation-ratio rank of the singular spectrum,
  divided by n. Effective rank.
- **ω** = λ_max((A + Aᵀ)/2) — numerical abscissa; worst-case initial
  growth rate. ω > 0 ⟺ G_max > 1 for stable A. Metric-relative.
- **G_max** = sup_{t ≥ 0} ‖e^{At}‖₂ — peak transient gain.

## Contact

tommy.terooatea@byu.edu
