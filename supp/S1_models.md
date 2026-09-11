# Supplementary Table S1 — Model equations, reductions, and source-paper labels

**Scope.** Every computation reported in Figure 2 and in the response letter
runs on the analytic Jacobian of a mechanistic ODE, evaluated at a stated
operating point, in a stated nondimensionalisation. This supplement records
each model in full and links it back to its source paper. Numerical values
of the Jacobian are the `A` matrix in `models/<key>.py`; parameter values
are the defaults in each `build()` call, printed here.

**Label provenance.** Prediction P1 promises that the functional label used
in placement is not something we imposed. For each model, the last three
rows give the exact phrase the source paper uses to describe the system's
behaviour, and where it appears. These entries are left blank in this
draft — the author will fill them from the source papers before submission,
so the assignment is auditable rather than paraphrased from memory.

---

## 1. TCR kinetic proofreading (`tcr`)

- **Domain**: immune signalling
- **Dimension**: n = 5
- **Region (this paper)**: amplifying, non-oscillatory
- **Source**: McKeithan (1995) *PNAS* 92:5042. Structure reused by
  Ganguli, Huh & Sompolinsky (2008) *PNAS* 105:18970.

**Dynamics.** Sequential phosphorylation states C₁,…,C_N of the TCR-pMHC
complex, driven by ligand L:

    dC₁/dt = k_on·L − (k_off + k_p)·C₁
    dC_i/dt = k_p·C_{i−1} − (k_off + k_p)·C_i     (2 ≤ i < N)
    dC_N/dt = k_p·C_{N−1} − k_off·C_N

**Reduction.** L drops out under linearisation; the intra-cascade Jacobian
is subdiagonal-only and does not depend on L. Nondimensionalisation
unnecessary — the cascade is intrinsically dimensionless (rates in 1/s).

**Parameters used.** N = 5, k_off = 0.05 s⁻¹, k_p = 0.5 s⁻¹.

**Operating point.** L-driven cascade equilibrium.

**Source-paper functional label**: `_____`
**Verbatim phrase**: `_____`
**Location in source**: `_____`

---

## 2. NF-κB / IκB negative feedback (`nfkb`)

- **Domain**: transcriptional regulation
- **Dimension**: n = 3
- **Region**: amplifying (and oscillating at Q = 3.86)
- **Source**: Hoffmann et al. (2002) *Science* 298:1241; minimal-loop
  form follows Nelson et al. (2004) *Science* 306:704 and Ashall et al.
  (2009) *Science* 324:242.

**Dynamics.** Three-variable linearisation around the tonic-signalling
steady state (x = free nuclear NF-κB, m = IκB mRNA, y = IκB protein):

    dx/dt = −γ_x·x − k·y
    dm/dt = α·x − γ_m·m
    dy/dt = ρ·m − γ_y·y

**Reduction.** Minimal-loop three-variable form; the biochemically full
Hoffmann model has additional species (IKK, IκBβ/ε) omitted here.
Linearisation at the equilibrium; every state scaled by its equilibrium
magnitude so entries of A are order 1 in the deviation coordinates.

**Parameters used.** γ_x = 0.10/min, γ_m = 0.03/min, γ_y = 0.02/min,
α = 0.05/min, ρ = 0.05/min, k = 0.10/min. Loop gain ≈ 2.5·10⁻⁴, safely
below the Hopf threshold (~7.8·10⁻⁴ for these decay rates); the operating
point is a damped focus.

**Source-paper functional label**: `_____`
**Verbatim phrase**: `_____`
**Location in source**: `_____`

---

## 3. Murphy–Miller balanced excitation–inhibition (`murphy_miller`)

- **Domain**: neural
- **Dimension**: n = 2
- **Region**: amplifying and oscillating
- **Source**: Murphy & Miller (2009) *Neuron* 61:635.

**Dynamics.** Two-population rate model, threshold-linear at the balanced
steady state:

    dE/dt = (−1 + w_EE)·E − w_EI·I
    dI/dt =        w_IE·E + (−1 − w_II)·I

**Reduction.** Linearisation at the balanced steady state at unit rate;
weights are dimensionless, rates in 1/τ_m.

**Parameters used.** w_EE = 4.7, w_EI = 4.5, w_IE = 5.0, w_II = 4.0
(Murphy & Miller Fig 2 example).

**Source-paper functional label**: `_____`
**Verbatim phrase**: `_____`
**Location in source**: `_____`

---

## 4. Neubert–Caswell reactive predator–prey (`predator_prey`)

- **Domain**: ecology
- **Dimension**: n = 2
- **Region**: amplifying, non-oscillatory
- **Source**: Neubert & Caswell (1997) *Ecology* 78:653.

**Dynamics.** Logistic prey with linear functional response:

    dx/dt = r·x·(1 − x/K) − a·x·y
    dy/dt = b·a·x·y − m·y

**Reduction.** Linearisation at the coexistence equilibrium
(x* = m/(b·a), y* = (r/a)(1 − x*/K)). Each state scaled by its
equilibrium magnitude before descriptor computation.

**Parameters used.** r = 1, K = 10, a = 0.1, b = 0.5, m = 0.4;
(x*, y*) = (8, 2). Neubert & Caswell Fig 3 example.

**Source-paper functional label**: `_____`
**Verbatim phrase**: `_____`
**Location in source**: `_____`

---

## 5. Repressilator, symmetric (`repressilator`)

- **Domain**: synthetic circuit
- **Dimension**: n = 3
- **Region**: oscillating, not amplifying
- **Source**: Elowitz & Leibler (2000) *Nature* 403:335.

**Dynamics.** Cyclic Hill repression of three genes:

    dx_i/dt = α / (1 + x_{i−1}ⁿ) − x_i,    i = 1..3 cyclic

**Reduction.** Linearisation at the *symmetric fixed point*, deliberately
in the **sub-Hopf, damped-focus regime** (β < 2). Not the limit cycle:
above the Hopf the fixed point is unstable and its linearisation
describes escape rather than the orbit, and Floquet analysis on the
orbit produces a different object; the manuscript notes this
explicitly.

**Parameters used.** α = 5.0, n = 2.0 → x* ≈ 1.516, β ≈ 1.394.

**Symmetry-artifact note.** The circulant three-variable Jacobian at the
symmetric fixed point is *exactly normal* (ν = 0) — a property of that
reduction, not of oscillators. The `repressilator_het` entry probes this.

**Source-paper functional label**: `_____`
**Verbatim phrase**: `_____`
**Location in source**: `_____`

---

## 6. Repressilator, heterogeneous (`repressilator_het`)

- **Same source and dynamics as (5)**, with per-gene decay rates
  (1 − s, 1, 1 + s) for spread s = 0.6. The Jacobian is no longer
  circulant, so ν > 0 (0.27), while the oscillation quality Q is
  preserved (4.35 vs 3.98 in the symmetric case).
- Included as a symmetry-artifact probe rather than as an independent
  source; no additional label provenance to record.

---

## 7. Goodwin oscillator (`goodwin`)

- **Domain**: transcriptional
- **Dimension**: n = 3
- **Region**: oscillating, not amplifying
- **Source**: Goodwin (1965) *Adv. Enzyme Regul.* 3:425.

**Dynamics.** Transcription–translation–repression loop:

    dx/dt = V / (Kⁿ + zⁿ) − a·x        (mRNA, repressed by z)
    dy/dt = b·x − c·y                    (protein)
    dz/dt = d·y − e·z                    (active/nuclear form)

**Reduction.** Symmetric steady state z* + z*^{n+1} = V with a = b =
c = d = e = 1 and V = K = 1. **Sub-Hopf**, damped focus: β = V·n·z*^{n−1}
/ (Kⁿ + z*ⁿ)² < 8 at n = 7 (β ≈ 1.32).

**Parameters used.** n = 7, V = 1, K = 1 → z* ≈ 0.812, β ≈ 1.32.

**Source-paper functional label**: `_____`
**Verbatim phrase**: `_____`
**Location in source**: `_____`

---

## 8. Selkov glycolytic oscillator (`selkov`)

- **Domain**: metabolism
- **Dimension**: n = 2
- **Region**: amplifying and oscillating
- **Source**: Selkov (1968) *Eur. J. Biochem.* 4:79.

**Dynamics.** Autocatalytic phosphofructokinase activation:

    dx/dt = −x + a·y + x²·y      (ADP)
    dy/dt =  b − a·y − x²·y      (F6P)

**Reduction.** Fixed point x* = b, y* = b / (a + b²). **Sub-Hopf**,
damped focus (Tr < 0 with negative discriminant); above the Hopf the
fixed point is unstable inside a limit cycle. Each state scaled by its
equilibrium magnitude before descriptor computation.

**Parameters used.** a = 0.1, b = 0.3 → x* = 0.30, y* ≈ 1.58.

**Why included.** NF-κB, the repressilator and Goodwin are three
parameterisations of one mechanism family (delayed negative-feedback
transcription). Selkov gives the oscillating region a biochemically
distinct representative — autocatalytic product activation of PFK.

**Source-paper functional label**: `_____`
**Verbatim phrase**: `_____`
**Location in source**: `_____`

---

## 9. Toggle switch, symmetric (`toggle`)

- **Domain**: haematopoiesis (PU.1 / GATA1)
- **Dimension**: n = 2
- **Region**: neither (contractive at the stable attractor)
- **Source**: Huang, Guo, May & Enver (2007) *Dev. Biol.* 305:695;
  parameterisation follows Chickarmane, Enver & Peterson (2009)
  *PLoS Comput. Biol.* 5:e1000268.

**Dynamics.** Mutual repression with self-activation:

    dx/dt = α·xⁿ/(Kⁿ + xⁿ) + β/(Kⁿ + yⁿ) − γ·x
    dy/dt = α·yⁿ/(Kⁿ + yⁿ) + β/(Kⁿ + xⁿ) − γ·y

**Reduction.** Linearisation at the **stable PU.1-hi/GATA1-lo attractor**,
not at the saddle between attractors. This choice is stated as a limit
of the illustration in the main text: the decision structure lives at
the saddle, not at either stable attractor, so the low-effective-rank
axis of the framework is not exercised by the linearisations reported
here.

**Parameters used.** α = β = γ = 1, n = 4, K = 0.5. At the symmetric
attractor A is coincidentally symmetric under this scaling, giving
ν = 0.

**Source-paper functional label**: `_____`
**Verbatim phrase**: `_____`
**Location in source**: `_____`

---

## 10. Toggle switch, heterogeneous (`toggle_het`)

- **Same source and dynamics as (9)**, with the GATA1-side parameters
  scaled by (1 − s) for s = 0.2. Broken exchange symmetry gives
  ν = 0.19 while ω remains at −0.74 and G_max at 1.00 — the ν = 0
  assignment was a symmetry artifact, the contractivity assignment was
  not. Probe entry, no independent label provenance.

---

## 11. Diffusion–decay relay, symmetric (`relay`)

- **Domain**: signalling (relay)
- **Dimension**: n = 4
- **Region**: neither (contractive)
- **Source**: reference model, standard reaction–diffusion form; see
  Murray, *Mathematical Biology Vol. I*, Springer 2002. Not a specific
  published biological circuit.

**Dynamics.** Chain of N nodes with symmetric nearest-neighbour
diffusion and per-node decay (Neumann boundaries):

    dx_i/dt = D·(x_{i−1} − 2·x_i + x_{i+1}) − γ·x_i,   i = 1..N

**Reduction.** Linear model — the operator is the Jacobian.

**Parameters used.** N = 4, D = 0.5, γ = 0.2.

**Symmetry-artifact note.** A per-node decay perturbation preserves the
symmetry of A (a diagonal shift), so ν stays at 0 under this probe.
Breaking non-normality would require directional coupling
(k_forward ≠ k_backward), which turns the relay into a cascade — a
different mechanism, already represented by TCR.

**Source-paper functional label**: `_____` *(reference model; not a
published biological finding)*

---

## 12. Diffusion–decay relay, heterogeneous (`relay_het`)

- **Same as (11)** with per-node decay spread 0.3. As above, ν stays
  ≈ 0 because the perturbation is diagonal; this is reported as a
  finding, not concealed.

---

## Summary

| Region                       | Model                    | n | ω (baseline) | Q     | ν      | G_max  |
|------------------------------|--------------------------|---|--------------|-------|--------|--------|
| Amplifying, non-oscillatory  | TCR                      | 5 | +0.075       | 0.00  | 0.672  | 1.302  |
| Amplifying, non-oscillatory  | Neubert–Caswell          | 2 | +0.012       | 0.00  | 0.655  | 1.007  |
| Amplifying & oscillating     | Murphy–Miller            | 2 | +3.707       | 2.91  | 0.951  | 2.988  |
| Amplifying & oscillating     | Selkov                   | 2 | +0.294       | 3.45  | 0.803  | 1.616  |
| Amplifying & oscillating     | NF-κB                    | 3 | +0.009       | 3.86  | 0.466  | 1.072  |
| Amplifying & oscillating     | Repressilator (het)      | 3 | +0.095       | 4.35  | 0.270  | 1.040  |
| Oscillating, not amplifying  | Repressilator (sym)      | 3 | −0.303       | 3.98  | 0.000  | 1.000  |
| Oscillating, not amplifying  | Goodwin                  | 3 | −0.341       | 2.10  | 0.140  | 1.000  |
| Neither                      | Toggle (sym)             | 2 | −0.876       | 0.00  | 0.000  | 1.000  |
| Neither                      | Toggle (het)             | 2 | −0.740       | 0.00  | 0.192  | 1.000  |
| Neither                      | Relay (sym)              | 4 | −0.200       | 0.00  | 0.000  | 1.000  |
| Neither                      | Relay (het)              | 4 | −0.193       | 0.00  | 0.000  | 1.000  |

Values match `results/descriptors.csv` (rows with variant = `A`).
