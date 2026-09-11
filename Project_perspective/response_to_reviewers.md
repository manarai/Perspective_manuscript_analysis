# Response to reviewers — Perspective route

**Submission ID** c8b832ae-cd82-4dd8-bcbe-3f9e53192a94
**Manuscript** Operator Descriptors: A Comparative Language for Biological Dynamics *(retitled; see "On the manuscript title")*
**Article type (revised)** Perspective (reclassification granted 25 Aug 2026)

We thank Dr Jolly and Professor Albert for reclassifying the manuscript, and both reviewers for careful reading. Reviewer 1's rejection rested on a mismatch between article type and content, which reclassification resolves. Reviewer 2's five major criticisms bear on framework clarity and falsifiability; each has driven a specific change.

**A note on where things live.** The Perspective's main text stays within the 3,000-word limit and carries the descriptor framework, the mechanism-independence claim, and its limits. Computational detail — model equations and reduction assumptions (S1), the full metric-sensitivity sweep (S2), and the concordance protocol (S3) — sits in supplementary material. Where this letter cites a number, the main text reports it and the supplementary material shows how it was obtained.

**Four things we got wrong or overstated**, all corrected in this revision and stated plainly rather than absorbed into the text:

1. **ν does not by itself predict transient amplification.** Non-normality (ν > 0) is a necessary precondition — for stable normal A, ‖e^{At}‖ decreases monotonically from ‖I‖ = 1, so G_max ≡ sup_{t≥0} ‖e^{At}‖ = 1 — but not sufficient. Goodwin (ν = 0.14, ω = −0.34, G_max = 1.00) is a counterexample in our own panel. Amplification is therefore located by the functional descriptor ω, with ν as a structural precondition.

2. **The descriptors are metric-relative, not invariant.** This is the overstatement Reviewer 2's question about the numerical abscissa uncovered. ν, ω and G_max are anchored to a stated metric; Q alone is exactly invariant. S2 tabulates ranges under three per-state rescalings for every computation. Under our stated metric the six amplifier assignments (TCR, Neubert–Caswell, Murphy–Miller, Selkov, NF-κB, heterogeneous repressilator) hold across all three rescalings; Selkov's ω remains positive throughout (min +0.294, baseline +0.294, max +0.937). Two marginal non-amplifiers (symmetric repressilator, Goodwin) flip sign, and NF-κB drops to ω = 0.001 without flipping. We report all of it.

3. **Some ν = 0 results were symmetry artifacts.** The symmetric three-variable repressilator's Jacobian is circulant and therefore exactly normal, a property of that reduction rather than of oscillators. We added heterogeneous variants of the repressilator, toggle and relay to probe this, and report the mixed outcome: the repressilator's ν = 0 is an artifact (heterogeneous variant, ν = 0.27), the toggle's likewise (0.000 → 0.192) although its contractivity assignment survives unchanged (ω = −0.74, G_max = 1.00), and the relay's ν = 0 persists under per-node decay heterogeneity because a diagonal shift preserves symmetry. The last is a property of modelling the relay as reversible diffusion–decay; real relays are typically directional, and a directional relay is non-normal, as the TCR cascade in the panel shows.

4. **The Fisher-memory anchor for "information processing" is withdrawn.** The Ganguli–Huh–Sompolinsky (2008) Fisher memory curve J(t) = vᵀ(e^{At})ᵀC⁻¹e^{At}v uses the stationary covariance C solving AC + CAᵀ + Σ = 0, not a metric one chooses. Setting C = I forces Σ = −(A + Aᵀ), which requires ω ≤ 0: exactly the non-amplifying regime. We verified this numerically (for Murphy–Miller, ω = 3.71, the required Σ has a −7.4 eigenvalue and is not a valid noise model). J is coordinate-invariant and is the natural next step, but computing it requires a noise structure per system that this illustration does not supply. We cite Ganguli et al. as the coordinate-invariant target, treat it as future work, and have removed the claim, the figure annotation carrying it, and the corresponding phrase from the title.

---

## Reviewer 1

> "I recommend rejection... the authors themselves state four times that the work is a conceptual perspective containing no data or analysis..."

We agree with the diagnosis. The manuscript was prepared and submitted as a Perspective, and the reviewer's reading is an accurate description of a Perspective assessed against research-article criteria. The editors reclassified it on 25 August 2026. The self-descriptions that prompted the objection are removed; they functioned as apologies for the absence of analysis rather than as statements of scope.

Although a Perspective does not require original analysis, we have added an illustrative computation, because Reviewer 2's objections could not be answered by rewording.

---

## Reviewer 2 — Major issues

### M1. Classes vague; "non-normal" and "spectral" too broad to correspond to specific information processing

We agree with all three parts of this criticism, and it drove the largest structural change in the revision.

**"Spectral" is dropped entirely**, replaced by the scalar Q = |Im λ₁| / |Re λ₁|: the phase advance of the leading mode per unit decay time.

**Structural and functional descriptors are separated.** The reviewer's point that non-normal structure "could correspond to any information processing" is right, and the fix is to make non-normality a precondition rather than a signature:

| Kind | Descriptor | Meaning |
|---|---|---|
| Structural | ν (Henrici) | departure from normality; necessary for transient amplification |
| Structural | r_PR / n | effective rank of the singular spectrum |
| Functional | Q | phase advance of leading mode per decay time |
| Functional | ω = λ_max((A + Aᵀ)/2) | initial growth rate of ‖e^{At}δ‖ under stated metric; ω > 0 ⟺ G_max > 1 for stable A |
| Functional | G_max | peak of ‖e^{At}‖₂ |

**On thresholds.** ω = 0 and Q = 1 are different kinds of boundary and the manuscript now says so.

- ω = 0 is the exact boundary of a theorem: ω > 0 ⟺ G_max > 1 for stable A, with monotone norm decrease on one side and transient growth on the other. What it lacks is an effect-size floor — ω = 0.009 and ω = 3.71 score alike under a sign test. We handle this by reporting the paired G_max wherever the sign appears, so the reader can see whether amplification is functional (Murphy–Miller, G_max = 2.99) or negligible (predator–prey, 1.007), rather than by imposing an arbitrary numerical floor.
- Q = 1 is a convention, and we now give its rationale: Q is radians of phase advance per e-folding, so Q > 1 means the leading mode completes at least 1/2π ≈ 16% of a cycle before decaying by a factor e. Our oscillating systems sit at Q = 2.1 to 4.4, that is 0.33 to 0.70 cycles per e-fold.

Region names are informal groupings of coordinate positions. A system's descriptor values, not its label, carry the claim, and the values appear wherever the label does.

**Why ν is the wrong functional proxy.** TCR and Neubert–Caswell predator–prey have nearly identical ν (0.672 vs 0.655) but excess gain (G_max − 1) differing more than forty-fold (0.302 vs 0.007), while Goodwin has ν = 0.14 with G_max = 1.00 exactly. Grouping by ν alone would separate functionally comparable systems and unite functionally different ones. This is now a stated finding in the main text rather than a defect the illustration conceals.

**On "what kind of information processing".** The coordinate-invariant answer is the Fisher memory curve (Ganguli et al. 2008), which we cite as the target and defer, for the reason given in item 4 above. We do not claim that ‖e^{At}δ‖² *is* Fisher information under any assumption on our metric.

**On contractive versus decision.** Contractivity is dynamical (G_max = 1, monotone norm decrease); low-dimensional decision is structural (a slow subspace dominates). At its stable attractor the symmetric toggle is exactly normal with r_PR/n = 0.996 — the decision signature lives at the saddle between attractors, not at either attractor. Our illustration linearises only at stable operating points, so it contains no positive example of the low-effective-rank axis, and the toggle and relay are consequently not separable by (Q, ν, ω). We report this as a degeneracy that identifies the missing axis, and name saddle linearisation as the work that would close it.

### M2. Figure confusing

We agree, and we thank the reviewer for enumerating the problems specifically enough to act on.

Fig 1 was rebuilt from scratch: (a) three unrelated systems with distinct microscopic dynamics; (b) local operator estimation, with the metric choice shown as an explicit step; (c) the (Q, ω) plane with the four descriptor regions labelled. No metaphorical illustrations; all terminology matches the text term by term. The specific complaints — mountain pass, "partial state inform", horizontal arrows, contractive dynamics drawn with diverging arrows, "Reduced coarse-grain observables" placement, figure–text mismatch in the region names — refer to elements no longer present. The diverging-arrow contractive panel was simply an error.

Fig 2 carries the illustration: (a) the model computations placed in (Q, ω) with sensitivity range bars; (b) each model paired with its spectrum-matched surrogate, showing the ν drop to zero by construction and the resulting change in G_max; (c) ‖e^{At}δ‖ against t for one representative per region, annotated as transient gain only — the earlier Fisher-memory annotation was removed with the withdrawn claim.

### M3. Only one example per archetype; need more to show mechanism-independence

The reviewer is right that one example per category cannot demonstrate convergence. Twelve computations across nine source systems in eight domains (Fig 2a; equations, dimensions and reduction assumptions in S1):

| Region | Model | ω | Q | ν | G_max | Domain | Mechanism |
|---|---|---|---|---|---|---|---|
| Amplifying, non-oscillatory | TCR proofreading | +0.075 | 0.0 | 0.672 | 1.30 | immune signalling | kinetic proofreading cascade |
| | Neubert–Caswell predator–prey | +0.012 | 0.0 | 0.655 | 1.007 | ecology | consumer–resource |
| Amplifying and oscillating | Murphy–Miller balanced E–I | +3.71 | 2.9 | 0.951 | 2.99 | neural | balanced excitation–inhibition |
| | Selkov glycolysis | +0.29 | 3.5 | 0.80 | 1.62 | metabolism | autocatalytic phosphofructokinase |
| | NF-κB / IκB | +0.009 | 3.9 | 0.466 | 1.07 | transcription | negative feedback with delay |
| | Repressilator (heterogeneous) | +0.095 | 4.4 | 0.270 | 1.04 | synthetic | negative feedback with delay |
| Oscillating, not amplifying | Repressilator (symmetric) | −0.30 | 4.0 | 0.000 | 1.00 | synthetic | negative feedback with delay |
| | Goodwin | −0.34 | 2.1 | 0.140 | 1.00 | transcription | negative feedback with delay |
| Not amplifying, not oscillating | Toggle (symmetric / heterogeneous) | −0.88 / −0.74 | 0.0 | 0.000 / 0.192 | 1.00 | hematopoiesis | mutual repression |
| | Relay (symmetric / heterogeneous) | −0.20 / −0.19 | 0.0 | 0.000 | 1.00 | signalling | diffusion–decay chain |

Three points, which are what the reviewer asked us to demonstrate.

- **The amplifying region holds four mechanistically unrelated systems**: kinetic proofreading (immune), balanced excitation–inhibition (neural), autocatalytic glycolysis (metabolism), consumer–resource dynamics (ecology). No shared components, topology, or characteristic timescale.
- **The oscillating region holds a metabolic oscillator alongside three transcriptional ones.** We added Selkov specifically because NF-κB, the repressilator and Goodwin are one mechanism — delayed negative feedback of a regulator on its own production — in three parameterisations, and a claim of mechanism-independence resting on those three alone would rest on a family resemblance. Selkov's autocatalytic product activation of phosphofructokinase is biochemically distinct.
- **Within a region, effect sizes vary widely.** G_max across the amplifying region runs from 1.007 to 2.99. We do not claim tight co-location; "occupy the same qualitative region despite unrelated mechanism" is the defensible statement and the one the manuscript makes.

**The toggle / relay degeneracy.** All four computations land at (Q = 0, ω < 0) and are not separable by (Q, ν, ω): a bistable switch and a passive relay give the same signature at their stable attractors. This is a real degeneracy, and it points at the axis the illustration does not exercise — effective rank at the saddle between toggle attractors rather than at either attractor. Addressing it requires linearising at the saddle, which we name as future work rather than claim the current axes distinguish the two.

### M4. "From Data to Archetypes" abstract; "dominant mode replicates across related datasets" unclear

The section was rewritten in full and renamed **"From Data to Descriptors"**, the "Archetypes" going with the term itself (see m3). Every remaining step now has an explicit procedure or is identified as a deliberate open choice.

"Replicates" now means the maximum principal angle between leading Koopman subspaces across bootstrap resamples, with rank fixed *a priori* as the largest r with σ_r/σ₁ > 0.1, and concordance at an angle below 30°. The threshold's margin is established rather than asserted: S3 reports a permutation null on shuffled snapshots, against which random subspaces sit near orthogonality, so 30° is not an easy bar.

All Fig 2 computations use scipy on the analytic Jacobians in S1; no claim in this manuscript depends on the scjdo preprint's data-driven pipeline, which is cited as the implementation of the follow-up work and as where a reader can see the concordance protocol running on real data.

### M5. "Testable Predictions" hedged so heavily as to be unfalsifiable

We agree, and we took the reviewer's diagnosis rather than his suggested remedy: the section is kept but reduced to one prediction, stated as a prediction for future testing rather than as a test scored on this panel. The consistency conditions move out from under the Predictions heading.

**Prediction P1** (not scored here) — *Descriptor values will predict functional labels assigned independently in source papers better than they predict molecular domain, on a panel large enough that domain grouping is not defeated by singletons.*

Of the twelve computations, the source-paper label agrees with descriptor placement for eight: three carry an amplifier label and have ω > 0 (TCR, Murphy–Miller, Neubert–Caswell); five carry an oscillator label and have Q > 1 (NF-κB, Goodwin, Selkov, both repressilator regimes). The groups do not overlap. Murphy–Miller and NF-κB have both descriptors positive and their source papers describe both behaviours; Selkov also has both positive while its source paper labels only the oscillation, a case where the descriptors say more than the label rather than disagreeing with it. The remaining four computations are the toggle and relay variants, indistinguishable from each other rather than misplaced.

Two limits on what this shows, both stated in the manuscript. These models were chosen as canonical examples of their labels, not held out, so the agreement is illustration rather than test. And at nine source systems across eight domains, domain grouping is almost all singletons, so the descriptor-versus-domain comparison the reviewer asked about cannot be run informatively here; it requires a curated panel with multiple representatives per domain, which is the target of the data-driven follow-up. S1 records the verbatim functional label from each source paper and where it appears, so the assignment is auditable rather than taken on trust.

**Consistency conditions.** Two identities hold analytically for stable A: G_max > 1 implies ν > 0, and sign(ω) = sign(G_max − 1). All twelve computations satisfy both, which verifies our code rather than the framework; we report them in one sentence as such.

---

## Reviewer 2 — Minor issues

**m1. What defines the "operator level"?** Defined operationally at first use: descriptors are functions of the Jacobian at the stated operating point (analytic case) or of the estimated Koopman generator over a stated observable set (data-driven case). The reviewer is right that an infinite-dimensional operator can represent anything; that is why the specification of the rank-r reduction and the observable set is part of each descriptor's definition rather than an implementation detail.

**m2. Formatting of the Koopman equation.** Reformatted in standard notation.

**m3. Is an "operator class" the same as an "archetype"?** We dropped both terms, since each implied a categorical assignment the framework does not support. The manuscript now uses "descriptor region" for a subset of (Q, ω) space and describes systems by coordinate values. Regions are not mutually exclusive. The two boundaries are of different kinds and both are stated as such: ω = 0 tracks a theorem, Q = 1 is a convention with a phase-per-decay-time rationale (M1). The title changed accordingly.

**m4. What does it mean for perturbations to "rotate"?** Phrase removed, replaced by Q or by direct reference to the leading complex-conjugate eigenpair.

**m5. In what sense is the numerical abscissa an invariant?** It is not, and the question is what uncovered the error. ω = λ_max((A + Aᵀ)/2) is the initial growth rate of ‖e^{At}δ‖ for worst-case unit δ under the stated metric; it is metric-relative and satisfies ω > 0 ⟺ G_max > 1 for stable operators. S2 tabulates its sensitivity, with ν and G_max, across three rescalings for every computation.

**m6. Pseudospectrum undefined.** Now defined at first use: σ_ε(A) is the set of complex λ with ‖(λI − A)⁻¹‖ > 1/ε, quantifying the sensitivity of the eigenvalue set to perturbation. ν is a scalar summary of the departure from normality that pseudospectral analysis characterises in the complex plane; Trefethen & Embree (2005) is cited for the full picture.

---

## On the manuscript title

From "Dynamical Operator Archetypes in Biological Information Flow" to **"Operator Descriptors: A Comparative Language for Biological Dynamics."**

"Archetypes" goes with the term itself, following the reviewer's concern that categorical language obscured the claim. "Information Flow" goes with the withdrawn Fisher-memory anchor: without it the manuscript makes no information-theoretic claim, and a title should not carry one. The new title names what the Perspective is for rather than committing to assignments the illustration does not always support.

---

## Summary of changes

- Reclassified and reformatted as a Perspective (≤3,000-word main text, ≤70-word abstract, ≤70 references); computational detail in S1–S3.
- Retitled; "archetype" and "class" replaced throughout by "descriptor region", including the renamed section "From Data to Descriptors".
- Structural (ν, r_PR/n) and functional (Q, ω, G_max) descriptors separated; ν stated as a necessary precondition for amplification, not a signature of it.
- ω = 0 distinguished as a theorem boundary, with G_max reported alongside every sign as the effect-size check; Q = 1 given an explicit phase-per-decay-time rationale.
- Invariance overstatement withdrawn; metric-relativity stated openly with per-descriptor ranges across three rescalings in S2, including the two assignments that flip sign.
- Fisher-metric claim withdrawn, with the figure annotation and title phrase that carried it; Ganguli et al. 2008 cited as the coordinate-invariant next step.
- Fig 1 rebuilt from scratch; Fig 2 added with twelve computations across nine source systems in eight domains and spectrum-matched surrogates.
- Selkov glycolytic oscillator added so the oscillating region carries mechanism diversity rather than three variants of delayed transcriptional feedback.
- Heterogeneous variants of repressilator, toggle and relay added to probe whether ν = 0 assignments are symmetry artifacts; mixed result reported.
- One prediction stated for future testing, not scored on the panel it was built from; source-paper labels recorded verbatim in S1; consistency conditions compressed to a code check.
- All Fig 2 computations run with scipy on analytic Jacobians; notebook, Jacobians and sensitivity sweep deposited with a Zenodo DOI (minted at resubmission).
- Self-descriptions as a conceptual perspective removed from the body.

We are grateful for the reviewers' engagement, and to Reviewer 2 in particular: one of these questions led us to an error we would otherwise have published.
