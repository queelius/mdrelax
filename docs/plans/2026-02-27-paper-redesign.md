# Paper Redesign: Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption

**Date:** 2026-02-27
**Status:** Approved design

## Central Claim

The standard C1-C2-C3 model for masked series system data assumes non-informative masking (C2). In practice, the masking mechanism is unknown and potentially informative. We characterize the bias that arises when C2 is violated but assumed to hold, and map the breakdown boundary through simulation.

## Audience

Statistical methodologists working on masked data / competing risks.

## Scope Decisions

- **C2 only.** C3 violations barely affect estimates (existing simulation evidence). C1 violations are a different paper.
- **Misspecification theorem only.** No FIM comparisons, no KL-divergence models.
- **One simulation study, sweep over violation severity.** Both exponential and Weibull.
- **P(C|K) as perturbation tool, not proposed model.** We don't claim the P matrix is the right masking model. It's a device for generating controlled C2 violations.

## What's Cut from Current Paper

- All C3 relaxation (power-weighted masking, alpha parameter, Definitions 3.7-3.9)
- All C1 relaxation (exp_series_relaxed_c1)
- KL-divergence constrained models (Definition 3.6)
- Rank-based masking model (Definition 3.3)
- Fisher information comparisons and derivations
- All current simulation Studies 1-5
- The `identifiability.tex` section (8 theorems)
- The `discussion.tex` as written

## What's Kept/Adapted

- Background: series systems, masked data, C1-C2-C3 likelihood (exponential + Weibull)
- General likelihood under C1 alone (existing Theorem 3.1)
- P(C|K) Bernoulli model (existing Definition 3.4, Theorem 3.5) — reframed as perturbation tool
- Existing simulation infrastructure (scenarios 3, 4b for exponential; W3, W4 for Weibull)
- The non-identifiability result for joint (theta, P) estimation (existing §3.6.4)

## Section Design

### 1. Introduction (~1.5 pages)

The masked data problem and C1-C2-C3 framework. C2 is the strongest, least justifiable assumption: the masking mechanism is unknown, complex, potentially non-stationary. Prior work has built alternative masking models (Lin & Guess 1994, Guttman et al. 1995, Mukhopadhyay & Basu 2006, Craiu & Reiser 2010). We ask the complementary question: how sensitive is standard C1-C2-C3 inference to C2 violations?

Contributions:
1. Misspecification theorem characterizing the pseudo-true parameter under C2 violation.
2. Systematic simulation sweep quantifying bias/RMSE as a function of violation severity for exponential and Weibull systems.
3. The practical finding: system-level reliability is robust to C2 violations; individual component rates are not.

### 2. Background (~2 pages)

- Series system model: hazard, survival, masked data setup.
- The C1-C2-C3 likelihood for exponential and Weibull components.
- Brief prior work review: who has relaxed C2 and how. Frame the gap: no systematic sensitivity analysis exists.

### 3. Sensitivity Framework (~3 pages)

**3.1 General likelihood under C1 alone.**
Keep existing Theorem 3.1. This shows the masking probability pi_{k,c}(t) weights the hazard contributions.

**3.2 The P(C|K) Bernoulli model as perturbation tool.**
Keep existing Definition 3.4 (P matrix) and Theorem 3.5 (likelihood). Reframe: we don't claim this is how masking works. It's a device that lets us generate data with known, controlled C2 violations.

**3.3 Measuring violation severity.**
Define a scalar measure of C2 violation from the P matrix. Candidates:
- Max column asymmetry: max_{j,k,k'} |P[j,k] - P[j,k']|
- Frobenius distance from nearest C2-satisfying P matrix
- Something simpler: the range of off-diagonal entries in a single column
This becomes the x-axis of the simulation sweep.

**3.4 Misspecification theorem.**
State and prove: when C2 is violated (P asymmetric) but you fit C1-C2-C3, the MLE converges to a pseudo-true parameter theta-dagger that satisfies E[score_wrong | theta-dagger] = 0. Show that theta-dagger differs from theta-star unless P satisfies C2. Show that sum(theta) remains approximately correct (system hazard robust). The bias in individual components depends on the asymmetry of P.

**3.5 The identifiability trap.**
Existing result: joint estimation of (theta, P) fails — theta and P are confounded. Total hazard identifiable, individual rates not. This means you can't "just estimate P from the data." Implication: sensitivity analysis is the right tool.

### 4. Simulation Study (~4 pages)

**4.1 Design.**

Exponential setup: m=3, theta=(1, 1.5, 2), tau=3, n=200, B=200 replications.
Weibull setup: m=2, shapes=(2, 1.5), scales=(3, 4), tau=8, n=200, B=100 replications.

P matrix sweep: parameterize a family of P matrices indexed by a scalar "severity" parameter s in [0, 1]:
- s=0: P satisfies C2 (all columns identical, off-diagonal = 0.5)
- s=1: maximally asymmetric (e.g., P[1,2]=0.9, P[2,1]=0.1)
- Intermediate values interpolate linearly

For each severity level: generate data from P(s), fit C1-C2-C3 model, compute bias, RMSE, coverage for each parameter.

**4.2 Results.**

Primary figures:
- Bias curves: bias(theta_j) vs severity s. One panel per component. Show individual components diverge, total hazard stays flat.
- RMSE curves: RMSE(theta_j) vs severity s.
- Coverage curves: 95% CI coverage vs severity s. Show when coverage breaks down.

Both exponential and Weibull on same figures (different line styles or panels).

**4.3 Interpretation.**

Key findings (expected from existing evidence):
- System-level quantities (total hazard, system MTTF) are robust to C2 violations.
- Individual component parameters show bias proportional to the asymmetry of P.
- For moderate violations (s < ~0.3?), bias is within one SE — C1-C2-C3 is "good enough."
- Weibull shape-scale interaction may amplify or dampen sensitivity relative to exponential.

### 5. Discussion (~1.5 pages)

**5.1 When does C2 matter?**
Individual component reliability estimation: yes. System-level reliability: mostly no.

**5.2 Practical guidance.**
Since the masking mechanism is unknowable, run sensitivity analysis: fit C1-C2-C3 under several plausible P matrices. If conclusions are stable, C2 doesn't matter for your application.

**5.3 Limitations.**
- Bernoulli P(C|K) model is one perturbation family; other violation structures could behave differently.
- Time-dependent masking P(C|K,T) not studied.
- Limited to exponential and Weibull; other lifetime distributions may behave differently.

### 6. Conclusion (~0.5 pages)

### Appendix
- Proof of misspecification theorem.
- Full simulation tables.
- Software note: mdrelax R package.

## Implementation Notes

### What needs to be written from scratch
- New introduction reframed around sensitivity
- §3.3 violation severity measure
- §3.4 misspecification theorem (adapt from existing, simplify for C2 only)
- §4 simulation study (new sweep design; can reuse simulation infrastructure)
- New discussion
- New conclusion

### What can be adapted from existing paper
- §2 background (trim the current background.tex)
- §3.1-3.2 general likelihood + P matrix (trim relaxed_models.tex, cut C3/C1 material)
- §3.5 identifiability trap (extract from existing relaxed_models.tex §3.6.4)
- Appendix proofs

### LaTeX files to modify
- `main.tex` — new title, new abstract, updated section inputs
- `sections/introduction.tex` — rewrite
- `sections/background.tex` — trim, add prior work review
- `sections/relaxed_models.tex` — gut and rewrite as sensitivity_framework.tex
- `sections/simulations.tex` — replace entirely with new sweep study
- `sections/discussion.tex` — rewrite
- `sections/conclusion.tex` — rewrite
- `sections/identifiability.tex` — delete (theorem moves to §3)
- `sections/appendix.tex` — update
- `refs.bib` — add Lin & Guess, Guttman et al., Mukhopadhyay & Basu, Craiu & Reiser

### R code needed
- New simulation script: `paper/run_sensitivity_sweep.R`
  - Uses existing `rexp_series_md_c1_c3()` and `mle_exp_series()`
  - Uses existing `rwei_series_md_c1_c3()` and `mle_wei_series()`
  - Sweeps P matrix severity parameter
  - Generates bias/RMSE/coverage curves

### Estimated page count
12-15 pages including references and appendix.
