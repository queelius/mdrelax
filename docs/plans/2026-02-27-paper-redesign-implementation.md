# Paper Redesign Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Rewrite the mdrelax paper as a focused sensitivity analysis: "How sensitive is C1-C2-C3 inference to violations of the non-informative masking assumption (C2)?"

**Architecture:** Replace the current 7-section, multi-model paper with a 6-section, single-question paper. The P(C|K) Bernoulli model is reframed as a perturbation tool (not a proposed model). One simulation study sweeps C2 violation severity for exponential + Weibull. One misspecification theorem. No C3, no C1, no FIM, no KL.

**Tech Stack:** LaTeX (pdflatex + bibtex), R (existing mdrelax package functions)

**Design doc:** `docs/plans/2026-02-27-paper-redesign.md`

---

### Task 1: Create the simulation sweep script

This is the R script that generates all empirical results for the paper. Do this first so figures/tables are available when writing LaTeX.

**Files:**
- Create: `paper/run_sensitivity_sweep.R`

**Step 1: Write the sweep script**

The script should:

1. Load the package with `devtools::load_all()`.
2. Define a function `make_sweep_P(m, s, base_p)` that creates a P matrix parameterized by severity `s in [0,1]`:
   - `s=0`: uniform off-diagonal = `base_p` (C2 satisfied)
   - `s=1`: maximally asymmetric. For m=3, one design: column k gets off-diagonal entries `base_p + s * spread_k` where spread varies per column, creating asymmetry. A clean approach: for column k (true failed component = k), component j's inclusion probability is `base_p + s * (offset_jk)` where offsets sum to zero per column and are clamped to [0.05, 0.95]. Example for m=3: create a fixed "direction matrix" D with zero diagonal, column sums = 0, and interpolate `P = P_uniform + s * D`.
3. For the exponential sweep: `m=3, theta=(1, 1.5, 2), tau=3, n=200, B=200`. Sweep `s` over `seq(0, 1, by=0.1)` (11 levels). At each level: generate data with `rexp_series_md_c1_c3(n, theta, P(s), tau)`, fit with `mle_exp_series(t, C, delta)`, record `theta_hat`, `se`, `converged`.
4. For the Weibull sweep: `m=2, shapes=(2, 1.5), scales=(3, 4), tau=8, n=200, B=100`. Sweep `s` similarly. Generate with `rwei_series_md_c1_c3(...)`, fit with `mle_wei_series(...)`.
5. Compute per-severity-level: mean bias, RMSE, 95% CI coverage (using `theta_hat +/- 1.96*se`), total hazard bias.
6. Save results to `paper/data/sensitivity_sweep_exp.rds` and `paper/data/sensitivity_sweep_wei.rds`.
7. Generate 3 PDF figures to `inst/simulations/figures/`:
   - `fig_bias_vs_severity.pdf` — bias for each parameter vs s, plus total hazard. Two panels: exponential (left), Weibull (right).
   - `fig_rmse_vs_severity.pdf` — same layout for RMSE.
   - `fig_coverage_vs_severity.pdf` — same layout for coverage.

For the "direction matrix" D (step 2), a good choice for m=3:
```
D = | 0    0.3  -0.3 |
    | -0.2  0    0.2 |
    | 0.2  -0.3  0.1 |
```
This means: when component 1 fails (col 1), component 3 is MORE likely included (+0.2) and component 2 is LESS likely (-0.2). Different columns create different asymmetries. Scale by `s` and clamp.

For m=2 (Weibull), a simpler direction:
```
D = | 0    0.4 |
    | -0.4  0  |
```
When component 1 fails, component 2 is LESS likely included. When component 2 fails, component 1 is MORE likely. Scale by `s`.

**Step 2: Run the script in quick mode (B=20) to verify it works**

Run: `Rscript paper/run_sensitivity_sweep.R --quick` (or add a `quick_mode` flag at the top of the script).

Expected: Script completes without error, produces .rds files and 3 PDF figures.

**Step 3: Run full simulations**

Run: `Rscript paper/run_sensitivity_sweep.R`

This will take ~15-30 minutes. Run in background.

**Step 4: Commit**

```bash
git add paper/run_sensitivity_sweep.R inst/simulations/figures/fig_bias_vs_severity.* inst/simulations/figures/fig_rmse_vs_severity.* inst/simulations/figures/fig_coverage_vs_severity.*
git commit -m "Add sensitivity sweep simulation for paper redesign"
```

---

### Task 2: Rewrite main.tex — new title, abstract, section structure

**Files:**
- Modify: `paper/main.tex`

**Step 1: Update title, abstract, and section inputs**

New title: "Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption"

New abstract (~150 words): The standard likelihood for masked series system data rests on three conditions (C1-C2-C3). Among these, C2 — non-informative masking — is the strongest and least verifiable: it asserts that the candidate set formation reveals nothing about which component failed. We study the sensitivity of maximum likelihood estimation to violations of C2. Using a Bernoulli perturbation model P(C|K) to generate controlled violations, we derive the pseudo-true parameter under C2 misspecification and show that the total system hazard remains well-identified while individual component parameters are biased. A simulation sweep over violation severity for both exponential and Weibull series systems maps the breakdown boundary. For moderate violations, the standard model remains adequate for system-level inference; for individual component estimation, sensitivity analysis over plausible masking structures is recommended.

Section inputs become:
```latex
\input{sections/introduction}
\input{sections/background}
\input{sections/sensitivity_framework}
\input{sections/simulations}
\input{sections/discussion}
\input{sections/conclusion}
```

Remove `\input{sections/identifiability}` and `\input{sections/relaxed_models}`.

**Step 2: Verify compilation**

Run: `cd paper && pdflatex main.tex`
Expected: Compiles (may have warnings about missing refs — that's fine at this stage).

**Step 3: Commit**

```bash
git add paper/main.tex
git commit -m "Rewrite main.tex: new title, abstract, section structure"
```

---

### Task 3: Rewrite introduction.tex

**Files:**
- Modify: `paper/sections/introduction.tex`

**Step 1: Rewrite the introduction**

Structure:
- Para 1: The masked data problem in series systems. Cite Usher 1988, Lin 1993.
- Para 2: C1-C2-C3 conditions. C1 (failed component in candidate set) is natural. C3 (parameter-independent masking) is often reasonable. C2 (non-informative masking) is the strongest assumption: it requires that knowing which component failed provides no information about the candidate set beyond C1.
- Para 3: Why C2 is suspect. The masking mechanism is a diagnostic process — unknown, complex, potentially non-stationary. Even well-intentioned diagnostic procedures may preferentially include certain components. Prior work (Lin & Guess 1994, Guttman et al. 1995, Mukhopadhyay & Basu 2006, Craiu & Reiser 2010) has proposed alternative masking models. But no one has systematically studied: what happens to the standard MLE when C2 is wrong?
- Para 4: Our contributions (3 items from design doc).
- Para 5: Outline of the paper.

Keep it to ~1.5 pages.

**Step 2: Commit**

```bash
git add paper/sections/introduction.tex
git commit -m "Rewrite introduction for sensitivity analysis framing"
```

---

### Task 4: Trim background.tex

**Files:**
- Modify: `paper/sections/background.tex`

**Step 1: Read current background.tex**

Read the full file first. Keep: series system model, hazard/survival, masked data setup, C1-C2-C3 likelihood for exponential, Weibull likelihood. Cut: anything about relaxed models (that moves to §3). Add: brief prior work subsection on dependent masking.

**Step 2: Add prior work subsection**

After the C1-C2-C3 likelihood section, add a subsection "Prior Work on Dependent Masking" (~0.5 pages). Cite:
- Lin & Guess (1994): proportional masking probabilities, 2-component, closed-form MLE.
- Guttman et al. (1995): Bayesian inference for 2-component dependent masking.
- Mukhopadhyay & Basu (2006): EM-based MLE for general m-component systems, masking depends on cause.
- Craiu & Reiser (2010): conditional masking probability models, identifiability.
- Frame the gap: all prior work proposes alternative masking models. None studies what happens when you don't use them — i.e., the sensitivity of the standard model to masking violations.

**Step 3: Commit**

```bash
git add paper/sections/background.tex
git commit -m "Trim background, add prior work on dependent masking"
```

---

### Task 5: Create sensitivity_framework.tex (§3)

**Files:**
- Create: `paper/sections/sensitivity_framework.tex`
- Delete: `paper/sections/relaxed_models.tex` (content absorbed)
- Delete: `paper/sections/identifiability.tex` (theorem moves here)

**Step 1: Write §3.1 — General Likelihood Under C1 Alone**

Adapt existing Theorem 3.1 and Remark 3.1 from `relaxed_models.tex`. Keep the proof. This is the anchor: under C1 alone, the likelihood has pi_{k,c}(t) weighting each hazard contribution. Under C2+C3, these weights factor out.

**Step 2: Write §3.2 — The Bernoulli Perturbation Model**

Adapt existing Definition 3.4 (P matrix), equation for pi_k(c), and Theorem 3.5 (likelihood). Reframe: "We do not claim this is how masking works in practice. Rather, the P matrix provides a device for generating data with known, controlled C2 violations, parameterized by the departure of P from its C2-satisfying form."

Keep the Remark about C2 in terms of P (constant columns = C2 holds).

**Step 3: Write §3.3 — Measuring Violation Severity**

New content. Define the severity measure used in simulations. Something like:

```latex
\begin{definition}[C2 Violation Severity]
For a P matrix, define the severity of C2 violation as
s(P) = max_{j != k, k' != k} |P[j,k] - P[j,k']| / max possible.
When s=0, all columns identical (C2 holds). When s=1, maximally asymmetric.
\end{definition}
```

Or more concretely, describe the linear interpolation family used in simulations: `P(s) = P_uniform + s * D` where D is a fixed direction matrix.

**Step 4: Write §3.4 — Misspecification Theorem**

Adapt from the existing C2 misspecification theorem (Theorem 4.7 in identifiability.tex). Simplify to focus on C2 only. The key result:

When C2 is violated but we fit C1-C2-C3, the misspecified score equates `sum_{k in c} h_k(t; theta_k)` to the weighted sum `sum_{k in c} h_k(t; theta_k) * pi_k(c)`. The pseudo-true parameter theta-dagger satisfies E[score_wrong(theta-dagger)] = 0 under the true P(C|K) model.

For exponential series: show that sum(theta) is approximately preserved (the misspecified score for the total hazard has the right expectation when masking is symmetric *on average*), but individual theta_j are biased toward the "apparent" rates that absorb the masking asymmetry.

Put proof in appendix.

**Step 5: Write §3.5 — The Identifiability Trap**

Adapt existing result from relaxed_models.tex §3.6.4. Joint estimation of (theta, P) fails: theta and P are confounded. Total hazard identifiable, individual components not. This means sensitivity analysis (not model estimation) is the right response.

**Step 6: Commit**

```bash
git add paper/sections/sensitivity_framework.tex
git rm paper/sections/relaxed_models.tex paper/sections/identifiability.tex
git commit -m "Create sensitivity framework section, remove old model/identifiability sections"
```

---

### Task 6: Rewrite simulations.tex (§4)

**Files:**
- Modify: `paper/sections/simulations.tex` (full rewrite)

**Step 1: Write §4.1 — Design**

Describe:
- The sweep: severity parameter s from 0 to 1 in steps of 0.1.
- Exponential setup: m=3, theta=(1, 1.5, 2), tau=3, n=200, B=200.
- Weibull setup: m=2, shapes=(2, 1.5), scales=(3, 4), tau=8, n=200, B=100.
- The P(s) family: `P(s) = P_uniform + s * D`.
- At each s: generate data from relaxed C2 with P(s), fit with standard C1-C2-C3 model.
- Metrics: bias, RMSE, 95% CI coverage. Also track total hazard sum(theta).

**Step 2: Write §4.2 — Results**

Include figures (from Task 1):
- Figure 1: Bias vs severity (exponential + Weibull panels)
- Figure 2: RMSE vs severity
- Figure 3: Coverage vs severity

Include a summary table with selected severity levels (s=0, 0.3, 0.5, 0.7, 1.0).

Populate table values from actual simulation results (read from .rds files).

**Step 3: Write §4.3 — Interpretation**

Key findings (verify against actual results):
- Total hazard: bias stays near zero across all severity levels.
- Individual components: bias grows approximately linearly with severity.
- Coverage: degrades smoothly; find the approximate threshold where coverage drops below 90%.
- Exponential vs Weibull comparison: note whether Weibull is more or less sensitive.

**Step 4: Commit**

```bash
git add paper/sections/simulations.tex
git commit -m "Rewrite simulations section with sensitivity sweep results"
```

---

### Task 7: Rewrite discussion.tex and conclusion.tex

**Files:**
- Modify: `paper/sections/discussion.tex`
- Modify: `paper/sections/conclusion.tex`

**Step 1: Rewrite discussion**

Three subsections:
- §5.1 When Does C2 Matter? — Individual component rates: yes. System hazard / MTTF: mostly no.
- §5.2 Practical Guidance — Run sensitivity analysis: fit C1-C2-C3 under several plausible P matrices. If substantive conclusions don't change, C2 violations are not a concern for your problem.
- §5.3 Limitations — Bernoulli P(C|K) is one perturbation family. Time-dependent masking not studied. Exponential + Weibull only.

~1.5 pages total.

**Step 2: Rewrite conclusion**

Short (~0.5 pages). Three contributions: (1) misspecification theorem, (2) simulation sweep, (3) practical finding about total hazard robustness. Mention mdrelax R package.

**Step 3: Commit**

```bash
git add paper/sections/discussion.tex paper/sections/conclusion.tex
git commit -m "Rewrite discussion and conclusion for sensitivity framing"
```

---

### Task 8: Update refs.bib and appendix.tex

**Files:**
- Modify: `paper/refs.bib`
- Modify: `paper/sections/appendix.tex`

**Step 1: Add missing bibliography entries**

Add to refs.bib:
```bibtex
@article{LinGuess1994,
    author = {Lin, D. K. J. and Guess, F. M.},
    title = {System Life Data Analysis with Dependent Partial Knowledge on the Exact Cause of System Failure},
    journal = {Microelectronics and Reliability},
    volume = {34},
    pages = {535--544},
    year = {1994}
}

@article{Guttman1995,
    author = {Guttman, I. and Lin, D. K. J. and Reiser, B. and Usher, J. S.},
    title = {Dependent Masking and System Life Data Analysis: {Bayesian} Inference for Two-Component Systems},
    journal = {Lifetime Data Analysis},
    volume = {1},
    pages = {87--100},
    year = {1995}
}

@article{Mukhopadhyay2006,
    author = {Mukhopadhyay, C. and Basu, A. P.},
    title = {Maximum Likelihood Analysis of Masked Series System Lifetime Data},
    journal = {Journal of Statistical Planning and Inference},
    volume = {136},
    pages = {803--838},
    year = {2006}
}

@article{CraiuReiser2010,
    author = {Craiu, R. V. and Reiser, B.},
    title = {About Conditional Masking Probability Models},
    journal = {Statistics \& Probability Letters},
    volume = {80},
    pages = {1536--1542},
    year = {2010}
}
```

Remove unused entries (KL-divergence, bootstrap, etc.) — but verify none are still cited first.

**Step 2: Update appendix**

- Keep proof of misspecification theorem (from §3.4).
- Add full simulation results tables (all 11 severity levels for both exponential and Weibull).
- Add software note: mdrelax package description and availability.
- Remove any C3/C1 appendix material.

**Step 3: Commit**

```bash
git add paper/refs.bib paper/sections/appendix.tex
git commit -m "Update bibliography and appendix for redesigned paper"
```

---

### Task 9: Final compilation and verification

**Files:**
- All paper files

**Step 1: Full LaTeX compilation**

```bash
cd paper && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

Expected: Compiles cleanly. Check for:
- No undefined references
- No missing citations
- No overfull hboxes (minor ones OK)
- Page count in target range (12-15)

**Step 2: Run package tests**

```bash
Rscript -e 'devtools::test()'
```

Expected: All 1276 tests pass (package code unchanged).

**Step 3: Review paper PDF**

Read through the compiled PDF. Verify:
- Figures render correctly
- Table numbers make sense
- Cross-references work
- No leftover references to C3, C1 relaxation, KL-divergence, etc.

**Step 4: Commit final state**

```bash
git add -A paper/
git commit -m "Complete paper redesign: sensitivity to non-informative masking assumption"
```
