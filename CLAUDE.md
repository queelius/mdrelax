# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R package (`mdrelax`) for relaxed candidate set models in masked series system failure data. Implements likelihood-based inference for exponential and Weibull series systems under relaxed C1/C2/C3 conditions, with minimal dependencies suitable for Monte Carlo simulation studies.

The R package is a self-contained methodology library. The companion **paper in `paper/`** is a separate artifact: a Technometrics submission titled *"Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption"*. The paper uses the package, but the package's API is not driven by the paper's narrative; keep them logically separate when making changes.

**Version:** 1.0.0 (per DESCRIPTION) — all core models + Weibull relaxed models implemented and tested.

**Related Projects:**
- `maskedcauses`: General masked data likelihood framework (validated code migrates there)
- `wei.series.md.c1.c2.c3`: Reference Weibull implementation
- Master's thesis: `/home/spinoza/github/papers/reliability-estimation-in-series-systems/`

## Development Commands

```r
devtools::load_all()                           # Load for development
devtools::test()                               # Run all tests (~240 test_that blocks across 11 files)
devtools::test(filter="c1-c2-c3")              # C1-C2-C3 model tests
devtools::test(filter="c1-c3")                 # Relaxed C2 model tests
devtools::test(filter="c1-c2$")                # Relaxed C3 model tests
devtools::test(filter="relaxed-c1")            # Relaxed C1 model tests
devtools::test(filter="wei-series$")           # Weibull distribution tests
devtools::test(filter="wei-series-c1-c2-c3")   # Weibull C1-C2-C3 tests
devtools::test(filter="wei-series-c1-c3")      # Weibull relaxed C2 tests
devtools::test(filter="wei-series-c1-c2$")     # Weibull relaxed C3 tests
devtools::test(filter="wei-simulations")       # Weibull simulation tests
devtools::test(filter="relaxed-expectations")  # Weibull TDD oracle (theoretical properties)
devtools::document()                           # Regenerate man/ and NAMESPACE
devtools::check()                              # Full R CMD check
```

### Paper Compilation

The `paper/Makefile` is the canonical build path. Direct `pdflatex` invocation still works but is no longer the recommended entry.

```bash
cd paper && make            # PDF: pdflatex, bibtex (auto-detected), pdflatex x2
cd paper && make html       # HTML version via tex2html into paper/html_paper/
cd paper && make clean      # Strip aux/log/bbl/out and html_paper/
```

### Simulation Pipeline

The paper has two simulation runners. The **sensitivity sweep** is the centerpiece of the current paper; the older Studies 1 to 5 framework remains for reference and figures used in earlier sections.

```r
# Paper centerpiece: 5-component Weibull C2 sensitivity sweep (~30-60 min, B=200)
# Generates fig_bias_vs_severity.pdf, fig_rmse_vs_severity.pdf, fig_coverage_vs_severity.pdf
Rscript paper/run_sensitivity_sweep.R              # full
Rscript paper/run_sensitivity_sweep.R --quick      # smoke test

# C3 sensitivity sweep on same 5-component Weibull system (~30-60 min, B=200)
# Tests first-order total-hazard preservation for Weibull under C3 violation
# Generates fig_c3_bias_vs_alpha.pdf, fig_c3_rmse_vs_alpha.pdf, fig_c3_coverage_vs_alpha.pdf
Rscript paper/run_sensitivity_sweep_c3.R           # full
Rscript paper/run_sensitivity_sweep_c3.R --quick   # smoke test

# All three sweeps accept --n=N and --B=B for size-aware ablation runs.
# Output filenames carry an _n{N} suffix when N != 500 to preserve baseline outputs.
Rscript paper/run_sensitivity_sweep_c1.R --n=2000 --B=100   # larger-sample C1 sweep

# Studies 4-5 (Exponential C3 misspec + Weibull baseline scenarios, ~20 min)
Rscript paper/run_paper_simulations.R [--quick] [--study=4|5|all]

# Earlier exploratory studies (Studies 1-3: ~30 min)
Rscript inst/simulations/run_all.R

# Figures and LaTeX tables for Studies 1-3 (only)
Rscript inst/simulations/generate_figures.R
Rscript inst/simulations/generate_latex_tables.R
```

**Output convention:** all paper figures land in `inst/simulations/figures/` because `main.tex` sets `\graphicspath{{../inst/simulations/figures/}}`. The sensitivity sweep writes `.rds` artifacts to `paper/data/` (gitignored). Do not introduce a parallel figure directory under `paper/`.

## Code Architecture

### Model Hierarchy

```
C1-C2-C3 (standard)
├── Relaxed C2 (informative masking, P matrix)
├── Relaxed C3 (parameter-dependent masking, alpha power weight)
└── Relaxed C1 (P(K in C) < 1, sensitivity analysis)
```

Each model tier has both Exponential and Weibull implementations.

### Source Files (R/)

| File | Model | Key Pattern |
|------|-------|-------------|
| `exp_series_c1_c2_c3.R` | Standard exponential | `loglik_exp_series`, `score_exp_series`, `fim_exp_series`, `mle_exp_series`, `rexp_series_md` |
| `exp_series_c1_c3.R` | Exponential relaxed C2 | P matrix masking; `make_P_matrix`, `compute_pi`, `satisfies_C2` |
| `exp_series_c1_c2.R` | Exponential relaxed C3 | Power weights: `p_j(θ) = base_p * θ_j^α / max(θ^α)` |
| `exp_series_relaxed_c1.R` | Exponential relaxed C1 | Likelihood sums over ALL k=1,...,m |
| `wei_series.R` | Weibull distributions | d/p/q/r/hazard/surv functions |
| `wei_series_c1_c2_c3.R` | Weibull standard | Same API pattern as exponential |
| `wei_series_c1_c3.R` | Weibull relaxed C2 | P matrix with Weibull hazards |
| `wei_series_c1_c2.R` | Weibull relaxed C3 | Power weights: `p_j(θ) = base_p * (k_j/λ_j)^α / max(...)` |
| `simulation_study.R` | Simulation framework | Scenarios 1-6b (exp) + W1-W7 (Weibull) |
| `robustness_intervals.R` | Applied robustness tool | `ri_simulation()`, `ri_first_order()`: max severity at which estimand stays within tolerance |

**Naming convention:** Each model file exports `loglik_*`, `score_*`, `fim_*`, `mle_*`, `r*_md` (data gen), and `*_df` (data frame wrappers). All `loglik_*` functions return **closures** `function(theta)`, not values.

### Simulation Scenarios

**Exponential:** 1 (baseline), 2/2b (overfit C2), 3 (C2 data→C1-C2-C3 misspec), 4/4b (C2 correct), 5 (overfit C3), 6 (C3 data→C1-C2-C3 misspec), 6b (C3 correct)

**Weibull:** W1 (baseline), W2 (overfit C2), W3 (C2 misspec), W4 (C2 correct), W5 (overfit C3), W6 (C3 misspec), W7 (C3 correct)

### Paper Structure

The paper was redesigned in February 2026 as a focused sensitivity analysis. Two `.tex` files in `paper/sections/` (`relaxed_models.tex`, `identifiability.tex`) are pre-redesign remnants. They still exist on disk but are NOT `\input` from `main.tex`; do not edit them expecting changes to land in the PDF.

```
paper/
├── Makefile                    # Canonical build (pdf + html targets)
├── main.tex                    # Master document; \input order = paper section order
├── refs.bib                    # Bibliography
├── run_sensitivity_sweep.R     # Paper centerpiece: 5-component Weibull sweep
├── run_paper_simulations.R     # Studies 4-5 (exp C3 misspec + Weibull baseline)
├── data/                       # Simulation results (.rds, gitignored)
├── html_paper/                 # tex2html build output (gitignored)
└── sections/
    ├── introduction.tex
    ├── background.tex
    ├── sensitivity_framework.tex   # Likelihood under C1; Bernoulli model; misspec theorem; non-identifiability
    ├── simulations.tex             # Sensitivity sweep results (severity vs bias/RMSE/coverage)
    ├── application.tex             # Guo et al. turbine engine real-data application
    ├── discussion.tex
    ├── conclusion.tex
    ├── appendix.tex                # Score derivations, non-identifiability evidence, software
    ├── relaxed_models.tex          # NOT INCLUDED in main.tex (pre-redesign)
    └── identifiability.tex         # NOT INCLUDED in main.tex (pre-redesign)

inst/simulations/
├── run_all.R                   # Older Studies 1-3 master runner
├── sim_kl_efficiency.R         # Study 1: KL-divergence and efficiency
├── sim_misspecification.R      # Study 2: Misspecification bias
├── sim_identifiability.R       # Study 3: Identifiability via FIM
├── generate_figures.R          # fig1-fig6 from saved Studies 1-3 results
├── generate_latex_tables.R     # LaTeX tables for Studies 1-3
├── figures/                    # ALL paper figures (fig1..fig10 + fig_*_severity)
├── tables/                     # Generated LaTeX tables
└── results/                    # Studies 1-3 results (.rds)
```

**Real-data application:** `data/guo_weibull_series_md.rda` and `guo_weibull_series_mle.rda` are `LazyData` package datasets used in `paper/sections/application.tex`. Source PDFs of the Guo et al. and Usher reference papers live in `inst/*.pdf` (not used at runtime, kept for citation reference).

**Project state file:** `.papermill.md` tracks paper metadata (thesis, venue, prior-art status, review history). Update it when the paper's central claim, venue, or experimental scope changes.

## Key Formulas

**Exponential C1-C2-C3:**
```
l(θ) = -Σᵢ tᵢ·Σⱼ θⱼ + Σᵢ:δᵢ=1 log(Σⱼ∈Cᵢ θⱼ)
```

**Relaxed C2 (both exp/Weibull):** adds `πⱼ(cᵢ)` weights from P matrix inside the log-sum.

**Relaxed C3 (both exp/Weibull):** adds `log π_c(θ)` term where inclusion probs depend on θ.

**Weibull C1-C2-C3:**
```
l(θ) = Σᵢ [-Σⱼ (tᵢ/λⱼ)^kⱼ] + Σᵢ:δᵢ=1 log(Σⱼ∈Cᵢ hⱼ(tᵢ))
```

## Conventions

- **Censoring:** `delta=1` = observed failure, `delta=0` = right-censored at `tau`
- **Candidate set:** `C[i,j]=TRUE` means component j is in observation i's candidate set
- **Exponential parameters:** `theta` = rate vector `(λ₁, ..., λₘ)`
- **Weibull parameters:** `theta = c(k1, β1, k2, β2, ...)` — interleaved shape-scale pairs (not grouped)
- **P matrix:** `P[j,k] = P(j in C | K = k)`. Diagonal = 1 under C1. Rows = included component, cols = true cause.
- **Functions return closures:** `loglik_*(t, C, delta, ...)` returns `function(theta)`, not a value
- **Only `stats` in Imports:** Package is deliberately dependency-minimal for Monte Carlo use

## Gotchas

- **Weibull theta ordering:** Parameters interleave as `c(k1, β1, k2, β2, ...)`, NOT `c(k1, k2, ..., β1, β2, ...)`. Getting this wrong silently produces wrong likelihoods.
- **Exponential identifiability:** From system-level data alone, only `sum(θ)` is identifiable — individual rates require masking information. Tests check `sum(theta_hat)` not individual components.
- **Score = numerical gradient:** Test pattern verifies analytical score against `numDeriv::grad`. All model files follow this cross-validation pattern.
- **L-BFGS-B bounds:** Weibull MLE uses `lower = 1e-6` bounds. If optimization hits bounds, the result may be unreliable (check `converged` flag).
- **Simulation `test-simulations.R`:** 18 tests are `skip()`-ed (legacy framework). The passing test counts exclude these.
- **Theoretical-expectations file:** `test-wei-series-relaxed-expectations.R` is a TDD oracle, not a regression suite. Each `test_that` documents a statistical property the implementation MUST satisfy (score zero at MLE, FIM positive-definite, etc.). Treat failures here as evidence the math is wrong, not as flaky tests.
