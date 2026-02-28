# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R package (`mdrelax`) for relaxed candidate set models in masked series system failure data. Implements likelihood-based inference for exponential and Weibull series systems under relaxed C1/C2/C3 conditions, with minimal dependencies suitable for Monte Carlo simulation studies.

**Version:** 1.0.0 (per DESCRIPTION) — all core models + Weibull relaxed models implemented and tested.

**Related Projects:**
- `maskedcauses`: General masked data likelihood framework (validated code migrates there)
- `wei.series.md.c1.c2.c3`: Reference Weibull implementation
- Master's thesis: `/home/spinoza/github/papers/reliability-estimation-in-series-systems/`

## Development Commands

```r
devtools::load_all()                           # Load for development
devtools::test()                               # Run all tests (1276 tests)
devtools::test(filter="c1-c2-c3")              # C1-C2-C3 model tests
devtools::test(filter="c1-c3")                 # Relaxed C2 model tests
devtools::test(filter="c1-c2$")                # Relaxed C3 model tests
devtools::test(filter="relaxed-c1")            # Relaxed C1 model tests
devtools::test(filter="wei-series$")           # Weibull distribution tests
devtools::test(filter="wei-series-c1-c2-c3")   # Weibull C1-C2-C3 tests
devtools::test(filter="wei-series-c1-c3")      # Weibull relaxed C2 tests
devtools::test(filter="wei-series-c1-c2$")     # Weibull relaxed C3 tests
devtools::test(filter="wei-simulations")       # Weibull simulation tests
devtools::document()                           # Regenerate man/ and NAMESPACE
devtools::check()                              # Full R CMD check
```

### Paper Compilation

```bash
cd paper && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

### Simulation Pipeline

```r
# Full paper simulations (Studies 4-5: ~20 min)
Rscript paper/run_paper_simulations.R

# Earlier studies (Studies 1-3: ~30 min)
Rscript inst/simulations/run_all.R

# Generate figures and LaTeX tables from saved results
Rscript inst/simulations/generate_figures.R
Rscript inst/simulations/generate_latex_tables.R
```

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

**Naming convention:** Each model file exports `loglik_*`, `score_*`, `fim_*`, `mle_*`, `r*_md` (data gen), and `*_df` (data frame wrappers). All `loglik_*` functions return **closures** `function(theta)`, not values.

### Simulation Scenarios

**Exponential:** 1 (baseline), 2/2b (overfit C2), 3 (C2 data→C1-C2-C3 misspec), 4/4b (C2 correct), 5 (overfit C3), 6 (C3 data→C1-C2-C3 misspec), 6b (C3 correct)

**Weibull:** W1 (baseline), W2 (overfit C2), W3 (C2 misspec), W4 (C2 correct), W5 (overfit C3), W6 (C3 misspec), W7 (C3 correct)

### Paper Structure

```
paper/
├── main.tex                    # Master document
├── refs.bib                    # Bibliography
├── run_paper_simulations.R     # Studies 4-5 runner (saves to data/, figures to inst/)
├── data/                       # Simulation results (.rds, gitignored)
└── sections/
    ├── introduction.tex
    ├── background.tex
    ├── relaxed_models.tex      # C2/C3 relaxation theory
    ├── identifiability.tex     # Theorems 4.1-4.8
    ├── simulations.tex         # Studies 1-5 results
    ├── discussion.tex
    ├── conclusion.tex
    └── appendix.tex

inst/simulations/
├── run_all.R                   # Studies 1-3 master runner
├── generate_figures.R          # fig1-fig6 from saved results
├── generate_latex_tables.R     # LaTeX tables from saved results
├── figures/                    # fig1-fig10 (PDF + PNG)
├── tables/                     # Generated LaTeX tables
└── results/                    # Studies 1-3 results (.rds)
```

**Figure path:** `\graphicspath{{../inst/simulations/figures/}}` in main.tex — all figures (including paper-specific fig7-fig10) go in `inst/simulations/figures/`.

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
- **Simulation `test-simulations.R`:** 18 tests are `skip()`-ed (legacy framework). The 1276 passing tests exclude these.
