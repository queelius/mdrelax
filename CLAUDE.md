# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R package (`mdrelax`) for relaxed candidate set models in masked series system failure data. Implements likelihood-based inference for exponential and Weibull series systems under relaxed C1/C2/C3 conditions, with minimal dependencies suitable for Monte Carlo simulation studies.

**Status:** v1.1.0 — all core models + Weibull relaxed models implemented and tested (1276 tests).

**Related Projects:**
- `likelihood.model.series.md`: General likelihood framework (will receive validated code after stabilization)
- `wei.series.md.c1.c2.c3`: Reference Weibull implementation
- Master's thesis: `/home/spinoza/github/papers/reliability-estimation-in-series-systems/`

## Development Commands

```r
devtools::load_all()                           # Load for development
devtools::test()                               # Run all tests (1276 tests)
devtools::test(filter="c1-c2-c3")              # C1-C2-C3 model tests (153)
devtools::test(filter="c1-c3")                 # Relaxed C2 model tests (180)
devtools::test(filter="c1-c2$")                # Relaxed C3 model tests (132)
devtools::test(filter="relaxed-c1")            # Relaxed C1 model tests (228)
devtools::test(filter="wei-series$")           # Weibull distribution tests (31)
devtools::test(filter="wei-series-c1-c2-c3")   # Weibull C1-C2-C3 tests (144)
devtools::test(filter="wei-series-c1-c3")      # Weibull relaxed C2 tests (162)
devtools::test(filter="wei-series-c1-c2$")     # Weibull relaxed C3 tests (177)
devtools::test(filter="wei-series-relaxed")    # Theoretical expectations (28)
devtools::test(filter="wei-simulations")       # Weibull simulation tests (41)
devtools::document()                           # Regenerate man/ and NAMESPACE
devtools::check()                              # Full R CMD check
```

## Code Architecture

### Model Hierarchy

```
C1-C2-C3 (standard)
├── Relaxed C2 (informative masking, P matrix)
├── Relaxed C3 (parameter-dependent masking, alpha power weight)
└── Relaxed C1 (P(K in C) < 1, sensitivity analysis)
```

Each model tier: Exponential + Weibull.

### R/exp_series_c1_c2_c3.R — Standard Model (C1, C2, C3)

Dependency-free exponential series system under full conditions.

| Function | Description |
|----------|-------------|
| `loglik_exp_series(t, C, delta)` | Log-likelihood (returns closure) |
| `score_exp_series(t, C, delta)` | Score / gradient function |
| `fim_exp_series(t, C, delta)` | Fisher information matrix |
| `mle_exp_series(t, C, delta, theta0)` | MLE (returns list: theta, se, vcov, loglik, converged, fim) |
| `rexp_series_md(n, theta, p, tau)` | Data generation (Bernoulli masking) |
| `as_dataframe(sim)` | Convert sim output to data frame |
| `decode_matrix(df, prefix)` | Extract x1,x2,... columns into matrix |
| `encode_matrix(mat, prefix)` | Convert matrix to x1,x2,... columns |
| `*_df()` variants | Data frame wrappers for all above |

### R/exp_series_c1_c3.R — Relaxed C2 (Informative Masking)

General Bernoulli model where `P[j,k] = P(j in C | K = k)` can vary with k.

| Function | Description |
|----------|-------------|
| `loglik_exp_series_c1_c3(t, C, delta, P)` | Log-likelihood (P=NULL for joint estimation) |
| `score_exp_series_c1_c3(t, C, delta, P)` | Analytical score (P must be known) |
| `fim_exp_series_c1_c3(t, C, delta, P)` | Analytical FIM (P must be known) |
| `mle_exp_series_c1_c3(t, C, delta, ..., fixed_P)` | MLE with optional joint P estimation |
| `rexp_series_md_c1_c3(n, theta, P, tau)` | Data generation with P matrix |
| `make_P_matrix(m, type, p, values)` | Create P matrix ("uniform", "symmetric", "full") |
| `satisfies_C2(P)` | Check if P satisfies condition C2 |
| `compute_pi(c, k, P)` | Compute P(C=c \| K=k) |
| `compute_pi_all(c, P)` | All pi_k values for candidates in c |

### R/exp_series_c1_c2.R — Relaxed C3 (Parameter-Dependent Masking)

Power-weighted hazard model: `p_j(theta) = base_p * theta_j^alpha / max(theta^alpha)`.

| Function | Description |
|----------|-------------|
| `loglik_exp_series_c1_c2(t, C, delta, alpha, base_p)` | Log-likelihood (alpha=NULL for joint) |
| `score_exp_series_c1_c2(t, C, delta, alpha, base_p)` | Score (numerical masking gradient) |
| `fim_exp_series_c1_c2(t, C, delta, alpha, base_p)` | Numerical FIM |
| `mle_exp_series_c1_c2(t, C, delta, alpha, base_p)` | MLE with fixed or joint alpha |
| `rexp_series_md_c1_c2(n, theta, alpha, base_p, tau)` | Data generation |

### R/exp_series_relaxed_c1.R — Relaxed C1 (P(K in C) < 1)

Failed component may not be in candidate set. Likelihood sums over ALL k=1,...,m.

| Function | Description |
|----------|-------------|
| `loglik_exp_series_relaxed_c1(t, C, delta, P)` | Log-likelihood |
| `score_exp_series_relaxed_c1(t, C, delta, P)` | Analytical score |
| `fim_exp_series_relaxed_c1(t, C, delta, P)` | Analytical FIM |
| `mle_exp_series_relaxed_c1(t, C, delta, P)` | MLE |
| `rexp_series_md_relaxed_c1(n, theta, P, tau)` | Data generation (diagonal of P < 1) |

### R/wei_series.R — Weibull Distribution Functions

R-style d/p/q/r functions for Weibull series systems.

| Function | Description |
|----------|-------------|
| `dwei_series(t, shapes, scales)` | PDF |
| `pwei_series(q, shapes, scales)` | CDF |
| `qwei_series(p, shapes, scales)` | Quantile (Newton's method) |
| `rwei_series(n, shapes, scales)` | Random generation |
| `hazard_wei_series(t, shapes, scales)` | Hazard function |
| `surv_wei_series(t, shapes, scales)` | Survival function |
| `wei_series_mttf(shapes, scales)` | Mean time to failure |

### R/wei_series_c1_c2_c3.R — Weibull Masked Data (C1, C2, C3)

Parameter convention: `theta = c(k1, beta1, k2, beta2, ..., km, betam)`.

| Function | Description |
|----------|-------------|
| `loglik_wei_series(t, C, delta)` | Log-likelihood |
| `score_wei_series(t, C, delta)` | Analytical score |
| `fim_wei_series(t, C, delta)` | Numerical FIM (4-point Hessian) |
| `mle_wei_series(t, C, delta, theta0)` | MLE via L-BFGS-B |
| `rwei_series_md(n, shapes, scales, p, tau)` | Data generation |

### R/wei_series_c1_c3.R — Weibull Relaxed C2 (Informative Masking)

Weibull series with general Bernoulli model: `P[j,k] = P(j in C | K = k)` can vary with k.

| Function | Description |
|----------|-------------|
| `loglik_wei_series_c1_c3(t, C, delta, P)` | Log-likelihood (weighted hazard) |
| `score_wei_series_c1_c3(t, C, delta, P)` | Analytical score |
| `fim_wei_series_c1_c3(t, C, delta, P)` | Numerical FIM (4-point Hessian) |
| `mle_wei_series_c1_c3(t, C, delta, P, theta0)` | MLE via L-BFGS-B |
| `rwei_series_md_c1_c3(n, shapes, scales, P, tau)` | Data generation with P matrix |
| `*_df()` variants | Data frame wrappers |

### R/wei_series_c1_c2.R — Weibull Relaxed C3 (Parameter-Dependent Masking)

Power-weighted model: `p_j(θ) = base_p * (k_j/λ_j)^α / max((k_l/λ_l)^α)`.

| Function | Description |
|----------|-------------|
| `wei_power_weights(shapes, scales, alpha)` | Compute normalized power weights |
| `loglik_wei_series_c1_c2(t, C, delta, alpha, base_p)` | Log-likelihood (alpha=NULL for joint) |
| `score_wei_series_c1_c2(t, C, delta, alpha, base_p)` | Score (numerical masking gradient) |
| `fim_wei_series_c1_c2(t, C, delta, alpha, base_p)` | Numerical FIM |
| `mle_wei_series_c1_c2(t, C, delta, alpha, base_p)` | MLE with fixed or joint alpha |
| `rwei_series_md_c1_c2(n, shapes, scales, alpha, base_p, tau)` | Data generation |
| `*_df()` variants | Data frame wrappers |

### R/simulation_study.R — Simulation Framework

#### Exponential Simulation Functions

| Function | Description |
|----------|-------------|
| `sim_config(n_sim, n_obs, theta, tau)` | Create configuration |
| `run_simulation_scenario(config, scenario, ...)` | Run single scenario |
| `run_all_simulations(config, scenarios, ...)` | Run multiple scenarios |
| `summarize_simulation(results)` | Compute bias, variance, MSE, coverage |
| `print_simulation_summary(summary_df)` | Print formatted results |
| `quick_simulation_study(...)` | Quick convenience wrapper |

**Exponential Scenarios:** 1 (C1-C2-C3 baseline), 2 (C1-C2-C3 data, relaxed C2 joint), 2b (known P), 3 (relaxed C2 data, C1-C2-C3 misspecified), 4 (relaxed C2 data, relaxed C2 joint), 4b (known P), 5 (C1-C2-C3 data, relaxed C3 misspecified), 6 (relaxed C3 data, C1-C2-C3 misspecified), 6b (known alpha).

#### Weibull Simulation Functions

| Function | Description |
|----------|-------------|
| `wei_sim_config(n_sim, n_obs, shapes, scales, tau)` | Create Weibull configuration |
| `run_wei_replication(config, scenario, ...)` | Run single Weibull replication |
| `run_weibull_simulation_study(...)` | Run full Weibull simulation study |

**Weibull Scenarios:** W1 (baseline), W2 (C1-C2-C3→relaxed C2, overfit), W3 (relaxed C2→C1-C2-C3, misspec), W4 (relaxed C2→relaxed C2, correct), W5 (C1-C2-C3→relaxed C3, overfit), W6 (relaxed C3→C1-C2-C3, misspec), W7 (relaxed C3→relaxed C3, correct).

## Key Formulas

**Log-likelihood (C1-C2-C3, exponential):**
```
l(theta) = -sum_i t_i * sum_j theta_j + sum_{i:delta_i=1} log(sum_{j in C_i} theta_j)
```

**Log-likelihood (C1-C3, relaxed C2):**
```
l(theta) = -sum_i t_i * sum_j theta_j + sum_{i:delta_i=1} log(sum_{j in C_i} theta_j * pi_j(c_i))
```
where `pi_j(c_i) = P(C=c_i | K=j)` from the Bernoulli model with P matrix.

**Log-likelihood (Weibull C1-C2-C3):**
```
l(theta) = sum_i [-sum_j (t_i/beta_j)^k_j] + sum_{i:delta_i=1} log(sum_{j in C_i} h_j(t_i))
```

## Key Formulas (Weibull Relaxed Models)

**Log-likelihood (Weibull C1-C3, relaxed C2):**
```
l(θ) = Σᵢ [-Σⱼ (tᵢ/λⱼ)^kⱼ] + Σᵢ:δᵢ=1 log(Σⱼ∈Cᵢ hⱼ(tᵢ) · πⱼ(cᵢ))
```
where `πⱼ(c) = P(C=c|K=j)` from the general Bernoulli model with P matrix.

**Log-likelihood (Weibull C1-C2, relaxed C3):**
```
l(θ) = Σᵢ [-Σⱼ (tᵢ/λⱼ)^kⱼ] + Σᵢ:δᵢ=1 [log(Σⱼ∈Cᵢ hⱼ(tᵢ)) + log(π_c(θ))]
```
where inclusion prob `pⱼ(θ) = base_p · (kⱼ/λⱼ)^α / max((kₗ/λₗ)^α)`.

## Tests

Run with `devtools::test()` (1276 tests total). Test files:

| File | Tests | Model |
|------|-------|-------|
| `test-exp-series-c1-c2-c3.R` | 153 | Standard exponential |
| `test-exp-series-c1-c3.R` | 180 | Exponential relaxed C2 |
| `test-exp-series-c1-c2.R` | 132 | Exponential relaxed C3 |
| `test-exp-series-relaxed-c1.R` | 228 | Exponential relaxed C1 |
| `test-wei-series.R` | 31 | Weibull distributions |
| `test-wei-series-c1-c2-c3.R` | 144 | Weibull C1-C2-C3 |
| `test-wei-series-c1-c3.R` | 162 | Weibull relaxed C2 |
| `test-wei-series-c1-c2.R` | 177 | Weibull relaxed C3 |
| `test-wei-series-relaxed-expectations.R` | 28 | Theoretical expectations |
| `test-wei-simulations.R` | 41 | Weibull simulation framework |
| `test-simulations.R` | 18 (skip) | Legacy simulation framework |

Tests verify: likelihood/score/FIM correctness, score matches numerical gradient, FIM equals negative Hessian, MLE convergence and unbiasedness (Monte Carlo), model nesting (exponential as Weibull special case, relaxed→standard reduction), misspecification bias, edge cases.

## Conventions

- **Censoring:** `delta=1` = observed failure, `delta=0` = right-censored
- **Candidate set:** `C[i,j]=TRUE` means component j is in observation i's candidate set
- **Exponential parameters:** `theta` = rate parameters (lambda_1, ..., lambda_m)
- **Weibull parameters:** `theta = c(k1, beta1, k2, beta2, ...)` shape-scale pairs
- **P matrix:** `P[j,k] = P(j in C | K = k)`. Diagonal = 1 under C1.
- **Functions return closures:** `loglik_*(...)` returns `function(theta)`, not a value

## Roadmap

1. Done: Exponential C1-C2-C3 (153 tests)
2. Done: Relaxed C2 — General Bernoulli (180 tests)
3. Done: Relaxed C3 — Power-weighted masking (132 tests)
4. Done: Relaxed C1 — P(K in C) < 1 (228 tests)
5. Done: Weibull series systems (175 tests)
6. Done: Simulation framework with all scenarios
7. Done: Weibull + relaxed C2 model (162 tests)
8. Done: Weibull + relaxed C3 model (177 tests)
9. Done: Weibull simulation scenarios W1-W7 (41 tests)
10. Done: Theoretical expectation tests (28 tests)
11. Pending: Migrate validated code to likelihood.model.series.md
