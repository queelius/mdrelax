# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

This R package (`mdrelax`) extends the work from the master’s thesis
“Reliability Estimation in Series Systems” (located at
`/home/spinoza/github/papers/reliability-estimation-in-series-systems/`)
by implementing relaxed candidate set models that go beyond the
traditional C1, C2, C3 conditions.

**Project Status:** Early draft (v0.9.1), currently unstable. The
package extends masked data analysis by allowing various candidate set
structures beyond strict C1/C2/C3 conditions.

**Related Projects (use as references):** - **Original Paper**
(`/home/spinoza/github/papers/reliability-estimation-in-series-systems/`):
Master’s thesis with extensive simulation studies for Weibull series
systems, BCa confidence intervals, and comprehensive validation -
**Model Selection Paper**
(`/home/spinoza/github/rlang/reliability-estimation-in-series-systems-model-selection/`):
Follow-up work on reduced models (homogeneous shape parameters) and
likelihood ratio tests - **FIM Precursor**
(`/home/spinoza/github/papers/expo-masked-fim/`): Closed-form Fisher
Information Matrix derivations for exponential series systems

## Development Commands

``` r
devtools::load_all()           # Load for interactive development
devtools::document()           # Generate documentation from roxygen
devtools::check()              # Run R CMD check
devtools::test()               # Run tests (tests need to be created)
covr::package_coverage()       # Test coverage report
pkgdown::build_site()          # Build documentation site
```

Data generation scripts in `data-raw/`:

``` r
source("data-raw/data.R")                    # Main data generation
source("data-raw/exp_series_stats_1.R")      # Exponential series statistics
```

## Core Dependencies

- `md.tools` (github::queelius/md.tools) - Masked data utilities for
  encoding/decoding
- `algebraic.mle` (github::queelius/algebraic.mle) - MLE framework with
  Fisher information support
- `cubature` - Numerical integration for complex likelihoods

## Mathematical Background

### Series Systems

A series system fails when ANY component fails. For $m$ components: -
System lifetime: $T_{i} = \min\{ T_{i1},T_{i2},\ldots,T_{im}\}$ - System
reliability:
$R(t;\theta) = \prod_{j = 1}^{m}R_{j}\left( t;\theta_{j} \right)$ -
System hazard:
$h(t;\theta) = \sum_{j = 1}^{m}h_{j}\left( t;\theta_{j} \right)$

### Masked Data Problem

**Observed per system:** 1. Right-censored system lifetime
$S_{i} = \min\left( T_{i},\tau_{i} \right)$ with indicator $\delta_{i}$
2. Candidate set $C_{i}$ when system fails (components that could have
caused failure) 3. Component failure times $T_{ij}$ are NOT observed

**Key challenge:** Failed component $K_{i}$ is latent but constrained to
$K_{i} \in C_{i}$

### Traditional Conditions (C1, C2, C3)

**Condition 1 (C1):** Failed component always in candidate set:
$\Pr\{ K_{i} \in C_{i}\} = 1$

**Condition 2 (C2):** Non-informative masking within candidate set
(equal probabilities for all candidates given system failure time)

**Condition 3 (C3):** Masking independent of system parameters θ

Under C1, C2, C3, the likelihood contribution for a failed system at
time $t$ with candidate set $c$ is:
$$L_{i}(\theta) \propto R(t;\theta) \times \sum\limits_{j \in c}h_{j}\left( t;\theta_{j} \right)$$

### Relaxed Models (This Package)

This package relaxes these conditions, particularly: - **Relaxed C1:**
Only enforce failed component in candidate set - **Informed masking:**
Probability depends on component failure ranks (simulation-friendly) -
**KL-divergence constraints:** Control “informativeness” relative to
baseline Bernoulli model

## Code Architecture

### Core Functions (R/md_candidate_set_models.R)

**Candidate Set Models:** -
[`md_bernoulli_cand_C1_C2_C3()`](https://queelius.github.io/mdrelax/reference/md_bernoulli_cand_C1_C2_C3.md) -
Traditional uninformed Bernoulli model -
[`md_bernoulli_cand_C1_kld()`](https://queelius.github.io/mdrelax/reference/md_bernoulli_cand_C1_kld.md) -
Relaxed model with KL-divergence from baseline -
[`md_cand_sampler()`](https://queelius.github.io/mdrelax/reference/md_cand_sampler.md) -
Generate candidate sets from probability vectors -
[`md_block_candidate_m3()`](https://queelius.github.io/mdrelax/reference/md_block_candidate_m3.md) -
Block model demonstrating non-identifiability

**Estimation:** -
[`md_mle_exp_series_C1_C2_C3()`](https://queelius.github.io/mdrelax/reference/md_mle_exp_series_C1_C2_C3.md)
/ `md_mle_weibull_series_C1_C2_C3()` - MLE computation -
[`md_loglike_exp_series_C1_C2_C3()`](https://queelius.github.io/mdrelax/reference/md_loglike_exp_series_C1_C2_C3.md)
/ `md_loglike_weibull_series_C1_C2_C3()` - Log-likelihood -
[`md_score_exp_series_C1_C2_C3()`](https://queelius.github.io/mdrelax/reference/md_score_exp_series_C1_C2_C3.md)
/ `md_score_weibull_series_C1_C2_C3()` - Score functions -
[`md_fim_exp_series_C1_C2_C3()`](https://queelius.github.io/mdrelax/reference/md_fim_exp_series_C1_C2_C3.md) -
Fisher information matrix

### Informed Masking (R/informed_candidate_set_utils.R)

- `informative_masking_by_rank(ts, alpha, beta)` - Rank-based
  probabilities
  - α → 0: uniform distribution over non-failed components
  - α → ∞: only failed component and rank-2 component included
  - β: maximum weight for rank-2 component
- `generate_relaxed_cand_C1(d, ts, p, ...)` - Find Q with target
  KL-divergence from P
- `kl_divergence_bernoulli(p, q)` - KL divergence between Bernoulli
  vectors

### Distribution Functions

Standard R naming (d/p/q/r): - Exponential: `dexp_series_system()`,
`pexp_series_system()`, `qexp_series_system()`, `rexp_series_system()` -
Weibull: `dweibull_series()`, `pweibull_series()`, `qweibull_series()`,
`rweibull_series()` - Hazard/Survival: `hazard_exp_series_system()`,
`survival_exp_series_system()`, etc.

### Data Encoding Pattern (via md.tools)

``` r
# Extract component times from prefixed columns (t1, t2, ..., tm → matrix)
Tm <- md_decode_matrix(md, "t")

# Encode probabilities back to data frame (matrix → q1, q2, ..., qm)
md <- md %>% bind_cols(md_encode_matrix(Q, "q"))

# Mark columns as latent
md <- md %>% md_mark_latent(paste0("q", 1:m))
```

## Typical Analysis Workflow

1.  **Generate/Load masked data** with component times (t1, t2, …, tm)
2.  **Apply candidate model** to add masking probabilities (q1, q2, …,
    qm):
    - `md_bernoulli_cand_C1_C2_C3(md, p=0.3)` for traditional model
    - `md_bernoulli_cand_C1_kld(md, p=0.3, d=0.5)` for relaxed model
3.  **Sample candidate sets** using `md_cand_sampler(md)` → creates (x1,
    x2, …, xm)
4.  **Compute MLE** using `md_mle_exp_series_C1_C2_C3(md)` or Weibull
    variant
5.  **Analyze** using standard MLE methods (confidence intervals,
    hypothesis tests)

## Known Issues and Development Notes

- **No test suite exists** - tests/ directory needs to be populated
- **Optimization failures** - Functions use tryCatch; convergence not
  guaranteed
- **Bug in informative_masking_by_rank():** Uses undefined
  `alpha0`/`beta0` instead of `alpha`/`beta` parameters (line 56-57)
- **grad_descent not defined:** `generate_relaxed_cand_C1` calls
  `grad_descent` but utils.R defines `grad_ascent`
- **Identifiability:** Some candidate set structures produce non-unique
  MLEs (see `md_block_candidate_m3`)

## Reference Implementation Comparison

The original thesis (`wei.series.md.c1.c2.c3` package) provides
validated implementations. Key simulation findings from that work: - MLE
performs well even with significant masking (p up to 0.4) and censoring
(up to 40%) - BCa confidence intervals have good coverage probability -
Shape parameters harder to estimate than scale parameters - Larger
samples (n ≥ 250) provide reliable estimates
