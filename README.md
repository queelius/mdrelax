# mdrelax

Relaxed Candidate Set Models for Masked Data in Series Systems

## Overview

This R package implements likelihood-based inference for series systems with masked failure data under relaxed candidate set conditions. It provides a complete framework for analyzing reliability data when the standard assumptions (C1, C2, C3) about candidate set formation may not hold.

**Model hierarchy:**

```
C1-C2-C3 (standard)
  |-- Relaxed C2: Informative masking (P matrix)
  |-- Relaxed C3: Parameter-dependent masking (power-weighted alpha)
  '-- Relaxed C1: Failed component may not be in candidate set
```

Both exponential and Weibull component lifetime distributions are supported.

## Installation

```r
# Install from GitHub
remotes::install_github("queelius/mdrelax")
```

## Quick Start

```r
library(mdrelax)

# --- Exponential series system (standard C1-C2-C3) ---
set.seed(42)
sim <- rexp_series_md(n = 200, theta = c(1, 1.5, 2), p = 0.3, tau = 5)
fit <- mle_exp_series(sim$t, sim$C, sim$delta)
fit$theta  # Estimated rates
fit$se     # Standard errors

# --- Relaxed C2 (informative masking) ---
P <- make_P_matrix(3, "full", values = c(0.2, 0.4, 0.5, 0.3, 0.6, 0.35))
sim <- rexp_series_md_c1_c3(n = 200, theta = c(1, 1.5, 2), P = P, tau = 5)
fit <- mle_exp_series_c1_c3(sim$t, sim$C, sim$delta, fixed_P = P)

# --- Weibull series system ---
sim <- rwei_series_md(n = 200, shapes = c(2, 1.5), scales = c(3, 4),
                       p = 0.3, tau = 8)
fit <- mle_wei_series(sim$t, sim$C, sim$delta)
fit$shapes  # Estimated shape parameters
fit$scales  # Estimated scale parameters
```

## Background

In series systems with masked failure data:
- The system fails when any component fails
- The failed component is not directly observed
- A *candidate set* of possible failed components is reported

Traditional analysis assumes three conditions:
- **C1:** The failed component is always in the candidate set
- **C2:** Masking probabilities do not depend on which component failed (non-informative)
- **C3:** Masking probabilities do not depend on the model parameters

This package implements inference under each relaxation:

| Model | Conditions | Key parameters |
|-------|-----------|----------------|
| Standard | C1 + C2 + C3 | theta (rates) |
| Relaxed C2 | C1 + C3 | theta + P matrix |
| Relaxed C3 | C1 + C2 | theta + alpha (power weight) |
| Relaxed C1 | None required | theta + P matrix (diag < 1) |

## Simulation Studies

The package includes a framework for comparing model robustness:

```r
# Compare scenarios: what happens under model misspecification?
study <- quick_simulation_study(n_sim = 100, n_obs = 200,
                                 theta = c(1, 2), p = 0.3)
print_simulation_summary(study$summary)
```

Key scenarios test:
- Fitting the wrong model when conditions are violated
- Cost of using a more flexible model when simpler conditions hold
- Parameter recovery under each relaxation

## Key Features

- Maximum likelihood estimation for exponential and Weibull series systems
- Analytical score functions and Fisher information matrices
- Standard error computation via observed Fisher information
- Data generation for simulation studies under all model variants
- Data frame interface for interoperability
- Minimal dependencies (only `stats`)

## Documentation

See the [package website](https://queelius.github.io/mdrelax/) for full documentation.

## Related Work

- [wei.series.md.c1.c2.c3](https://github.com/queelius/wei.series.md.c1.c2.c3) - Weibull series systems under C1-C2-C3
- [likelihood.model.series.md](https://github.com/queelius/likelihood.model.series.md) - General likelihood framework

## License

GPL (>= 3)
