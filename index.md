# mdrelax

Relaxed Candidate Set Models for Masked Data in Series Systems

## Overview

This R package implements likelihood-based inference for series systems
with masked failure data when traditional conditions are relaxed. It
extends the standard C1-C2-C3 framework by allowing:

- **Informative masking (relaxed C2):** Candidate set probabilities can
  depend on which component failed
- **Parameter-dependent masking (relaxed C3):** Masking mechanism can
  depend on system parameters
- **Flexible candidate set models:** Bernoulli, rank-based, and
  KL-divergence constrained models

## Installation

``` r
# Install from GitHub
remotes::install_github("queelius/mdrelax")
```

## Key Features

- Maximum likelihood estimation for exponential and Weibull series
  systems
- Fisher information matrix computation for efficiency analysis
- Informative masking models (rank-based, KL-constrained)
- Identifiability analysis tools
- Simulation utilities for Monte Carlo studies

## Quick Start

``` r
library(mdrelax)

# Generate masked data with Bernoulli candidate sets
md <- md_bernoulli_cand_C1_C2_C3(data, p = 0.3)

# Sample candidate sets
md <- md_cand_sampler(md)

# Compute MLE for exponential series system
fit <- md_mle_exp_series_C1_C2_C3(md)

# Get Fisher information matrix
fim <- md_fim_exp_series_C1_C2_C3(md, params(fit))
```

## Background

In series systems with masked failure data: - The system fails when any
component fails - The failed component is not directly observed - A
*candidate set* of possible failed components is reported

Traditional analysis assumes: - **C1:** Failed component is always in
the candidate set - **C2:** Non-informative masking (uniform probability
within candidate set) - **C3:** Masking independent of system parameters

This package provides tools for inference when C2 and/or C3 are
violated.

## Documentation

See the [package website](https://queelius.github.io/mdrelax/) for full
documentation.

## Related Work

- [wei.series.md.c1.c2.c3](https://github.com/queelius/wei.series.md.c1.c2.c3) -
  Original implementation for Weibull series systems under C1-C2-C3

## License

GPL (\>= 3)
