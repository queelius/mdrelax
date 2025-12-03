# Generate masked data for exponential series system (C1,C2,C3)

Generates synthetic masked failure data for an exponential series system
using the Bernoulli candidate set model satisfying C1, C2, C3
conditions.

## Usage

``` r
md_exp_series_system_bernoulli_cand_C1_C2_C3(n, rates, p, tau)
```

## Arguments

- n:

  Sample size (number of systems).

- rates:

  Vector of exponential rate parameters for each component.

- p:

  Probability that each non-failed component is in the candidate set.

- tau:

  Right-censoring time. Can be scalar (same for all) or vector.

## Value

A masked data frame with columns:

- t1, t2, ..., tm: Component failure times (latent)

- t: System failure time

- k: Index of failed component (latent)

- delta: Censoring indicator (TRUE = censored)

- q1, ..., qm: Candidate set probabilities

- x1, ..., xm: Candidate set indicators

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
md <- md_exp_series_system_bernoulli_cand_C1_C2_C3(
  n = 100,
  rates = c(1, 1.5, 2),
  p = 0.3,
  tau = 2.0
)
} # }
```
