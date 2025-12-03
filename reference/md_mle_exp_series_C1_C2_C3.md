# Maximum likelihood estimator for exponential series system (C1,C2,C3)

Computes the MLE of component failure rates for an exponential series
system with masked data under C1, C2, C3 candidate set conditions.

## Usage

``` r
md_mle_exp_series_C1_C2_C3(
  md,
  theta0 = NULL,
  sysvar = "t",
  setvar = "x",
  deltavar = "delta",
  use_annealing = FALSE,
  control = list(),
  lower = 1e-10,
  upper = Inf,
  hessian = TRUE
)
```

## Arguments

- md:

  Masked data frame containing:

  - System lifetime column (default: "t")

  - Candidate set columns (default: "x1", "x2", ..., "xm")

  - Optional right-censoring indicator (default: "delta")

- theta0:

  Initial parameter values for optimization. If NULL, uses the method of
  moments estimator based on total system hazard.

- sysvar:

  Column name for system lifetime. Default is "t".

- setvar:

  Column prefix for candidate set (Boolean matrix). Default is "x".

- deltavar:

  Column name for right-censoring indicator. Default is "delta". If NULL
  or column doesn't exist, assumes no censoring.

- use_annealing:

  Logical. If TRUE, uses simulated annealing to find a good starting
  point before local optimization. Default is FALSE.

- control:

  List of control parameters passed to optim(). Default uses L-BFGS-B
  with reasonable settings.

- lower:

  Lower bounds for parameters. Default is 1e-10.

- upper:

  Upper bounds for parameters. Default is Inf.

- hessian:

  Logical. If TRUE, compute Hessian at the MLE. Default is TRUE.

## Value

An object of class `mle` (from `algebraic.mle` package) containing:

- `theta.hat`: The MLE of rate parameters

- `sigma`: Estimated variance-covariance matrix

- `loglike`: Log-likelihood at the MLE

- `info`: Observed Fisher information matrix

- `nobs`: Number of observations

- `converged`: Whether optimization converged

## Details

Uses L-BFGS-B optimization with optional simulated annealing for finding
a good starting point.

The MLE for exponential series systems with masked data has a
closed-form solution only in special cases. In general, numerical
optimization is required. This function uses the L-BFGS-B algorithm with
analytically computed gradients for efficiency.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate masked data
library(tibble)
set.seed(42)
n <- 100
theta_true <- c(1, 1.5, 2)
# ... (data generation code)

mle <- md_mle_exp_series_C1_C2_C3(md, theta0 = theta_true)
print(mle)
confint(mle)
} # }
```
