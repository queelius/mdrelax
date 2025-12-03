# Log-likelihood function for exponential series system with masked data (C1,C2,C3)

Returns a log-likelihood function for an exponential series system with
masked component cause of failure under candidate sets satisfying
conditions C1, C2, and C3.

## Usage

``` r
md_loglike_exp_series_C1_C2_C3(
  md,
  sysvar = "t",
  setvar = "x",
  deltavar = "delta"
)
```

## Arguments

- md:

  Masked data frame containing:

  - System lifetime column (default: "t")

  - Candidate set columns (default: "x1", "x2", ..., "xm")

  - Optional right-censoring indicator (default: "delta")

- sysvar:

  Column name for system lifetime. Default is "t".

- setvar:

  Column prefix for candidate set (Boolean matrix). Default is "x".

- deltavar:

  Column name for right-censoring indicator. Default is "delta". If NULL
  or column doesn't exist, assumes no censoring.

## Value

A function `f(theta)` that computes the log-likelihood given rate
parameters `theta = (lambda_1, ..., lambda_m)`.

## Details

For a series system with m components having exponential lifetimes with
rate parameters `theta = (lambda_1, ..., lambda_m)`:

- System reliability: `R(t) = exp(-sum(theta) * t)`

- System hazard: `h(t) = sum(theta)`

For observation i with system lifetime `s_i`, right-censoring indicator
`delta_i`, and candidate set `C_i`:

- If uncensored (delta_i = FALSE):
  `L_i(theta) = R(s_i) * sum_{j in C_i} h_j(s_i)`

- If censored (delta_i = TRUE): `L_i(theta) = R(s_i)`

The log-likelihood is: \$\$\ell(\theta) = \sum_i \log L_i(\theta)\$\$

For uncensored observations: \$\$\log L_i(\theta) = -s_i \sum\_{j=1}^m
\lambda_j + \log\left(\sum\_{j \in C_i} \lambda_j\right)\$\$

For censored observations: \$\$\log L_i(\theta) = -s_i \sum\_{j=1}^m
\lambda_j\$\$

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate some masked data
md <- tibble::tibble(
  t = c(0.5, 1.2, 0.8),
  x1 = c(TRUE, TRUE, FALSE),
  x2 = c(TRUE, FALSE, TRUE),
  x3 = c(FALSE, TRUE, TRUE),
  delta = c(FALSE, FALSE, TRUE)
)
ll <- md_loglike_exp_series_C1_C2_C3(md)
ll(c(1, 1.5, 2))  # Evaluate at theta = (1, 1.5, 2)
} # }
```
