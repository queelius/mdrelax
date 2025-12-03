# Score function (gradient of log-likelihood) for exponential series system (C1,C2,C3)

Returns the score function (gradient of log-likelihood with respect to
theta) for an exponential series system with masked data under C1, C2,
C3 conditions.

## Usage

``` r
md_score_exp_series_C1_C2_C3(
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

A function `g(theta)` that computes the score (gradient) vector at
parameter values `theta = (lambda_1, ..., lambda_m)`.

## Details

The score for component j is: \$\$\frac{\partial \ell}{\partial
\lambda_j} = -\sum_i s_i + \sum\_{i: \delta_i=0} \frac{I(j \in
C_i)}{\sum\_{k \in C_i} \lambda_k}\$\$

where `I(j in C_i)` is 1 if component j is in candidate set C_i, 0
otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
md <- tibble::tibble(
  t = c(0.5, 1.2, 0.8),
  x1 = c(TRUE, TRUE, FALSE),
  x2 = c(TRUE, FALSE, TRUE),
  x3 = c(FALSE, TRUE, TRUE)
)
score <- md_score_exp_series_C1_C2_C3(md)
score(c(1, 1.5, 2))
} # }
```
