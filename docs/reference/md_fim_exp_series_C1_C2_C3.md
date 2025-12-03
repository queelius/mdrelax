# Observed Fisher information matrix for exponential series system (C1,C2,C3)

Returns a function that computes the observed Fisher information matrix
(FIM) for an exponential series system with masked data under C1, C2, C3
conditions.

## Usage

``` r
md_fim_exp_series_C1_C2_C3(md, sysvar = "t", setvar = "x", deltavar = "delta")
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

A function `I(theta)` that computes the m x m observed Fisher
information matrix at parameter values `theta`.

## Details

The observed FIM is the negative Hessian of the log-likelihood, which
under regularity conditions converges to the expected FIM and can be
inverted to obtain the asymptotic variance-covariance matrix of the MLE.

The (j,k) element of the observed FIM is: \$\$I\_{jk}(\theta) =
\sum\_{i: \delta_i=0} \frac{I(j \in C_i) I(k \in C_i)}{(\sum\_{l \in
C_i} \lambda_l)^2}\$\$

Note: The survival contribution to the Hessian is zero (second
derivative of linear term is zero), so only uncensored observations
contribute.

## Examples

``` r
if (FALSE) { # \dontrun{
md <- tibble::tibble(
  t = c(0.5, 1.2, 0.8),
  x1 = c(TRUE, TRUE, FALSE),
  x2 = c(TRUE, FALSE, TRUE),
  x3 = c(FALSE, TRUE, TRUE)
)
fim <- md_fim_exp_series_C1_C2_C3(md)
fim(c(1, 1.5, 2))
} # }
```
