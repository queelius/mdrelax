# Exponential Series Systems with Masked Data

## Introduction

The `mdrelax` package provides tools for likelihood-based inference in
series systems with masked failure data. This vignette demonstrates the
basic workflow using exponential component lifetimes.

## Setup

``` r
library(mdrelax)
```

## Simulating Masked Data

We simulate a series system with 3 components having exponential
lifetimes with rates $\lambda = (1.0,1.5,2.0)$.

``` r
set.seed(42)

# True parameters
lambda_true <- c(1.0, 1.5, 2.0)
n <- 100  # sample size
m <- 3    # number of components

# Generate component failure times
T_mat <- matrix(rexp(n * m, rate = rep(lambda_true, each = n)),
                nrow = n, ncol = m)

# System failure time is minimum of component times
t_sys <- apply(T_mat, 1, min)

# Identify which component failed (latent)
k_failed <- apply(T_mat, 1, which.min)

# Create data frame with component times
md <- data.frame(
  t = t_sys,
  t1 = T_mat[, 1],
  t2 = T_mat[, 2],
  t3 = T_mat[, 3],
  delta = rep(1, n)  # all observed (no censoring): delta=1 means failure observed
)
```

## Applying Bernoulli Masking (C1-C2-C3)

Under the standard C1-C2-C3 conditions, we apply Bernoulli masking where
each non-failed component is included in the candidate set with
probability $p$.

``` r
# Apply Bernoulli masking with p = 0.3
md <- md_bernoulli_cand_C1_C2_C3(md, p = 0.3)

# Sample candidate sets
md <- md_cand_sampler(md)

# View first few rows
head(md[, c("t", "x1", "x2", "x3")])
#>             t    x1    x2    x3
#> 1 0.198336812  TRUE FALSE FALSE
#> 2 0.064822223  TRUE FALSE  TRUE
#> 3 0.124212025  TRUE FALSE  TRUE
#> 4 0.038191898  TRUE FALSE  TRUE
#> 5 0.250735989 FALSE FALSE  TRUE
#> 6 0.000824189  TRUE  TRUE  TRUE
```

The columns `x1`, `x2`, `x3` indicate which components are in each
candidate set (1 = in candidate set, 0 = not in candidate set).

## Maximum Likelihood Estimation

We compute the MLE for the component rate parameters:

``` r
fit <- md_mle_exp_series_C1_C2_C3(md)
print(fit)
#> Maximum likelihood estimator of type mle_exp_series_C1_C2_C3 is normally distributed.
#> The estimates of the parameters are given by:
#> [1] 1.226501 1.075383 1.525874
#> The standard error is  0.2588029 0.2361879 0.2839848 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>             2.5%    97.5%
#> param1 0.7192569 1.733746
#> param2 0.6124630 1.538303
#> param3 0.9692738 2.082474
#> The MSE of the individual components in a multivariate estimator is:
#>              [,1]         [,2]         [,3]
#> [1,]  0.066978964 -0.006206286 -0.013825185
#> [2,] -0.006206286  0.055784723 -0.008415402
#> [3,] -0.013825185 -0.008415402  0.080647362
#> The log-likelihood is  -39.0219 .
#> The AIC is  84.0438 .
```

## Comparing to True Values

``` r
# Extract point estimates using algebraic.mle accessor
estimates <- algebraic.mle::params(fit)

# Compare to true values
comparison <- data.frame(
  Component = paste0("lambda_", 1:m),
  True = lambda_true,
  Estimate = estimates,
  Bias = estimates - lambda_true
)
print(comparison)
#>   Component True Estimate       Bias
#> 1  lambda_1  1.0 1.226501  0.2265013
#> 2  lambda_2  1.5 1.075383 -0.4246173
#> 3  lambda_3  2.0 1.525874 -0.4741262
```

## Fisher Information Matrix

The Fisher Information Matrix provides information about estimation
precision:

``` r
# Create FIM function and evaluate at MLE
fim_fn <- md_fim_exp_series_C1_C2_C3(md)
fim <- fim_fn(estimates)
print(fim)
#>           [,1]      [,2]      [,3]
#> [1,] 15.737278  2.192326  2.926569
#> [2,]  2.192326 18.518158  2.308160
#> [3,]  2.926569  2.308160 13.142209

# Standard errors from inverse FIM
se <- sqrt(diag(solve(fim)))
cat("\nStandard errors:", round(se, 4), "\n")
#> 
#> Standard errors: 0.2588 0.2362 0.284
```

## Summary

This vignette demonstrated:

1.  Simulating exponential series system data
2.  Applying Bernoulli masking under C1-C2-C3 conditions
3.  Computing the MLE for component rate parameters
4.  Evaluating estimation precision via the Fisher Information Matrix

For more advanced usage including informative masking models and relaxed
conditions, see the package documentation.
