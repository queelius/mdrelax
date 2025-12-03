# Component failure probability under conditions C1 and C2

Computes the probability that component j caused system failure given
the candidate set C and system failure time t.

## Usage

``` r
md_series_component_failure_probability_C1_C2(theta, C, j = NULL)
```

## Arguments

- theta:

  Rate parameters for components.

- C:

  Candidate set (Boolean vector or matrix).

- j:

  Component index to compute probability for.

## Value

Probability that component j caused the failure.

## Details

Under conditions C1 and C2, this probability is: \$\$P(K=j \| C, T=t) =
\frac{\lambda_j}{\sum\_{k \in C} \lambda_k}\$\$
