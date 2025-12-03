# kl_divergence_bernoulli

KL divergence between two Bernoulli random vectors. p = (p1,p2,...,pn)
and q = (q1,q2,...,qn).

## Usage

``` r
kl_divergence_bernoulli(p, q)
```

## Arguments

- p:

  probability of success for the first distribution

- q:

  probability of success for the second distribution

## Value

KL divergence between the two distributions

## Examples

``` r
kl_divergence_bernoulli(c(0.1,0.1), c(0.1,0.1)) = 0
#> Error in kl_divergence_bernoulli(c(0.1, 0.1), c(0.1, 0.1)) = 0: target of assignment expands to non-language object
```
