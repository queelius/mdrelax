# informative_masking_by_rank

Returns `m=length(ts)` probabilities, where the j-th probability is the
probability that the j-th component will be in the candidate set.
probabilities are a function of rank (rather than by component failure
times).

## Usage

``` r
informative_masking_by_rank(ts, alpha, beta)
```

## Arguments

- ts:

  component failure times for the series system

- alpha:

  rate of change with respect to rank, as alpha -\> 0 probabilities go
  to discrete uniform distribution

- beta:

  max weight (for ranked 2 component)

## Details

The shape of the masking is given by two parametrs, `alpha`, which must
be non-negative, and `beta`, which must be between 0 and 1. As `alpha`
goes to 0, the probability vector approaches (beta, beta, ..., 1, beta,
..., beta) where the 1 is for the failed component and as `alpha` goes
to infinity, it assigns probability 1 and `beta` respectively to rank 1
(failed component) and rank 2 components (and the remaining
probabilities are zero).
