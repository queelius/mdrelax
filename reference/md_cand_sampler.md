# md_cand_sampler

Candidate set generator. Requires columns for component probabilities
e.g., `q1,...,qm` where `qj` is the probability that the jth component
will be in the corresponding candidate set generated for that
observation in the `md` table.

## Usage

``` r
md_cand_sampler(md, qvar = "q", setvar = "x")
```

## Arguments

- md:

  masked data.

- qvar:

  column prefix for component probabilities, defaults to `q`, e.g.,
  `q1, q2, ..., qm`.

- setvar:

  column prefix for candidate sets (as Boolean matrix), defaults to `x`,
  e.g., `x1, x2, ..., xm`.
