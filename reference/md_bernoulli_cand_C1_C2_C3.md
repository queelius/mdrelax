# md_bernoulli_cand_C1_C2_C3

Bernoulli candidate set model is a particular type of *uninformed*
model. Note that we do not generate candidate sets with this function.
See `md_cand_sampler` for that.

## Usage

``` r
md_bernoulli_cand_C1_C2_C3(
  md,
  p,
  compvar = "t",
  qvar = "q",
  deltavar = "delta"
)
```

## Arguments

- md:

  masked data.

- p:

  a vector of probabilities

- compvar:

  column name of the component lifetime variables, defaults to `t`,
  e.g., `t1, t2, ..., tm`.

- qvar:

  column prefix for component probabilities, defaults to `q`, e.g.,
  `q1, q2, ..., qm`.

- deltavar:

  column name of the right-censoring indicator variable, defaults to
  `delta`.

## Details

This model satisfies conditions C1, C2, and C3. The failed component
will be in the corresponding candidate set with probability 1, and the
remaining components will be in the candidate set with probability `p`.
