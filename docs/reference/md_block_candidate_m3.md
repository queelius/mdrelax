# md_block_candidate_m3

Block candidate model. It produces an MLE that is non-unique, i.e., if
estimating exponential series system, the MLE for λ₁, λ₂, and λ₃
satisfies `λ̂₁ + λ̂₂ = λ₁ + λ₂` and `λ̂₃ = λ₃`.

## Usage

``` r
md_block_candidate_m3(md, qvar = "q", compvar = "t")
```

## Arguments

- md:

  masked data

- qvar:

  column prefix for component probabilities, defaults to `q`, e.g.,
  `q1, q2, ..., qm`.

- compvar:

  Prefix-code for the component lifetimes, defaults to `t`, e.g.,
  `t1, t2, ..., tm`.

## Details

This illustrates one of the conditions for the MLE to be unique:

    (1) multiple components must not always show together. if they
    do, then the best we can learn is that they satisfy some hyper-plane
    constraint.
