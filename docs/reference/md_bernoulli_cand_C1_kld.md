# md_bernoulli_cand_C1_kld

For each row (observation) in `md`, a probability
`Q = (q1, q2, ..., qm)` is constructed such that `qj` is the probability
that the j-th component is in the candidate set, `qk = 1`, where `k` is
failed component.

## Usage

``` r
md_bernoulli_cand_C1_kld(
  md,
  p,
  d,
  eps = 1e-04,
  max_iter = 100000L,
  lr = 1,
  lambda = 1,
  alpha0 = 5,
  beta0 = 0.5,
  debug = F
)
```

## Arguments

- md:

  component failure times for the series system

- p:

  numeric, defines `P = (p, ..., p, 1, p, ..., p)`.

- d:

  numeric, the KL divergence from P = (p, p, ..., p, 1, p, ..., p) to
  try to obtain

- eps:

  numeric, stopping condition.

- max_iter:

  Integer, maximum number of iterations before giving up.

- lr:

  numeric, learning rate.

- lambda:

  numeric, controls how much the two constraints are weighted. Lower
  value specifies more enforcement of the KL-divergence constraint being
  closer to `d`. Defaults to 1.

- alpha0:

  numeric, initial guess for `alpha` parameter of
  `informative_masking_by_rank`.

- beta0:

  numeric, initial guess for `beta` parameter of
  `informative_masking_by_rank`.

- debug:

  Logical, whether to output debugging information while running

## Details

`Q` is an *informed* candidate model that uses
`informative_masking_by_rank` to assign higher probabilities to
components that failed earlier (which is something we typically only
know in, say, a simulation study).

The probabilities `Q` have two constraints on them. Let
`P = (p, ..., p, 1, p, ..., p)` be the bernoulli candidate model that
satisfies conditions C1, C2, and C3. Then, the KL-divergence between `P`
and `Q` is as close as possible to `d` while satisfying
`sum(P) == sum(Q)`.

For `d = 0`, `Q == P`. As `d` increases, `Q` becomes more informative
about the components. Given the structure of
`informative_masking_by_rank`, it may not be possible to satisfy every
`d` specified, but we get as close as we can, which should permit useful
experiments.
