# generate_relaxed_cand_C1

Generates a probability `Q = (q1, q2, ..., qm)` such that `qj` is the
probability that the j-th component is in the candidate set, `qk = 1`,
where `k` is failed component. `Q` is an *informed* candidate model that
uses `informative_masking_by_rank` to assign higher probabilities to
components that failed earlier (which is something we typically only
know in, say, a simulation study).

## Usage

``` r
generate_relaxed_cand_C1(
  d,
  ts,
  p,
  debug = F,
  eps = 1e-08,
  alpha0 = 1,
  beta0 = p,
  lambda = 1,
  max_iter = 10000L,
  lr = 1
)
```

## Arguments

- d:

  numeric, the KL divergence from P = (p, p, ..., p, 1, p, ..., p) to
  try to obtain

- ts:

  component failure times for the series system

- p:

  numeric, defines `P = (p, ..., p, 1, p, ..., p)`.

- debug:

  Logical, whether to output debugging information while running

- eps:

  numeric, stopping condition.

- alpha0:

  numeric, initial guess for `alpha` parameter of
  `informative_masking_by_rank`.

- beta0:

  numeric, initial guess for `beta` parameter of
  `informative_masking_by_rank`.

- lambda:

  numeric, controls how much the two constraints are weighted. Lower
  value specifies more enforcement of the KL-divergence constraint being
  closer to `d`. Defaults to 1.

- max_iter:

  Integer, maximum number of iterations before giving up.

- lr:

  numeric, learning rate.

## Details

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
