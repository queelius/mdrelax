# What masked failures can still tell you

**An exponential system's fitted total failure rate can remain identifiable
even when its component rates are not.** This worked result explains and checks
the factorization already given in Theorem 3.6 of the
[archived mdrelax paper, v0.2.0](https://doi.org/10.5281/zenodo.20468529).
It is a reproducibility companion to an existing result, not a new theorem,
general robustness guarantee, or complete verification of that paper.

## The question

A series system fails when any component fails. A diagnostic report may give
only a set of possible failed components. Can the system's failure rate still
be estimated when those sets are unreliable?

For the **standard fitted exponential likelihood**, yes in a precise sense:
changing nonempty candidate sets while holding the observed times and failure
indicators fixed does not change its unconstrained fitted total rate. It can
change the component rates, or leave them unidentified altogether.

## Claim and short argument

Let `t[i]` be retained observation times (including right-censored exposure),
`delta[i]` indicate an observed failure, and `C[i]` be a nonempty candidate set
for each failure. Set `d = sum(delta)` and `E = sum(t)`, with `d > 0`, `E > 0`.
For component rates `theta[j]`, write `S = sum(theta)` and `phi = theta / S`.
The standard fitted log-likelihood, up to terms independent of the rates, is

```text
ell(S, phi) = d log(S) - S E + sum_{i: delta[i]=1} log(sum_{j in C[i]} phi[j]).
```

The last term contains the candidate sets but not `S`. Thus the maximizer of
the total-rate part is `S_hat = d/E`, with second derivative `-d/S^2 < 0`.
The statement also applies when component proportions have only a boundary
supremum: every fixed admissible composition has the same maximizing total.
Separate user-imposed bounds on component rates can destroy this separation.

This is a **finite-sample fitted-likelihood identity**. Consistency for a true
population rate additionally requires that retention and censoring preserve
the appropriate exponential rate/exposure relation. Independent administrative
right censoring is one such setting. It does not imply finite-sample
unbiasedness: for uncensored independent exponential observations, `n/sum(t)`
has expectation `n*S/(n-1)` when `n > 1`.

## Reproduce from a clean installed package

From the root of this source checkout, with R installed:

```sh
sh paper/reproduction/exponential-total-rate/reproduce.sh /tmp/mdrelax-results
```

The script installs this checkout into a fresh temporary R library and runs R
with `--vanilla`. The package and verification script require only base R and
the recommended `stats` package. The optional output directory receives
structured CSV results. The temporary installed library is removed afterward.
Use the exact source commit linked by the
[central research record](https://github.com/queelius/research-records/tree/main/records/masked-exponential-total-rate)
to reproduce this verification rather than an evolving branch tip.

The fixed synthetic example has ten failures, two right-censored observations,
and total exposure 7.8. It produces, to numerical tolerance:

| Reports supplied to the fit | Component 1 | Component 2 | Component 3 | Total rate |
|---|---:|---:|---:|---:|
| Singleton sets, counts 5, 3, 2 | 0.641026 | 0.384615 | 0.256410 | 1.282051 |
| Labels 1 and 3 exchanged | 0.256410 | 0.384615 | 0.641026 | 1.282051 |
| Every component in every failure set | arbitrary allocation | arbitrary allocation | arbitrary allocation | 1.282051 |

The equal allocation returned by the optimizer in the fully masked case is
an initialization-dependent choice, not evidence for equal component rates.
The observed information has rank one in that case.

The script checks likelihood factorization, radial score, an independent scalar
optimizer, fitted rates and identification. Package regression tests also check
radial observed information. Absolute tolerances are `1e-10` for algebraic
script checks and `2e-5` for numerical fitted rates.

## Controls that prevent overinterpreting the result

- **Censoring:** throwing away censored exposure produces a different rate.
- **Selection:** a deterministic grid of exponential quantiles at rate 3 gives
  an overall fitted rate near 3.001. Keeping only times above the median gives
  about 1.773. The fitted identity still holds, but now estimates the selected
  sample's rate functional; its population limit is `3/(1 + log(2))`, not 3.
- **No observed failures:** the rate likelihood has its supremum at zero. An
  optimizer with positive lower bounds cannot return that boundary exactly.
- **Binding component bounds:** imposing small upper bounds prevents the
  optimizer from reaching `d/E`.
- **Component attribution:** reversing diagnostic labels reverses component
  rankings in the example despite leaving the total unchanged.

The synthetic controls establish implementation behavior, not field performance.
They do not justify extending the constant-hazard argument to general Weibull
models, inferring a causal repair benefit, or treating unobserved AI answer
correctness as a censored completion time.

## Evidence, citation and maintenance

[evidence.json](evidence.json) records the tested source hashes, environment,
test/coverage results, and the verified historical archive. The adjacent logs
record the reproduction commands' output. The original PDF was downloaded from
Zenodo and checked against its advertised checksum; its SHA-256 is also recorded.
The November 2026 checks are new supplementary evidence and are **not contained
in the June 2026 archived PDF**. Cite the exact paper DOI for the prior result
and the exact source commit for this companion.

Author: **Alexander Richard Towell**, ORCID
[0000-0001-6443-9897](https://orcid.org/0000-0001-6443-9897).
The exponential rate estimator and likelihood factorization use established
likelihood methods. The source paper credits the preceding masked-failure and
misspecified-likelihood literature; this companion makes no historical novelty
claim. Substantial OpenAI Codex assistance was used in the implementation,
verification and documentation; AI review is not human peer review or formal
proof verification. Code and this companion are GPL-3.0-or-later; the existing
paper archive has its separately recorded CC-BY-4.0 license.

See [maintenance scope](../../../MAINTENANCE.md). Stable central record ID:
`masked-exponential-total-rate`. The supported completion is this worked path;
the broader manuscript and every package in the ecosystem have not been
revalidated by this exercise.
