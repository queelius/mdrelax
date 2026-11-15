# Maintenance and verified scope

The package remains available as research software. The bounded November 2026
closeout supports the
[exponential total-rate worked result](paper/reproduction/exponential-total-rate/README.md),
its reproduction script and regression tests. Compatibility and correctness
reports affecting that path take priority. This states a maintenance scope,
not a response-time or indefinite support commitment.

The existing Weibull and relaxed-model implementations remain available, but
the closeout does not establish all of their theoretical claims or reproduce
the full paper's simulations. In particular, `ri_first_order()` and
`ri_simulation()` are exported experimental sensitivity tools. A first-order
extrapolation is not a certified interval of global robustness or a confidence
interval. Do not infer such a guarantee from the Guo application narrative.

No new features, application platform, or CRAN submission is planned as part of
this closeout. The other ecosystem packages retain their own maintenance and
distribution responsibilities. No repository has been archived or withdrawn.

Historical review notes are dated records, not current editorial decisions.
An AI review label such as “minor revision” does not establish journal acceptance
or human peer review. Report new issues through the owning repository.
