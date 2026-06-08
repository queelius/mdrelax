# Unified editorial review: mdrelax (2026-06-04)

Final pre-submission pass. Companion sensitivity paper to the masked-causes
foundational work. Primary target: Technometrics (backups IEEE TR, Lifetime
Data Analysis, JSS/R Journal).

This synthesis integrates seven specialist lenses (logic-checker,
methodology-auditor, novelty-assessor, citation-verifier, prose-auditor,
format-validator, literature-context). The orchestrator that launched them was
interrupted by a rate-limit before writing this synthesis; the per-specialist
files are complete and were written 2026-06-04. WebSearch proper was
unavailable to the specialists; prior-art freshness was checked via Crossref
and arXiv REST APIs plus DOI resolution.

## Verdict: minor-revision (verging on ready)

Both Criticals from the 2026-05-27 review are confirmed FIXED by direct
verification, not patched over. No new Critical. One carried Major
(under-powered C3 ablation, mitigated). The remaining items are minor and
submission-mechanical.

## The two prior Criticals: confirmed resolved

- **CRIT-1 (Table 1 did not match the data): FIXED, confirmed by re-run.** The
  methodology-auditor ran `paper/regenerate_tables.R` against `paper/data/*.rds`
  and diffed every cell against the current Table 1; they match. The worst-three
  C2 components changed correctly relative to the buggy version (now
  lambda_2 / k_3 / k_2), and the caption documents the regeneration hook, which
  is exactly the reproducibility signal a Technometrics referee wants. Tables 2
  and 3 re-reproduced and still correct.
- **CRIT-2 (robustness intervals named but never defined): FIXED.** Definition
  3.19 (local sensitivity index), Definition 3.20 (robustness interval), and
  Proposition 3.21 (coverage/remainder) are now in Section 3.8; the
  supplementary-material forward reference is gone. The logic-checker verified
  Proposition 3.21 is sound and submission-grade: the Lagrange-remainder
  second-derivative bound is computed correctly by the chain rule and captures
  both curvature sources, and the C^2 / bounded-gradient / bounded-Hessian
  hypotheses are sufficient for the stated deterministic pseudo-true-bias bound.
  Remark 3.23 honestly states the interval bounds asymptotic bias and does not
  absorb Monte Carlo uncertainty.

## Proof / numerical soundness

Submission-grade. The C2-sensitivity core (exact-preservation theorem for
exponential components, Cor 3.5; first-order Weibull robustness, Cor 3.11 plus
the non-significant Table 2 slope) is sound and consistent with the masked-causes
C2 definition. The new robustness-interval proposition is correct as stated.

## Major (1, carried and partially mitigated)

- **C3 n=2000 ablation is under-powered as a positive preservation claim**
  (B=100 for a slope magnitude ~0.008). The "preserved" verdict rests on the
  structural first-order argument and is corroborated, not established, by the
  non-significant slopes. Fix: one sentence in Section 4.6 stating this
  explicitly (and noting the C1 rejection at the same B confirms the test is
  powered to detect a real effect); optionally rerun C2/C3 n=2000 at B=200 to
  remove the asymmetry with the main run (cheap, everything converged).

## Minor (submission-mechanical; the first is the one to actually fix)

1. **Title vs metadata mismatch** (prose-auditor + format-validator). The
   on-page title (main.tex) uses the unified "...Coarsening Conditions on Masked
   Failure Data" framing; refs.bib, state.md, and the published Zenodo record
   say "...Non-Informative Masking Assumption". The PDF title and the DOI
   landing-page title will differ at submission. Pick one canonical title and
   reconcile (the Zenodo metadata is editable on the record; or align the paper
   to the deposited title). This is the one pre-submission consistency item.
2. **Robustness interval never demonstrated on data** (methodology-auditor +
   novelty-assessor). Contribution 3 is now defined and proved but never
   exercised; adding one robustness-interval row to the Guo application
   (currently a stability anecdote at n=30) would convert it into a
   demonstration of the paper's own new tool. Strongly recommended; it raises
   the apparent significance of the headline new contribution.
3. **One-sentence boundary vs the foundational companion** (novelty-assessor).
   Add after the contribution list: the foundational paper establishes when the
   C1-C2-C3 likelihood is valid; this paper quantifies what is lost when those
   conditions fail. Closes the only novelty-framing gap and helps the sibling
   papers that cite this one.
4. **Soften "multi-dimensional sensitivity geometry"** (contribution 1): the
   sweeps are three independent univariate axes, not a joint geometry. Rephrase
   to "calibrated severities along three independent violation axes," or add a
   genuine 2-D joint sweep.
5. **Citations**: cite the dead `Stefanski2002` entry in the sandwich-variance
   Limitations sentence (or remove it); add Hu-Huang-Shen 2023 (the most recent
   reliability-side masked-cause paper, closes a ~20-year gap in the reliability
   prior-work paragraph). Optional: Lo et al. 2022.
6. **Format**: co-locate figures into `paper/figures/` (drop the `../inst` path)
   and remove orphan section files for a clean submission tarball; drop the
   "1276 unit tests" vanity count in appendix.tex; retemplate to the ASA class
   and trim to the Technometrics page limit (move Software/Score-Functions
   appendices to supplementary; the 41-page double-spaced draft is over typical
   regular-article length).

## Prior-art freshness

Current and defensible. The concurrent Bakoyannis 2025 robustness-interval paper
(the main "are they scooped?" risk) is cited, correctly characterized as
concurrent, and differentiated (semiparametric Cox binary-missing vs parametric
set-valued masking). No 2024-2026 work supersedes the contribution. One
reasonable currency add (Hu-Huang-Shen 2023) noted above.

## Single most important remaining item

Reconcile the title-vs-Zenodo-metadata mismatch before submission (the PDF and
DOI landing-page titles currently differ). Everything else is optional polish
or a one-sentence honesty caveat.
