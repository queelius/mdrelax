# Methodology Auditor Report

**Date**: 2026-06-08
**Focus**: experimental design, statistical rigor, reproducibility. This pass re-runs the table-regeneration pipeline, independently reproduces the application headline, and audits the software-claim-vs-implementation gap for the paper's new robustness-interval tool.

## Reproducibility: re-verified by direct re-run

I ran `Rscript paper/regenerate_tables.R` against the committed `paper/data/*.rds` and diffed against the manuscript. **Tables 1, 2, and 3 reproduce exactly**, every cell:

- Table 1 (component fragility): C1 k_5 43%/0.46, lambda_1 37%/0.98, k_1 30%/0.21; C2 lambda_2 93%/0.98, k_3 37%/0.69, k_2 23%/0.44; C3 lambda_5 163%/0.96, lambda_3 146%/0.99, lambda_2 41%/1.00.
- Table 2 (n=500 slopes): C1 -0.32 (p<0.001, broken), C2 -0.04 (p=0.256, preserved), C3 +0.004 (p=0.689, preserved).
- Table 3 (n=2000 ablation): C1 -0.30 (p<0.001), C2 -0.01 (p=0.305), C3 -0.008 (p=0.204).
- Headline sys-hazard rel bias at max severity: -3.2% (C1), +2.6% (C2), -0.4% (C3). Reproduced.

This is a strong reproducibility posture: all three tables and the headline numbers come from one committed script over committed data. CRIT-1 from the 2026-05-27 review (Table 1 did not match the data) stays genuinely resolved.

## Design assessment (sound, unchanged)

5-component Weibull system, shapes (2.0,1.5,1.2,1.8,1.0), scales (3.0,4.0,5.0,3.5,4.5), k_5=1.0 giving an embedded exponential baseline. Three independent univariate sweeps (C1 11 levels, C2 11 levels, C3 9 levels), B=200 at n=500, B=100 at n=2000, all fits converged. Metrics bias / RMSE / 95% Wald coverage with delta-method SEs for the system hazard. Appropriate. The quadratic-linear-term slope test is a reasonable operationalization of "is first-order preservation broken?".

## NEW (this pass): software claim vs implementation gap for the robustness interval (MAJOR-leaning, see review.md calibration)

The paper's contribution 3 (robustness intervals) and contribution 7 (software) state the mdrelax package "provides ... robustness-interval construction" (introduction.tex item 7; conclusion item 5; appendix A4 / app:software). I checked the package against that claim:

- The implementation exists: `R/robustness_intervals.R` defines `ri_simulation` (Monte-Carlo robustness interval) and `ri_first_order` (closed-form ISNI-based interval, with the ISNI computed at line ~196). The roxygen header even gives the exact practitioner sentence the paper's application is missing: "my system-MTTF estimate is robust to C1 violations up to alpha = 0.27".
- BUT these functions are **not exported** (absent from NAMESPACE; `grep robustness|isni NAMESPACE` returns nothing) and have **no man pages** (`man/` has no ri_/robust/isni topic). They are de facto internal. A user installing mdrelax 1.0.0 cannot call the paper's headline new tool through the public API.
- The same tool is **never demonstrated on data** anywhere in the paper (carryover from 2026-06-04, methodology + novelty).

Net: the paper claims a software contribution that, as packaged, is unreachable, and the contribution is never exercised. Fix is small and high-value: (a) add `#' @export` + roxygen to `ri_first_order` (and `ri_simulation`), regenerate NAMESPACE/man; (b) add one robustness-interval row/sentence to the Guo application using `ri_first_order` on, e.g., system MTTF or a component scale. This single change simultaneously closes the "tool never demonstrated" gap (methodology), raises the apparent significance of contribution 3 (novelty), and makes the software claim true (format). Confidence: high (verified NAMESPACE and man/ directly).

## Application table: baseline scales do not match the package's own Guo MLE (MINOR, new this pass)

application.tex Table 4 (tab:guo-sensitivity) reports the s=0 standard fit as scales (1000, 897, 847). The package's committed canonical fit for the same dataset (`data-raw/guo_weibull_series_md.R`, object `guo_weibull_series_mle`) is (994.4, 908.9, 840.1), i.e. lambda_2 = 909 and lambda_3 = 840, not 897 and 847. The shapes (1.26, 1.17, 1.13) match the canonical (1.258, 1.164, 1.131) within rounding.

Important: this does NOT affect the paper's conclusion or its headline number. I computed h_T(249.5) from both scale sets and got 0.00307 either way (the system hazard is the robust estimand, so it is insensitive to the scale discrepancy). The issue is confined to the printed baseline component scales: a referee re-running the package will get lambda_2 = 909 / lambda_3 = 840, not the table's 897 / 847. The likely cause is that the application table was produced by a separate refit (the section says it refits under the relaxed-C2 model with known P(s); at s=0 that should coincide with the standard fit but may have used different start values / tolerances) rather than from `guo_weibull_series_mle`. Fix: regenerate the application table from a single committed code path tied to the canonical MLE (mirroring what regenerate_tables.R does for the simulation tables), and reconcile the s=0 row. Confidence: high on the discrepancy; medium on the cause.

## Carryover methodology items (still open)

### MAJOR (carryover, partially mitigated). C2/C3 n=2000 ablation underpowered as a positive preservation claim.
Table 3 C3 n=2000: linear -0.008, p=0.204, B=100; C2 n=2000: -0.01, p=0.305, B=100. A fail-to-reject at B=100 with a slope of 0.008 is a Type II non-detection, not a positive demonstration. Mitigation: the n=500 B=200 runs also fail to reject, the C1 case (which SHOULD reject) rejects strongly at both sizes (p<0.001), so the test demonstrably has power where a real effect exists; and the structural argument (exponential exact preservation + Weibull first-order via the same marginal-distribution mechanism) is the real basis. Fix: one sentence in Section 4.6 stating the preservation verdict rests on the structural argument and is corroborated (not established) by the non-significant slopes, with the C1 rejection at the same B as the power check. Optionally rerun C2/C3 n=2000 at B=200 (cheap; everything converged).

### MINOR (carryover). Application n=30 "2% range" is at the noise floor.
application.tex reports the Guo system hazard "varying by less than 2% across the full severity range" at n=30, m=3. At n=30 a 2% movement is within Monte Carlo / optimizer noise. The qualitative point (system hazard stable, component scales move) is correct; present "2%" as illustrative, not as evidence.

### MINOR (carryover). Wald coverage under misspecification.
Coverage uses model-based Fisher-information Wald intervals (assume correct specification). The paper acknowledges this in Limitations and names the sandwich estimator as the right tool, deferred to software. Acceptable disclosure; the coverage-collapse phenomenon is reported as a feature, not hidden.

## Reproducibility scorecard

- Tables 1/2/3: regenerable from committed data via one committed script. Strong.
- Application Table 4: NOT regenerated from a committed path; baseline scales drift from the package's canonical MLE (see above). The one reproducibility soft spot.
- Robustness-interval tool: implemented but unexported, undocumented, and never demonstrated. The main methodological gap.

## Bottom line (methodology)

Simulation reproducibility is strong and re-verified. The two items worth fixing before submission are: (1) export + demonstrate the robustness interval (closes a software-claim gap AND the never-demonstrated gap in one move); (2) reconcile the Guo application baseline scales to the package's canonical MLE via a committed code path. Plus the carryover one-sentence C2/C3 power caveat. None is a correctness blocker.
