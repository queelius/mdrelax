# Methodology Auditor Report

**Date**: 2026-06-04
**Focus**: experimental design, statistical rigor, reproducibility. Priority: re-verify Table 1 against the data (the prior Critical) and re-reproduce Tables 2 and 3.

## Verdict on CRIT-1 (Table 1 did not match the data): FIXED, confirmed by direct re-run

I ran the new `paper/regenerate_tables.R` against `paper/data/*.rds` (invoked from the repo root so the script's `normalizePath("paper")` fallback resolves) and diffed every cell against the current `simulations.tex` Table 1 (tab:component-fragility, lines 120-136).

**Every entry matches exactly.**

| Violation | Param | Manuscript | Regenerated | Severity (ms / regen) | Match |
|---|---|---|---|---|---|
| C1 | k_5 | 43%, cov 0.46 | 43%, 0.46 | a=0.50 / 0.50 | yes |
| C1 | lambda_1 | 37%, cov 0.98 | 37%, 0.98 | a=0.50 / 0.50 | yes |
| C1 | k_1 | 30%, cov 0.21 | 30%, 0.21 | a=0.50 / 0.50 | yes |
| C2 | lambda_2 | 93%, cov 0.98 | 93%, 0.98 | s=0.5 / 0.5 | yes |
| C2 | k_3 | 37%, cov 0.69 | 37%, 0.69 | s=1.0 / 1.0 | yes |
| C2 | k_2 | 23%, cov 0.44 | 23%, 0.44 | s=1.0 / 1.0 | yes |
| C3 | lambda_5 | 163%, cov 0.96 | 163%, 0.96 | a=1.25 / 1.25 | yes |
| C3 | lambda_3 | 146%, cov 0.99 | 146%, 0.99 | a=1.75 / 1.75 | yes |
| C3 | lambda_2 | 41%, cov 1.00 | 41%, 1.00 | a=2.00 / 2.00 | yes |

Note the table membership also changed correctly relative to the buggy version: the worst-three C2 components are now lambda_2 / k_3 / k_2 (regenerated ranking), not the old lambda_2 / k_3 / lambda_3; and C3 is lambda_5 / lambda_3 / lambda_2. The new caption explicitly says "Numbers are regenerated directly from the sweep RDS files via paper/regenerate_tables.R," which is exactly the reproducibility hook a Technometrics referee wants. **The prior Critical is genuinely resolved, not patched over.**

## Tables 2 and 3 re-reproduced

`quadratic_linear_term()` (regress bias ~ sev + I(sev^2) on the sys_hazard subset) reproduces the manuscript verbatim:
- Table 2 (n=500): C1 -0.32 (p<0.001, broken); C2 -0.04 (p=0.256, preserved); C3 +0.004 (p=0.689, preserved).
- Table 3 (n=2000): C1 -0.30 (p<0.001); C2 -0.01 (p=0.305); C3 -0.008 (p=0.204).

All match. The "theory match" column is internally consistent (C1 -> cor:hierarchy item 3 broken; C2 -> prop:hazard-robust; C3 -> cor:hierarchy item 3 preserved).

## Headline numeric claims re-verified against data

- simulations.tex lines 91-92: system-hazard relative bias at max severity = -3.2% (C1), +2.6% (C2), -0.4% (C3). I recomputed directly from the rds files: -3.2%, +2.6%, -0.4%. Exact match.
- simulations.tex line 311 ablation sentence: "lambda_5 under C3: +163% at n=500 vs +149% at n=2000." Recomputed: n=500 max = 162.6% -> 163%; n=2000 max = 149.4% -> 149% (at severity 2.0). Exact match.

The numerical claims in the paper are now fully reproducible from the committed data via a single committed script. This is a strong reproducibility posture.

## Design assessment (unchanged from prior review, still sound)

- 5-component Weibull system, shapes (2.0,1.5,1.2,1.8,1.0), scales (3.0,4.0,5.0,3.5,4.5), with k_5=1.0 giving an embedded exponential baseline. Sensible.
- Three independent univariate sweeps (C1: 11 levels; C2: 11 levels; C3: 9 levels), B=200 at n=500, B=100 at n=2000. All fits converged.
- Metrics: bias, RMSE, 95% Wald coverage, with delta-method SEs for the system hazard. Appropriate.
- The slope-test battery (full-range, small-alpha, quadratic-linear-term) is a reasonable operationalization of "is first-order preservation broken?" The quadratic-linear-term is the cleanest of the three and is what the tables report.

## Open methodology issues

### MAJ (carryover, partially mitigated). C3 n=2000 ablation is underpowered as a positive preservation claim.
- Table 3 C3 at n=2000: linear coef -0.008, p=0.204, B=100. A fail-to-reject with B=100 and a slope of 0.008 is a Type II non-detection, not a positive demonstration of preservation. Same concern applies to C2 n=2000 (-0.01, p=0.305).
- Mitigation present: the n=500 run (B=200) for C2/C3 also fails to reject, and the n=500 C1 result (the one that SHOULD reject) does reject strongly at both sample sizes (-0.32 and -0.30, p<0.001). So the test demonstrably has power where a real effect exists, which substantially blunts the Type II worry. The structural argument (exponential exact preservation + Weibull first-order via the same marginal-distribution mechanism) is the real basis for the preservation claim; the simulation is corroborative.
- Recommendation: add one sentence to Section 4.6 stating that the C2/C3 "preserved" verdict rests on the structural first-order argument and is corroborated (not established) by the non-significant slopes, and that the C1 rejection at the same B confirms the test is powered to detect a real first-order effect. Optionally rerun C2/C3 n=2000 at B=200 to remove the asymmetry with the main run; cheap given everything converged.

### MIN. Application n=30 is at the noise floor for the "2% range" claim.
- application.tex line 51 reports the Guo system hazard "varying by less than 2% across the full severity range" at n=30, m=3. At n=30 a 2% movement is within Monte Carlo / optimizer noise. The qualitative point (system hazard stable, component scales move) is correct and consistent with the simulation, but the "2%" should be presented as illustrative, not as evidence. This section is also still thin (MIN-10 carryover): adding a robustness-interval row (now that Definition 3.20 exists) would convert the application from a stability anecdote into a demonstration of the paper's own new tool. Strongly recommended, since the new robustness-interval machinery is never actually exercised on data anywhere in the paper.

### MIN. Wald coverage under misspecification.
- Coverage uses model-based Fisher-information Wald intervals, which assume correct specification. The paper acknowledges this in Limitations (discussion.tex lines 143-150) and notes the sandwich estimator is the right tool, deferred to software. Acceptable disclosure; the coverage-collapse phenomenon (lambda_1 under C3 -> 0) is reported as a feature (confidently-wrong inference) rather than hidden.

### MINOR (carryover). State file drift.
- The 2026-05-27 review flagged that .papermill/state.md called the C1/C3 sweeps "planned"/"script-ready." The state.md I was given in-context now reflects the unified framework and refined thesis; if the experiments block still lists planned/script-ready statuses it should be updated to complete (B=200 main, B=100 ablation, dated). This is a hygiene item, not a paper defect.

## Reproducibility scorecard

- Tables 1, 2, 3: regenerable from committed data via one committed script. Strong.
- Figures: produced by run_sensitivity_sweep*.R; figure path issue is a format concern (see format-validator), not a reproducibility one.
- The robustness-interval tool (the paper's new contribution) is defined and implemented but not demonstrated on data. This is the main methodological gap remaining: the paper proves the interval and ships the code but never shows one.

## Bottom line (methodology)

CRIT-1 is fixed and verified by direct re-run; all table and headline numbers are reproducible. Design is sound. Remaining items are: (a) the C2/C3 n=2000 power caveat (one sentence + optional B=200 rerun), and (b) actually exercise the new robustness interval on the Guo data. Neither is a blocker; both would materially strengthen a Technometrics submission.
