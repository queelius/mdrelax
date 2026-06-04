# Methodology Auditor Report

**Date**: 2026-05-27
**Paper**: Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption

## Scope

Verifies (a) state-file vs manuscript simulation status, (b) reproducibility of Tables 1, 2, 3 against actual data files in paper/data/, (c) ablation power, (d) Guo application adequacy, (e) coverage anomalies discussion.

## Critical Findings

### CRIT-M1. Table 1 (Component Fragility, simulations.tex line 105) numbers do not match data files for several entries.

Verified by recomputing per-component max relative bias and coverage from sensitivity_sweep_c1.rds, sensitivity_sweep.rds (the n=500 5-component C2 file), and sensitivity_sweep_c3.rds.

| Paper Table 1 (sim sweep) | Actual data |
|---|---|
| C1 k_5: 41%, cov 0.69, alpha=0.5 | 42.7%, cov 0.455 (alpha=0.5) |
| C1 lambda_1: 29%, cov 0.91, alpha=0.5 | 36.6%, cov 0.985 (alpha=0.5) |
| C1 k_4: 21%, cov 0.83, alpha=0.5 | 20.1%, cov 0.540 (alpha=0.45) max |
| C2 lambda_2: 32%, "near nominal", s=1.0 | 41% at s=1.0 OR 93% at s=0.5 (cov 0.98) |
| C2 k_3: 24%, cov 0.34, s=1.0 | 37.1%, cov 0.695 (s=1.0) |
| C2 lambda_3: 18%, cov 0.66, s=1.0 | 14.2%, cov 0.665 (s=1.0) close |
| C3 lambda_5: 149%, cov 0.72, alpha=2.0 | 161% at alpha=2.0 OR 162.6% at alpha=1.25 |
| C3 lambda_3: 123%, cov 1.00, alpha=1.5 | 142.8%, cov 1.000 (alpha=1.5) |
| C3 lambda_1: 30%, cov 0.00, alpha=2.0 | 29.1%, cov 0.000 (alpha=2.0) matches |

Eight out of nine table entries differ by more than rounding error, and the coverage figures for C1 k_5 (0.69 in table, 0.455 in data) and C1 lambda_1 (0.91 in table, 0.985 in data) are off by very large margins. This is a reproducibility-killing discrepancy that any Technometrics reviewer will catch on a code spot-check. Suggestion: regenerate Table 1 from current data files; verify "worst" definition (per-component max across full sweep range vs at fixed maximum severity).

### CRIT-M2. Tables 2 and 3 slope-test numbers DO match the data; this contrasts with Table 1's drift.

Verification:
- C1 n=500: paper -0.32 p<0.001 | data -0.3225 p=0 (matches)
- C2 n=500: paper -0.04 p=0.256 | data -0.0413 p=0.256 (matches)
- C3 n=500: paper +0.004 p=0.689 | data +0.0036 p=0.689 (matches)
- C1 n=2000: paper -0.30 p<0.001 | data -0.3005 p=0 (matches)
- C2 n=2000: paper -0.01 p=0.305 | data -0.0142 p=0.305 (matches)
- C3 n=2000: paper -0.008 p=0.204 | data -0.0079 p=0.204 (matches)

So Tables 2 and 3 were regenerated from the current data, but Table 1 appears to be from an older (perhaps the C2-only redesign era) run. The C3 row of Table 1 also looks like it may have been hand-edited from a stale source. Fix: re-extract Table 1 from the same data the slope tests use.

## Major Findings

### MAJ-M1. State file is stale (says C1 "planned," C3 "script-ready" with smoke-only B=10), but data files show C1 and C3 sweeps were run at B=200, n=500 on 2026-05-24 and 2026-05-27.

The actual data files:
- sensitivity_sweep_c1.rds (B=200, dated 2026-05-27) -- COMPLETE
- sensitivity_sweep_c1_n2000.rds (B=100, dated 2026-05-27) -- COMPLETE
- sensitivity_sweep_c3.rds (B=200, dated 2026-05-24) -- COMPLETE
- sensitivity_sweep_c3_n2000.rds (B=100, dated 2026-05-27) -- COMPLETE
- sensitivity_sweep.rds (B=200, n=500, 5-comp, dated 2026-03-14) -- COMPLETE
- sensitivity_sweep_n2000.rds (B=100, dated 2026-05-27) -- COMPLETE

Action: update state.md to reflect that all three sweeps (main and ablation) are complete. Mark experiment status as "complete" instead of "planned"/"script-ready."

### MAJ-M2. Sample-size ablation (B=100 at n=2000) under-powers the small-effect C3 case.

Quoted from Table 3: "C3 ... n=2000 B=100 ... linear coef -0.008, p = 0.204."

For a slope of magnitude 0.008 against typical Monte Carlo noise sigma ~ 0.1, with only B=100 replications per severity level and 9 levels, the power to detect the slope is well under 50 percent. The "preserved" verdict for C3 at p=0.204 (n=2000) is a fail-to-reject, not a positive demonstration of preservation. Suggestion: explicitly disclose this as an underpowered test, or rerun the n=2000 C3 ablation at B=200+ to match the main run.

### MAJ-M3. Guo turbine application uses n=30 and m=3; the sensitivity-analysis demonstration is qualitative and the system-hazard stability (0.00302 to 0.00308) is at scale comparable to Monte Carlo variation given n=30.

Quoted from application.tex line 50: "the system hazard is remarkably stable (varying by less than 2% across the full severity range)."

With n=30 the sampling SE of the system hazard estimate is large; the 2% range observed across severity sweeps is plausibly within the noise floor of a single dataset. The application is fine as a workflow demonstration, but the manuscript should explicitly say so and not lean on the 2% range as quantitative validation of system-hazard robustness. Suggestion: frame as "illustrative workflow" rather than "real-data validation."

### MAJ-M4. The exponential exact-preservation result (Theorem 3.4 and Cor 3.5) is the strongest theoretical guarantee, and the C2 sibling-paper citation use case lands here cleanly: for series-system reliability with exponential components under any masking violation, the system hazard is exactly preserved. For Weibull, the empirical bounds (4% bias under C1 at alpha=0.5 of severity range; first-order robust under C2/C3 with empirical slopes 0.04 and 0.004) are clean for citation.

## Minor Findings

### MIN-M1. Coverage anomalies subsection (simulations.tex line 167) adequately discusses lambda_1 collapse under C3 and lambda_2 inflation at s=0.5; this addresses the 2026-03-14 review's CRIT-C1.

### MIN-M2. Direction matrix D is described as "cyclic structure with rows summing to zero (the explicit form is given in the supplementary code)." A reproducibility-conscious reviewer would prefer the explicit 5x5 D matrix in the appendix; current reference to "supplementary code" routes through run_paper_simulations.R or run_sensitivity_sweep.R. Recommend adding the explicit matrix to appendix.

### MIN-M3. Bernoulli model with independent inclusion of each non-failed component is acknowledged as a parametric perturbation family in Limitations (discussion.tex section 5.3, item 1). Adequate.

### MIN-M4. The Wald-vs-sandwich variance issue is acknowledged in Limitations (discussion.tex item 5, "Variance estimation"). Adequate disclosure, even though the sandwich version is not implemented.

## Cross-Verification

I executed all simulation data inspections via Rscript directly on the data files. The verdict that Table 1 is out of sync with the data files but Tables 2 and 3 match is data-grounded.

For the sibling-paper question (can the five coarsening papers cite this for "the cost of C2 violations"): yes, for exponential components (Cor 3.5: exact preservation), and yes for Weibull C2 (slope -0.04, p=0.256 at n=500; -0.014, p=0.305 at n=2000, on a 5-component Weibull system with sample sizes 500 and 2000). The C2-robust claim for sibling-paper citation is solid.
