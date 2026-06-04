# Multi-Agent Review Report

**Date**: 2026-05-27
**Paper**: Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption (towell2026mdrelax)
**Target Venue**: Technometrics; backups IEEE TR, Reliability Engineering and System Safety, Lifetime Data Analysis
**Recommendation**: major-revision

## Summary

**Overall Assessment**: The post-redesign unified C1/C2/C3 framework is conceptually sound and delivers on its core promise: a robustness hierarchy with exact exponential preservation (Theorem 3.4) anchoring a tiered classification of estimands. The five sibling coarsening papers can cite this work as the canonical reference for "the cost of C2 violations," especially for exponential components where the preservation is exact and at all orders. However, three blockers remain before Technometrics submission: Table 1's simulation numbers do not match the actual data files, the abstract is the pre-redesign C2-only version, and two claimed contributions (robustness intervals and the specification test) are asserted but not delivered in the manuscript text.

**Strengths**:
1. Theorem 3.4 (Exact Total-Hazard Preservation for exponential via profile-likelihood factorization) is a clean new structural result that strengthens prior first-order claims (logic-checker, novelty-assessor).
2. Theorem 3.7 (C1 closed-form m=2 expansion) is algebraically correct and the result matches the empirical slope (-0.32) in Tables 2 and 3 (logic-checker).
3. Robustness hierarchy (Cor 3.11) is the natural sibling-paper citation target and is well-supported by the underlying theorems (novelty-assessor).
4. Slope tests in Tables 2 and 3 reproduce exactly from the data files (methodology-auditor).
5. Bibliography is well-targeted to the reframed thesis; CAR/ISNI/biostatistics missing-cause corpus is fully engaged (citation-verifier).
6. Build is clean: zero LaTeX warnings, all cross-references resolve, bibliography compiles (format-validator).
7. Non-nesting lemma (Lemma A.1) for set-valued vs point-cause masking is a useful original result that supports the differentiation from Bakoyannis 2025 and Moreno-Betancur 2015 (logic-checker, novelty-assessor).

**Weaknesses**:
1. Table 1 (component fragility) numbers do not match the data files for several entries (methodology-auditor).
2. Abstract is the pre-redesign C2-only version (prose-auditor).
3. Contribution count off-by-one: state file claims six contributions, manuscript enumerates five; specification test is in Future Work, not in the contribution list (prose-auditor, novelty-assessor).
4. Robustness intervals (contribution #3) named but never formally defined in the paper text (novelty-assessor).
5. Cross-package figure path `../inst/simulations/figures/` will break Technometrics camera-ready submission (format-validator).
6. ORCID missing from main.tex title block (format-validator).
7. State file (.papermill/state.md) is stale: claims C1 sweep "planned" and C3 "script-ready," but data files show both are complete with B=200 (methodology-auditor).
8. Page count at 39 pages double-spaced is over typical Technometrics limit (format-validator).

**Finding Counts**: Critical: 2 | Major: 9 | Minor: 11 | Suggestions: 5

## Critical Issues

### CRIT-1. Table 1 numbers do not match the data files (source: methodology-auditor; cross-verified via direct Rscript inspection of paper/data/*.rds)

- **Location**: simulations.tex line 105 (Table 1, tab:component-fragility)
- **Quoted from manuscript**:
```
C1 & k_5        & 41%   & 0.69          & alpha_C1 = 0.5
   & lambda_1  & 29%   & 0.91          & alpha_C1 = 0.5
   & k_4        & 21%   & 0.83          & alpha_C1 = 0.5
C2 & lambda_2  & 32%   & near nominal  & s = 1.0
   & k_3        & 24%   & 0.34          & s = 1.0
   & lambda_3  & 18%   & 0.66          & s = 1.0
C3 & lambda_5  & 149%  & 0.72          & alpha_C3 = 2.0
   & lambda_3  & 123%  & 1.00          & alpha_C3 = 1.5
   & lambda_1  & 30%   & 0.00          & alpha_C3 = 2.0
```
- **Actual data** (verified directly from sensitivity_sweep_c1.rds, sensitivity_sweep.rds, sensitivity_sweep_c3.rds, all B=200 n=500 5-component Weibull):

| Row | Paper | Data |
|---|---|---|
| C1 k_5 | 41%, cov 0.69 | 42.7%, cov 0.455 |
| C1 lambda_1 | 29%, cov 0.91 | 36.6%, cov 0.985 |
| C1 k_4 | 21%, cov 0.83 | 20.1%, cov 0.540 (at alpha=0.45) |
| C2 lambda_2 | 32%, cov "near nominal" | 41% at s=1.0 OR 93% at s=0.5 (cov 0.98) |
| C2 k_3 | 24%, cov 0.34 | 37.1%, cov 0.695 |
| C2 lambda_3 | 18%, cov 0.66 | 14.2%, cov 0.665 (matches) |
| C3 lambda_5 | 149%, cov 0.72 | 162.6% at alpha=1.25 OR 161% at alpha=2.0 |
| C3 lambda_3 | 123%, cov 1.00 | 142.8%, cov 1.000 (alpha=1.5) |
| C3 lambda_1 | 30%, cov 0.00 | 29.1%, cov 0.000 (matches) |

Eight of nine entries are off by more than rounding error; coverage for C1 k_5 is off by 0.24 (0.69 in paper vs 0.455 in data). This will fail a Technometrics reproducibility check.

- **Suggestion**: Regenerate Table 1 from current data using the same code path as Tables 2 and 3. Document the "worst" criterion (per-component max across full sweep range, or at fixed maximum severity). Note that Tables 2 and 3 DO match the data, so the issue is isolated to Table 1.
- **Cross-verified**: Yes, by direct Rscript inspection of all three data files.

### CRIT-2. Theorem 3.9 (C3 misspecification) and contribution #3 (Robustness intervals) defer key content to nonexistent supplementary material (source: novelty-assessor, logic-checker)

- **Location**: sensitivity_framework.tex line 532-548 (Theorem 3.9); conclusion.tex line 33-39 (contribution #3)
- **Quoted text (Theorem 3.9)**: "Closed-form expressions for the individual coefficients Delta^C3_j with m = 2 and general m are given in the supplementary material."
- **Quoted text (conclusion)**: "The implementation (Monte-Carlo and first-order variants) is available in the mdrelax R package."
- **Problem**: No supplementary material file exists in the repository. The robustness interval construction is named as a contribution but never formally defined in the paper text. For a Technometrics methodology paper, the methodology contribution must be defined in the paper itself, not delegated to a software package.
- **Suggestion**: Either (a) move the deferred derivations into Appendix B and define robustness intervals formally in Section 3.8, or (b) downgrade contribution #3 from "robustness intervals" to "robustness hierarchy informed by ISNI" and remove the supplementary-material references. Option (a) requires writing approximately 2-3 pages; option (b) requires only contribution-list editing.
- **Cross-verified**: Yes, by reading sensitivity_framework.tex Theorem 3.9 and conclusion.tex.

## Major Issues

### MAJ-1. Abstract is the pre-redesign C2-only version (source: prose-auditor)

- **Location**: main.tex lines 115-133
- **Quoted text**: "Maximum likelihood estimation for series systems with masked failure data relies on the non-informative masking assumption (Condition C2)... We characterize the bias that arises when C2 is violated but assumed to hold..."
- **Problem**: The body of the paper has been reframed (2026-05-04) to a unified C1/C2/C3 framework with five contributions and an exact-preservation theorem for exponential components. The abstract still describes only C2 violations on a 5-component Weibull system, with no mention of C1, C3, the exact-preservation theorem, the robustness hierarchy, or the mdrelax package. A reader (or sibling-paper citer) consulting the abstract will not see the unified framework that anchors the cross-series citation use.
- **Suggestion**: Rewrite the abstract to mention C1, C2, C3 unification, the profile-likelihood exact-preservation result for exponential, the first-order robustness for Weibull system-level estimands, and the robustness hierarchy.

### MAJ-2. Contribution count off-by-one (source: prose-auditor, novelty-assessor)

- **Location**: introduction.tex lines 68-113 (five enumerated contributions); conclusion.tex lines 11-56 (five enumerated contributions); state.md line 50 ("Six contributions")
- **Problem**: The state file describes six contributions including a "formal specification test for {C1, C2, C3}." Neither the introduction nor the conclusion lists six. The specification test appears only in discussion.tex Future Work (line 166-180): "substantial development is needed for a publication-ready version."
- **Suggestion**: Reconcile by either (a) developing the specification test into a Section 4 contribution and aligning all three lists at six, or (b) updating the state file to "Five contributions" and downgrading the specification-test claim. Sibling-paper authors who pre-cited Towell 2026 for the specification test based on the state file need to be informed.

### MAJ-3. Robustness intervals contribution not delivered in paper text (source: novelty-assessor)

- **Location**: introduction.tex lines 89-96, conclusion.tex lines 33-39, discussion.tex
- **Problem**: The robustness interval construction is referenced multiple times but never formally defined. The reader should expect: "Given tolerance epsilon, the robustness interval for estimand psi at violation type C_k is RI_psi(epsilon) = {alpha_C_k : |psi(theta^dagger(alpha_C_k)) - psi(theta^*)| < epsilon}" or analogous formulation. The paper does not provide this.
- **Suggestion**: Add a brief subsection in Section 3.8 (Local Sensitivity Indices and Robustness Hierarchy) formally defining robustness intervals and exhibiting them in the application section.

### MAJ-4. State file is stale relative to data files (source: methodology-auditor)

- **Location**: .papermill/state.md experiments section
- **Quoted from state**: "C1 sweep: status: planned"; "C3 sweep: status: script-ready" (smoke test only B=10)
- **Actual**: paper/data/sensitivity_sweep_c1.rds (B=200, dated 2026-05-27), sensitivity_sweep_c1_n2000.rds (B=100, 2026-05-27), sensitivity_sweep_c3.rds (B=200, 2026-05-24), sensitivity_sweep_c3_n2000.rds (B=100, 2026-05-27). All sweeps are complete.
- **Suggestion**: Update state.md experiment statuses to "complete" with the actual replication counts and timestamps.

### MAJ-5. Cross-package figure path will break Technometrics submission (source: format-validator)

- **Location**: main.tex line 18: `\graphicspath{{../inst/simulations/figures/}}`
- **Problem**: For ManuscriptCentral/ASA tarball submission, figures should be co-located with the paper source. The current path pulls from the R package's inst/ directory.
- **Suggestion**: Modify run_sensitivity_sweep*.R scripts to write figures to paper/figures/ instead of inst/simulations/figures/, then update graphicspath. State file already flags this.

### MAJ-6. Corollary 3.8 (general-m C1 hazard preservation) proof has a minor gap (source: logic-checker)

- **Location**: sensitivity_framework.tex lines 460-480, Corollary 3.8 proof
- **Problem**: The retention-independence argument relies specifically on the Bernoulli C1 parametrization (3.41) where off-diagonal probabilities are constant in k. If off-diagonal varied with k, the retention probability would be k-dependent and the marginal-preservation argument would fail. The proof should make this dependence explicit.
- **Suggestion**: Insert one line: "Because (3.41) leaves the off-diagonal probabilities constant in k, the retention event {|C| >= 1} is independent of both k and t, hence independent of T marginally."

### MAJ-7. Pre-redesign remnant section files clutter paper/sections/ (source: format-validator)

- **Location**: paper/sections/relaxed_models.tex (502 lines), identifiability.tex (382 lines), introduction.pre-unified.tex (74 lines), specification_test_sketch.tex (207 lines)
- **Problem**: These files are not included in main.tex and contain orphan labels (thm:ident-C123, eq:bernoulli-prob-general, eq:power-masking, eq:factorization). They do not affect the build but will confuse reviewers who inspect the source tarball.
- **Suggestion**: Move to paper/sections/_unused/ or delete before submission. State file already flags this.

### MAJ-8. C3 n=2000 ablation underpowered (source: methodology-auditor)

- **Location**: Table 3 (simulations.tex line 277)
- **Quoted**: "C3 ... n=2000 ... -0.008 ... 0.204"
- **Problem**: With B=100 replications and slope magnitude 0.008, the slope test is severely underpowered. The "preserved" verdict (fail to reject) is not a positive demonstration of preservation; it is a Type II failure-to-detect.
- **Suggestion**: Rerun n=2000 C3 ablation at B=200 to match the main run, or explicitly disclose the power limitation.

### MAJ-9. ORCID missing from main.tex title block (source: format-validator)

- **Location**: main.tex lines 97-101
- **Problem**: DESCRIPTION lists ORCID 0000-0001-6443-9897; main.tex author block does not. State file polish checklist flags this.
- **Suggestion**: Add `\thanks` with department and ORCID to the author block.

## Minor Issues

### MIN-1. C1/C2/C3 vocabulary alignment with cross-series canonical statement (prose-auditor)

The cross-series canonical C2 statement is verbal: "conditional on the candidate set, the masking distribution does not depend on which element is the true cause." The paper's Condition C2 (background.tex line 67-72) gives the equivalent equation. Insert a one-sentence verbal restatement to make the cross-series citation chain explicit.

### MIN-2. Notation `ell^wrong` (prose-auditor)

Replace globally with `ell^mis` or `ell^M`. The current notation is informal.

### MIN-3. Add `\citep{towell2026masked}` to background.tex (citation-verifier)

The foundational masked-causes-in-series-systems paper is the natural reference for C1-C2-C3 concepts. Currently the manuscript cites only towell2023reliability (master's thesis) and Lin-1993. Adding the foundational paper at background.tex lines 11 and 59 strengthens the citation chain across the series.

### MIN-4. Vanity count flagged (format-validator)

appendix.tex line 174: "1276 unit tests verifying..." Replace with: "The package's test suite verifies likelihood correctness, score-gradient agreement, FIM consistency, and MLE convergence properties." (Drop the integer.)

### MIN-5. Document class is article, not ASA template (format-validator)

State file polish checklist flags retemplating to Technometrics ASA template before submission.

### MIN-6. Page count over Technometrics limit (format-validator)

39 pages double-spaced, typically over Technometrics regular-article limit. Consolidate Limitations subsection; move Software appendix to supplementary materials.

### MIN-7. Hyperref Unicode warnings in main_build.log (format-validator)

Two "Token not allowed in a PDF string (Unicode)" warnings. Cosmetic; resolve with `\texorpdfstring`.

### MIN-8. Cor 3.11 item 3 mixes derived and empirical claims (logic-checker)

"empirical magnitude below 4% in the 5-component sweep of Section 4" is an empirical observation inside a corollary. Separate the theoretical statement (first-order coefficient is nonzero for Weibull C1) from the empirical bound (4% in this configuration).

### MIN-9. Cinelli-Hazlett "multi-dimensional sensitivity geometry" is aspirational (novelty-assessor)

The paper's three sweeps are independent univariate. The "geometry" claim is overstated; rephrase as "calibrated severities across three independent violation axes."

### MIN-10. Application section is brief (prose-auditor)

65 lines on the Guo et al. m=3 turbine analysis vs 339 lines of simulation. Add a robustness-interval table and a brief decision discussion ("how does a practitioner decide whether lambda_3's shift from 847 to 780 is a real change or sensitivity to C2 violation?").

### MIN-11. html_paper/ directory contains LaTeXML output (format-validator)

Move to a build artifacts directory or gitignore.

## Suggestions

1. Lead with Theorem 3.4 (exact preservation for exponential) as the headline novelty. This is the strongest novel structural finding and currently appears mid-Section 3 after the C2 misspecification discussion. Restructure Section 3 to make this theorem prominent.

2. Add a "Master Bernoulli Model" paragraph at the top of Section 3 noting that all three violation types live within the same Bernoulli substrate. This unifies the presentation and the reuse becomes a unification feature rather than an ad hoc convenience.

3. Once the C1 and C3 sweep data are documented as complete in state.md, all five sibling coarsening papers can confidently cite this paper for the cost of C2 violations using the language: "For series-system reliability under exponential components, the system hazard is exactly preserved under arbitrary masking violations (Towell 2026, Corollary 3.5). For Weibull components, the bias is bounded and first-order preserved under C2 (Towell 2026, Corollary 3.11)."

4. Consider adding a brief paragraph to the discussion clarifying that the towell2026mdrelax paper provides the sensitivity bounds that the sibling coarsening papers (scrna-coarsening, spatial-coarsening, dp-coarsening, weaksup-coarsening, phenotype-coarsening) rely on. This bidirectional citation makes the research-series structure visible to readers.

5. Run the formal specification test development now if possible; this would convert Future Work item #1 into a sixth contribution and would address the state file's "six contributions" framing without forcing a downgrade.

## Detailed Notes by Domain

### Logic and Proofs (logic-checker.md)

Core theorems are correct: Theorem 3.4 (exact preservation via profile factorization) is rigorous and clean; Theorem 3.7 (C1 m=2 closed form) is algebraically correct; Theorem 3.6 (non-identifiability) constructive proof verifies. The C2 misspecification Theorem 3.2 has the "approximately" issue addressed in the 2026-03-14 review. New issues: Corollary 3.8 (general-m C1) needs a one-line clarification of why retention is k-independent; Theorem 3.9 (C3) defers details to nonexistent supplementary material; Cor 3.11 item 3 mixes derived with empirical claims.

### Novelty and Contribution (novelty-assessor.md)

The robustness hierarchy (Cor 3.11) is the strongest organizational contribution and the natural sibling-paper citation anchor. Theorem 3.4 (exact preservation via profile factorization) is the strongest novel structural finding but is currently positioned after the C2 misspecification discussion, which understates its centrality. Two contributions (robustness intervals, specification test) are asserted but not delivered in the manuscript text; alignment with state file's "Six contributions" framing requires either developing the missing content or downgrading the contribution list.

### Methodology (methodology-auditor.md)

The simulation design is sound: 5-component Weibull, three independent sweeps, two sample sizes. The slope tests in Tables 2 and 3 reproduce exactly from the data. Table 1 numbers do not match the data (Critical 1). C3 n=2000 ablation is underpowered. The Guo et al. application (n=30, m=3) is appropriate as a workflow demonstration but the "2% range" is at the noise floor of n=30.

### Writing and Presentation (prose-auditor.md)

The post-redesign prose is generally clean and consistent. Major issues: abstract still C2-only; contribution count off-by-one between state and manuscript; cross-series C2 vocabulary alignment needs a one-sentence verbal restatement.

### Citations and References (citation-verifier.md)

All 35 citation keys resolve. Bibliography is well-targeted. Add `\citep{towell2026masked}` for the foundational paper. Self-citation to master's thesis (towell2023reliability) is acceptable but the published towell2026masked is preferred where applicable.

### Formatting and Production (format-validator.md)

Build is clean (zero warnings). Cross-package figure path will break submission; pre-redesign remnant files clutter source tree; ORCID missing; document class still article (not ASA). State file flags all of these as pending.

## Cross-Series Citation Readiness

**Can the five sibling coarsening papers (scrna-coarsening, spatial-coarsening, dp-coarsening, weaksup-coarsening, phenotype-coarsening) cite this paper as-is for "the cost of C2 violations"?**

Yes, with one minor caveat. The Theorem 3.4 + Corollary 3.5 chain provides exact preservation of the system hazard under any C2 violation for exponential components, which is the strongest possible sibling-paper citation. The Weibull case relies on Cor 3.11 item 3 (first-order robust under C2) plus Table 2 empirical slope of -0.04 (p=0.256). The caveat: Table 1's component-level fragility numbers do not match the data files, so sibling papers should cite Table 2 (slope tests) and Cor 3.11 (hierarchy) rather than Table 1's specific component-bias percentages.

## Literature Context Summary

The paper engages well with the reframed corpus: Heitjan-Rubin CAR foundations, Jacobsen-Keiding CAR generalization, Troxel-Ma-Heitjan and Zhang-Heitjan ISNI, Copas-Eguchi and Siannis-Copas-Lu local bias, Cinelli-Hazlett robustness value, Moreno-Betancur missing-cause pattern-mixture, Bakoyannis-Yiannoutsos cumulative incidence misclassification, Bakoyannis 2025 robustness intervals (concurrent work), Ebrahimi/Van Rompaye/Ha-Tsodikov point-cause misclassification (non-nested per Lemma A.1). The reliability dependent-masking lineage (Lin-Guess, Guttman, Mukhopadhyay, Craiu-Reiser, Flehinger) is also fully engaged. No major gaps relative to the paper's positioned scope.

## Review Metadata
- Agents involved: literature-context (auto-generated), logic-checker, methodology-auditor, prose-auditor, novelty-assessor, citation-verifier, format-validator
- Direct data inspection: Rscript on paper/data/*.rds to verify Tables 1, 2, 3 against the actual simulation outputs
- Cross-verifications performed: 4 (Table 1 numbers vs data files; slope tests vs data; contribution count vs state file; C2 verbal vs equation form vs cross-series canonical)
- Disagreements noted: 0 (specialists converged on the same blockers)
