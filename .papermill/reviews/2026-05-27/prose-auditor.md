# Prose Auditor Report

**Date**: 2026-05-27
**Paper**: Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption

## Critical Findings

### CRIT-P1. Abstract still says "C2 violation" but the body covers unified C1/C2/C3.

Quoted from main.tex line 115-133: "Maximum likelihood estimation for series systems with masked failure data relies on the non-informative masking assumption (Condition C2)... We characterize the bias that arises when C2 is violated but assumed to hold... Using a Bernoulli masking model that generates controlled C2 violations of varying severity, we conduct a systematic simulation sweep on a 5-component Weibull series system. The results establish a practical breakdown boundary..."

The abstract is the pre-redesign abstract preserved from the C2-only era. It does not mention C1 or C3 violations, the robustness hierarchy, the exponential exact-preservation theorem, robustness intervals, or the mdrelax R package. The introduction and body have all been rewritten to the unified framework; the abstract has not. Sibling-paper readers consulting the abstract will not see the unified-coarsening framing that drives the cross-series citation.

Suggested rewrite: "Maximum likelihood estimation for series systems with masked failure data relies on three coarsening conditions on the candidate-set diagnostic: that the true cause is in the set (C1), that masking is non-informative within the set (C2), and that masking is parameter-independent (C3). We develop a unified sensitivity framework that parameterizes each violation on a calibrated severity scale, derives closed-form first-order bias expansions of the misspecified MLE, and establishes a robustness hierarchy: system-level estimands (total hazard, MTTF) are preserved exactly for exponential components and at first order for Weibull components under C2 and C3, with a small bounded bias under C1; component-level estimands are fragile under all three. A 5-component Weibull simulation sweep confirms the hierarchy at sample sizes 500 and 2000. Methods are available in the mdrelax R package."

## Major Findings

### MAJ-P1. Contribution count mismatch.

The introduction (sections/introduction.tex lines 68-113) lists exactly five contributions:
1. Calibrated severity scales
2. Local sensitivity indices for set-valued masking
3. Robustness intervals
4. Robustness hierarchy across estimands
5. Software

The conclusion (sections/conclusion.tex lines 12-56) lists exactly five contributions with the same five categories (numbered identically).

The state file `.papermill/state.md` (line 50) describes "Six contributions" including "(4) formal specification test for {C1, C2, C3}, addressing the ad-hoc-GOF gap in Moreno-Betancur et al. (2015) and the screening-tool-not-test framing of ISNI." This contribution does not appear in either the introduction enumeration or the conclusion enumeration. The specification test is mentioned only in Future Work (discussion.tex line 166-180): "We have implemented and validated the singleton version of the test in our companion software; substantial development is needed for a publication-ready version."

Implication: the state file overcounts by one. Either (a) the manuscript should add the specification test as contribution #6 and bring it from Future Work into a proper section, or (b) the state file should be corrected to "Five contributions" and the specification-test claim in state.md should be downgraded.

Risk for Technometrics: a reviewer reading the paper expects five contributions; the abstract should be aligned to the five (currently it lists none explicitly). Sibling-paper authors who consulted the state file would have inferred six contributions; they should be informed of the actual contribution count.

### MAJ-P2. C1/C2/C3 vocabulary alignment with cross-series canonical phrasing.

The cross-series canonical C2 statement (per scrna-coarsening, phenotype-coarsening, spatial-coarsening, weaksup-coarsening, dp-coarsening translations) is verbal:

> "conditional on the candidate set, the masking distribution does not depend on which element is the true cause"

The mdrelax paper's Condition C2 (background.tex line 66-72) is given as an equation:

> P{C_i = c | T_i = t, K_i = j} = P{C_i = c | T_i = t, K_i = j'}, for all j, j' in c.

These are equivalent: the verbal statement says "conditional on C = c (and T = t), P{C = c | T, K} does not depend on K within {j in c}," which is exactly the equation. But the paper does NOT include a one-sentence verbal restatement that would let sibling-paper readers immediately match the equation to the canonical verbal form.

Suggestion: insert after Condition C2: "In words, conditional on the candidate set (and the failure time), the masking distribution does not depend on which element of the set is the true cause. This is the form of the condition cited in the broader coarsening literature." This makes the cross-series citation chain explicit.

### MAJ-P3. "Robustness hierarchy" usage is consistent. The term anchors the paper.

Used in: introduction (contribution #4), sensitivity_framework.tex (Cor 3.11 title and items), simulations.tex (validation), discussion.tex section 5.1 title, conclusion item #4. Cleanly consistent vocabulary.

## Minor Findings

### MIN-P1. Notation `ell^wrong` for misspecified log-likelihood.

The superscript "wrong" appears throughout sensitivity_framework.tex and appendix.tex. While intuitive, Technometrics style typically prefers `ell^mis`, `ell_M`, or `ell^star` for the misspecified model. "wrong" is informal and could be perceived as colloquial. Replace globally.

### MIN-P2. Application section (sections/application.tex, 65 lines) reads as briefly tacked on.

The 2026-03-14 review flagged "no real-data example" as a critical issue. The current Guo et al. application is added but spends only 65 lines on a single sensitivity sweep with five severity levels and no interpretive discussion of the component-level shifts. A reviewer comparing with the simulation section's 339 lines may judge the application underweight. Suggestion: expand to include a robustness-interval table for system hazard, a brief discussion of how a practitioner would decide whether to publish lambda_3's shift from 847 to 780 as "real" or "C2-sensitive."

### MIN-P3. Bernoulli model reuse for all three violation types is signposted adequately.

Definition 3.1 (Bernoulli Masking Model) is introduced in section 3.2 (in the C2 context). Section 3.6 introduces the "Bernoulli C1 violation model" by parametrizing eq (3.41) within the same Bernoulli family. Section 3.7 reuses the same family for the power-weight C3 model. The reuse is signalled but could be reinforced with a one-paragraph "Master Bernoulli model" subsection at the top of section 3 making the multi-violation structure explicit.

### MIN-P4. Future Work paragraph on specification test (discussion.tex line 166-180).

Quoted: "We have implemented and validated the singleton version of the test in our companion software; substantial development is needed for a publication-ready version."

This tense and voice mixes self-deprecation with implementation claim. Either remove the "we have implemented" clause (the contribution is properly Future Work) or strengthen it into a contribution section. Currently it falls awkwardly between.

### MIN-P5. "vlambda" macro defined but "lambda_j" used throughout simulations.

main.tex line 67 defines \vlambda. Used in introduction.tex line 67 (correct usage for vector). simulations.tex uses individual lambda_j (scalar components of vector). Notation is internally consistent; no issue.

### MIN-P6. Cross-references all resolve. Verified via grep for label and ref matches.

Labels referenced in main sections all exist within main sections; four labels (thm:ident-C123, eq:bernoulli-prob-general, eq:power-masking, eq:factorization) exist only in pre-redesign files (identifiability.tex, relaxed_models.tex, specification_test_sketch.tex) but are not referenced from main sections, so the dangling pre-redesign files do not produce build warnings. Build is clean (zero LaTeX warnings).

## Cross-Verification

Verified the contribution-count claim by reading both introduction.tex and conclusion.tex enumeration items directly. Verified the C2 statement consistency by reading background.tex Condition C2 and cross-checking against scrna-coarsening/sections/translation.tex line 52 and phenotype-coarsening/sections/translation.tex line 72. Verified abstract content against main.tex.
