# Logic Checker Report

**Date**: 2026-06-08
**Focus**: proof correctness, logical chain integrity, claim support. This pass re-verifies the load-bearing proofs against the current source and the abstract's quantitative claims against committed data, and checks internal consistency of scope statements.

## Headline: no critical or major logic defects

The proof core is sound and unchanged from the verified 2026-06-04 state. I re-read the two most load-bearing results (exact-preservation theorem, robustness-coverage proposition) line by line against the current files and independently reproduced every numeric claim the proofs are advertised to support.

## Re-verified proofs (current text)

- **Theorem 3.6 / thm:exact-hazard-exp (exact total-hazard preservation, exponential).** Sound. The profile-likelihood factorization ell(S, phi) = ell_S(S) + ell_phi(phi) (eq:like-factored) is correct: substituting theta_k = S phi_k gives log(sum_{k in c} S phi_k) = log S + log(sum_{k in c} phi_k), with exactly n_F failure terms contributing the n_F log S, so the S-score is phi-independent and S-hat = n_F / sum s_i. The one substantive hypothesis (retained T ~ Exp(S*)) is discharged by the corollaries. This is the structural headline and it holds.
- **Corollary 3.7 / cor:c2-hazard-exp, Corollary 3.14 / cor:c1-hazard-general.** Both correct. The C1 general-m retention probability Pr{|C|>=1 | K=k, T=t} = 1 - alpha_C1 (1-p_0)^{m-1} is independent of k and t, so retention is independent of T and Theorem 3.6 applies. The old MAJ-6 gap (general-m C1) stays closed.
- **Theorem 3.13 / thm:misspec-C1 (m=2 first-order expansion).** Algebra re-checked: the combined linear system forces Delta_1 + Delta_2 = 0 and yields Delta_1 = p_0(theta_2* - theta_1*)/(1-p_0), Delta_2 = -Delta_1 (eq:c1-delta). The amplification factor p_0/(1-p_0) and its divergence as p_0 -> 1 are consistent with Remark 3.15.
- **Theorem 3.18 / thm:misspec-C3.** Total preservation (sum_j Delta_j^C3 = 0, all-orders via thm:exact-hazard-exp) and the rank-(m-1) score-Jacobian system J Delta = r are internally consistent; the rank deficiency is correctly tied to the non-identifiability of Theorem 3.16.
- **Proposition 3.21 / prop:robustness-coverage.** Re-verified sound. f(alpha) = psi(theta-dagger(alpha)); the second derivative f'' = thetadot^T grad^2 psi thetadot + grad psi^T thetaddot is computed correctly and captures both curvature sources; the bound |f''| <= L_psi M_k + B_psi D_k^2 follows by Cauchy-Schwarz / operator-norm submultiplicativity; the linear term f'(0) = ISNI matches Definition 3.19. Honestly stated as a deterministic pseudo-true bias bound, not a frequentist coverage statement (Remark 3.23 makes this explicit).
- **Theorem 3.16 / thm:non-ident and Lemma C.2 / lem:non-nesting.** Both directions of each are correct (re-checked the m=2 identifiability construction and the m=4 point-cause-vs-set-valued counterexample). Lemma C.2 is the load-bearing differentiation result and it holds.

## Quantitative claims independently reproduced from committed data

Ran the committed `paper/regenerate_tables.R` and direct R computation on `paper/data/*.rds`:

- Abstract "system-hazard relative bias bounded under 4 percent throughout": confirmed. Max |sys-hazard rel bias| = 3.22% (C1, alpha=0.50), 2.95% (C2, s=0.90), 2.07% (C3, alpha=0.50). All < 4%.
- Abstract "relative bias exceeding 100 percent under C3": confirmed (lambda_5 = 163%, lambda_3 = 146%).
- Abstract "confidence-interval coverage collapsing to zero for some parameters": confirmed (lambda_1 under C3, coverage 0.00 at alpha=2.0).
- Tables 1/2/3: every cell reproduced exactly (43%/0.46, 93%/0.98, 163%/0.96; slopes C1 -0.32 p<0.001, C2 -0.04 p=0.256, C3 +0.004 p=0.689; n=2000 C1 -0.30, C2 -0.01, C3 -0.008).
- Application h_T(249.5) = 0.00307: reproduced exactly from the canonical Guo MLE.

No hallucinated or unsupported numeric claim found.

## Internal-consistency / scope finding (MINOR, new this pass)

The paper was reframed from a C2-only study to a unified C1/C2/C3 study, and two scope sentences did not migrate:

1. **background.tex line 84**: "This paper studies the consequences of C2 failure." Contradicts the abstract and introduction, which promise sensitivity to all three conditions. Quoted text confirmed present.
2. **sensitivity_framework.tex lines 6-7**: "We now develop the tools needed to study the sensitivity of C1-C2-C3 inference to C2 violations. We proceed in four steps:" The four-step plan that follows is C2-specific, even though the section later adds C1 (sec:c1-sensitivity) and C3 (sec:c3-sensitivity) subsections. The "to C2 violations" scope clause undersells the section's actual content.

Logically these are not errors in any theorem; they are stale framing that a referee will read as inconsistent with the unified contribution list. Fix: change line 84 to "...the consequences of C1, C2, and C3 failure" and reword the framework-intro scope clause to cover all three. Confidence: high (text quoted directly).

## Carryover logic nits (MINOR, unchanged)

1. Corollary 3.17 item 3 (cor:hierarchy) embeds a configuration-specific empirical number ("empirical magnitude below 4% in the 5-component sweep") inside a corollary. The theoretical content (Weibull C1 first-order coefficient is nonzero) and the empirical bound should be typographically separated; the 4% is data, not theorem.
2. Proposition 3.21 restates the L_psi, B_psi sup-definitions twice (in the hypotheses and again two lines later). Harmless duplication.
3. The word "Coverage" in the title of Proposition 3.21 may invite a referee to expect a 1-alpha frequentist coverage statement; the body is correct but a retitle ("Asymptotic bias bound for psi(theta-dagger)") would preempt the misread.
4. (Optional, carryover) The proposition would be cleaner with one sentence stating the uniform Hessian bound M_k exists whenever the expected information at theta* is non-singular (implicit function theorem applied to the score system), which holds generically; for Weibull the map is implicit and the bound is currently asserted rather than constructed.

## Bottom line (logic)

No critical or major logic defects; the proof core is submission-grade and every advertised numeric claim reproduces from committed data. The only new item is the residual C2-only scope language in Background and the Framework intro (MINOR, internal consistency). The carryover nits are all cosmetic.
