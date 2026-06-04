# Logic Checker Report

**Date**: 2026-06-04
**Focus**: proof correctness, logical chain integrity, claim support. Priority: the new Section 3.8 robustness-interval Definition + Proposition (unreviewed), and the C2-sensitivity core results.

## Verdict on the two previously-flagged Criticals

**CRIT-2 (robustness intervals named but never defined; supplementary-material forward reference) is FIXED.**
- Robustness intervals are now formally defined: Definition 3.20 (def:robustness-interval, sensitivity_framework.tex lines 610-630) and the local sensitivity index Definition 3.19 (def:isni, lines 587-601).
- Coverage/remainder is now a stated-and-proved Proposition 3.21 (prop:robustness-coverage, lines 632-681).
- The former Theorem 3.9 deferral to "supplementary material" is gone. The C3 result is now Theorem 3.18 (thm:misspec-C3), which presents the general-m linear system inline (eq:c3-linsys) rather than deferring. A full-text search for "supplementary"/"supplement" across all included sections returns nothing.

## Assessment of the NEW robustness-interval Proposition 3.21 (prop:robustness-coverage)

The argument is sound and submission-grade. Detail:

- Setup: f(alpha) = psi(theta-dagger(alpha)). The proof applies a first-order Taylor expansion with Lagrange remainder (1/2) f''(xi) alpha^2.
- The second derivative is computed correctly by the chain rule:
  f''(alpha) = thetadot^T grad^2 psi thetadot + grad psi^T thetaddot.
  This correctly captures BOTH curvature sources: the curvature of the pseudo-true map (thetaddot, bounded by M_k) and the curvature of the estimand (grad^2 psi contracted with thetadot-squared, bounded by B_psi D_k^2).
- The bound |f''(xi)| <= L_psi M_k + B_psi D_k^2 follows from the stated uniform bounds (||grad psi|| <= L_psi, ||grad^2 psi|| <= B_psi, ||thetaddot|| <= M_k, ||thetadot|| <= D_k) by Cauchy-Schwarz / operator-norm submultiplicativity. This is correct: the first term bounds grad psi^T thetaddot by L_psi M_k, the second bounds thetadot^T grad^2 psi thetadot by B_psi D_k^2.
- The linear term f'(0) = ISNI_{C_k}[psi] = grad psi(theta*)^T Delta^{C_k}(theta*) matches Definition 3.19 (eq:isni). Consistent.
- Conclusion: psi(theta-dagger(alpha)) lies within RI of half-width |ISNI| * sup(A_k) up to an O(sup(A_k)^2) remainder. This is honestly stated as a first-order interval, not an exact confidence interval.

**Are the assumptions (C^2, bounded gradient, bounded Hessian) sufficient for the stated bound?** Yes, for the stated conclusion, which is a deterministic bias bound on the pseudo-true value, not a finite-sample coverage statement. Two honesty caveats the paper already handles well, plus one residual gap:

1. (Handled) The proposition does NOT claim coverage of psi(theta*) by a random interval at level 1-alpha; it bounds the asymptotic pseudo-true bias. Remark 3.23 (rem:robustness-interval-usage) explicitly says the interval "bounds the asymptotic pseudo-true bias and does not absorb Monte Carlo uncertainty," and tells the practitioner to widen by the Wald half-width. This is the correct and honest framing; the title word "Coverage" in prop:robustness-coverage refers to covering the pseudo-true value, which the body makes clear.

2. (Handled) The C^2 / bounded-derivative assumptions on the map alpha -> theta-dagger(alpha) are asserted "by the regularity conditions of Theorems 3.7/3.13/3.18." Those theorems actually prove only first-order EXPANSIONS (existence of Delta_j and an O(alpha^2) error). The existence of a uniform Hessian bound M_k on a neighborhood is a slightly stronger statement than "the first-order expansion exists." For exponential components M_k is available in closed form (the pseudo-true map solves a smooth score system with non-singular Jacobian away from p_0 -> 1). 

3. (Residual, MINOR) The proposition would be cleaner if it stated where the uniform bounds M_k, D_k come from for the Weibull case. For exponential the score system is rational in theta with the information matrix non-singular on the open simplex, so smoothness and local boundedness are immediate. For Weibull the pseudo-true map is defined implicitly and its second derivative bound is asserted rather than constructed. This does not invalidate the proposition (the hypotheses are stated as assumptions), but a referee may ask for one sentence justifying that the implicit-function-theorem regularity holds at the standard-model fit (non-singular expected information), which it does generically. Suggested one-line addition after the hypotheses: "These bounds exist whenever the expected information at theta* is non-singular, by the implicit function theorem applied to the score system."

Net: the Proposition is correct as stated; the only actionable item is the optional sentence in (3).

## C2-sensitivity core results (priority verification)

All sound and mutually consistent.

- **Theorem 3.6 (thm:exact-hazard-exp), exact total-hazard preservation.** The profile-likelihood factorization ell(S, phi) = ell_S(S) + ell_phi(phi) (eq:like-factored) is correct: substituting theta_k = S phi_k into the exponential standard-model log-likelihood gives -S sum s_i + n_F log S + sum_{failures} log(sum_{k in c} phi_k), since log(sum_{k in c} S phi_k) = log S + log(sum_{k in c} phi_k) and there are exactly n_F failure terms. The S-score is then independent of phi, giving S-hat = n_F / sum s_i, the standard exponential-rate MLE. The only substantive hypothesis is that retained T ~ Exp(S*); the corollaries discharge it. This is the headline structural result and it is rigorous.
- **Corollary 3.7 (cor:c2-hazard-exp).** Correct: under any Bernoulli C2 model with C1 holding, |C| >= 1 a.s., no retention, T ~ Exp(S*), so Theorem 3.6 applies. Clean.
- **Proposition 3.4 (prop:hazard-robust).** Correctly stated as the exponential specialization; its proof correctly defers to thm:exact-hazard-exp. Scope is exponential-only and the text says so.
- **Theorem 3.13 (thm:misspec-C1), m=2 C1 first-order expansion.** The score-equation expansion is algebraically correct: setting Delta_1 + Delta_2 = 0 from the combined linear system and solving gives Delta_1 = p_0(theta_2* - theta_1*)/(1-p_0), Delta_2 = -Delta_1 (eq:c1-delta). The amplification factor p_0/(1-p_0) and its divergence as p_0 -> 1 are consistent with Remark 3.15. The empirical slope this predicts is consistent with the broken-first-order C1 verdict in Table 2.
- **Corollary 3.14 (cor:c1-hazard-general), general-m C1 total-hazard preservation.** This is the old MAJ-6 gap from the 2026-05-27 review, and it is now FIXED. The proof explicitly states Pr{|C| >= 1 | K=k, T=t} = 1 - alpha_{C1}(1-p_0)^{m-1}, independent of both k and t, so retention is independent of T and Theorem 3.6 applies. The one-line clarification the prior review asked for is present.
- **Theorem 3.18 (thm:misspec-C3).** The total-preservation item (sum_j Delta_j^{C3} = 0, all-orders via thm:exact-hazard-exp because the power-weight model makes |C| >= 1 deterministic) is correct. The general-m linear system J Delta = r with the rank-(m-1) observation tied to the non-identifiability of Theorem 3.16 is internally consistent. This is the result the prior review flagged as deferring to nonexistent supplementary material; the inline linear-system presentation resolves CRIT-2's logic half.

## Consistency with the foundational paper's C2 definition

Condition C2 in background.tex (cond:C2): for all c, all t, all j, j' in c, Pr{C=c | T=t, K=j} = Pr{C=c | T=t, K=j'}. Remark 3.2 (rem:C2-P) restates this as "each row of P has constant off-diagonal entries." This is the standard non-informative-masking / symmetry condition (Usher-Lin-Guess 1993 "symmetry" assumption) and is the same condition the foundational towell2026masked framework uses. The C1 and C3 sweep constructions are deliberately built to hold C2 fixed (C1 sweep: off-diagonal constant in k, comment at eq:c1-violation; C3 sweep: normalizer independent of k, comment at lines 510-514), which is exactly what "isolating one violation axis" requires. Logically clean.

## Other proofs spot-checked

- **Theorem 3.16 (thm:non-ident), non-identifiability of (theta, P).** The m=2 constructive family (theta_1' = theta_1*+eps, p_21' and p_12' chosen to preserve the three likelihood-determining functionals) is correct: S' = S* by construction, and the single-component contributions theta_1(1-p_21) and theta_2(1-p_12) are preserved by the explicit reparametrization. Supported empirically by Appendix B (joint MLE converges to (1.62,1.44) vs true (1.0,2.0) while sum-hat = 3.01 vs 3.00). Sound.
- **Lemma C.2 (lem:non-nesting), set-valued vs point-cause non-nesting.** Both directions are correct. Direction (i): Bernoulli set-valued admits empty C with positive probability; point-cause-plus-masking always has C containing K-obs, so C nonempty a.s. Direction (ii): the m=4 construction yields Pr{2 in C, 3 in C | K=1} = 0 under point-cause but = q_2 q_3 = 1/4 under independent-inclusion set-valued. This is the key differentiation result and it holds.

## Remaining logic-level nits (MINOR)

1. (Carryover MIN-8) Corollary 3.17 item 3 (cor:hierarchy) still embeds an empirical claim ("empirical magnitude below 4% in the 5-component sweep") inside a corollary statement. The theoretical content (Weibull C1 first-order coefficient is nonzero) and the empirical bound (under 4% in this configuration) should be typographically separated; the empirical number is configuration-specific, not a theorem.
2. Proposition 3.21 line "with L_psi = sup ... and B_psi = sup ..." repeats the definitions of L_psi, B_psi already given two lines above in the hypotheses. Harmless duplication; could be trimmed.
3. The word "Coverage" in the title of Proposition 3.21 (prop:robustness-coverage) may invite a referee to expect a 1-alpha frequentist coverage statement. Consider retitling "Asymptotic bias bound for psi(theta-dagger)" to preempt the misread, even though the body is correct.

## Bottom line (logic)

No critical or major logic defects. Both prior Criticals are resolved at the proof level. The new robustness-interval Definition 3.20 and Proposition 3.21 are correct and submission-grade, with one optional clarifying sentence (existence of the uniform Hessian bound via non-singular information) and one optional retitle. The C2 core (exact preservation, m=2 C1 expansion, general-m preservation, C3 linear system) is sound and consistent with the foundational C2 definition.
