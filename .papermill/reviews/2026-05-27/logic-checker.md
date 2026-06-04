# Logic and Proof Correctness Review

**Date:** 2026-05-27
**Paper:** "Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption"

## Summary

The paper's core theoretical results are correct but unevenly stated. The strongest result (Theorem 3.3, exact total-hazard preservation for exponential) is rigorous and the proof is clean. The C1 misspecification result (Theorem 3.5, m=2 case) is correct but the proof has a small but important gap in the "combining them yields a + br = 0" step. Theorem 3.7 (C3 misspecification) is presented essentially without proof and the contribution claims rest on a result whose details are deferred to nonexistent supplementary material. Several inference steps in the connecting prose also overreach what the formal results support.

## Critical Findings

### C-1: Theorem 3.7 (Pseudo-True Parameter Under C3 Misspecification) is asserted without a complete proof
- **Location:** sensitivity_framework.tex, lines 524-548 (Theorem 3.7).
- **Quoted text:** "Closed-form expressions for the individual coefficients $\Delta_j^{C_3}$ with $m = 2$ and general $m$ are given in the supplementary material."
- **Problem:** No supplementary material exists in the repository. The theorem statement asserts a first-order expansion with $\Delta_j^{C_3}$ "obtained as the solution of a linear system analogous to (c1-linsys-1, c1-linsys-2) but with right-hand-side terms determined by the derivatives of $p_j(k; \vtheta)$ with respect to $\vtheta$". This is a plausibility argument, not a proof. The construction is also asymmetric (paper says "with $p_k(k) = 1$ ... and the same functional form for all $k$ (so C2 also holds in a restricted sense)", line 517). The phrase "in a restricted sense" is undefined.
- **Suggestion:** Either provide the explicit derivation as Appendix B (alongside the existing score functions appendix), or weaken the contribution claim. The simulations in Section 4 only validate first-order preservation empirically via the slope test, which is consistent with a stated result but does not constitute a proof.
- **Cross-verified:** see also methodology-auditor finding M-1.

### C-2: The C2 contribution mismatch is structural, not just a phrasing issue
- **Location:** introduction.tex contribution 2, lines 80-87; sensitivity_framework.tex Theorem 3.4 (misspec-C2) and Proposition 3.4.
- **Quoted text (intro):** "We derive closed-form first-order expansions of the misspecified MLE in each $\alpha_{C_k}$ at the standard-model fit, specializing the ISNI construction... The expansions are explicit for exponential components and tractable for Weibull components."
- **Problem:** For C2, the only "expansion" in the main text is Proposition 3.4, which gives a sufficient condition for total-hazard robustness, not the first-order expansion of $\theta_j^\dagger$. Theorem 3.4 only invokes White-Huber to assert existence of a pseudo-true parameter; the closed-form first-order expansion of the individual $\theta_j^\dagger$ in $s$ is not given. For Weibull under C2, no expansion at all is given; only the empirical slope test in Section 4. The contribution claim is overstated relative to what is delivered.
- **Suggestion:** Either deliver the explicit first-order expansions for the exponential C2 case (analogous to Theorem 3.5 for C1), or rewrite the contribution to say "first-order expansions in closed form for the exponential case for C1 violations, with the exponential C2 case characterized via Proposition 3.4 and the C3 and Weibull cases characterized empirically".

## Major Findings

### M-1: Theorem 3.5 (Pseudo-True Parameter Under C1 Misspecification, m=2) proof has a gap in the combining step
- **Location:** sensitivity_framework.tex lines 449-458 (proof of Theorem 3.5).
- **Quoted text:** "Setting $a = \Delta_1 / \theta_1^*$, $b = \Delta_2 / \theta_2^*$, and $r = \theta_2^* / \theta_1^*$, the system (eq:c1-linsys-1), (eq:c1-linsys-2) reduces to two linear equations in $(a, b)$. Combining them yields $a + br = 0$, i.e., $\Delta_1 + \Delta_2 = 0$."
- **Problem:** "Combining them yields $a + br = 0$" hides the actual combination. Direct addition of (c1-linsys-1) and (c1-linsys-2) does not yield $a + br = 0$; it yields $p_0(\theta_2^* - \theta_1^*)^2/(\theta_1^*\theta_2^*)$ on the RHS, which is generally nonzero. The correct combination is $\theta_1^* \cdot$ (c1-linsys-1) $+ \theta_2^* \cdot$ (c1-linsys-2), which after rearrangement yields $(1-p_0)(\Delta_1 + \Delta_2) + p_0(\Delta_1 + \Delta_2) = 0$, hence $\Delta_1 + \Delta_2 = 0$. The result is correct but the proof step is underspecified.
- **Suggestion:** Rewrite the proof step explicitly: "Multiplying (c1-linsys-1) by $\theta_1^*$, (c1-linsys-2) by $\theta_2^*$, and adding the two resulting equations cancels the right-hand sides and yields $(1-p_0+p_0)(\Delta_1 + \Delta_2) = \Delta_1 + \Delta_2 = 0$. Substituting $\Delta_2 = -\Delta_1$ into (c1-linsys-1) gives $\Delta_1 = p_0(\theta_2^* - \theta_1^*)/(1-p_0)$."

### M-2: The "first-order expansion" rhetoric overreaches in the C1 corollary
- **Location:** sensitivity_framework.tex lines 461-481 (Corollary 3.6, c1-hazard-general).
- **Quoted text:** "In particular, the $m=2$ first-order coefficients in (eq:c1-delta) satisfy $\Delta_1^{C_1} + \Delta_2^{C_1} = 0$, and this preservation holds at all orders for general $m$."
- **Problem:** The corollary proof appeals to Theorem 3.3 (exact-hazard-exp) which requires the retained system lifetimes to satisfy $T \sim \mathrm{Exp}(S^*)$. The proof correctly verifies that retention probability under the Bernoulli C1 model is $1 - \alpha_{C_1}(1-p_0)^{m-1}$, independent of $t$ and $k$. Independence is correctly argued. However, the leap to "preservation holds at all orders for general $m$" is a consequence of the exact-preservation theorem, not of the m=2 first-order expansion in (c1-delta). The phrasing reads as though the first-order coefficients themselves extend to higher orders, when really the exact preservation makes the higher-order coefficients all vanish for total hazard. The result is correct but the prose conflates two distinct results.
- **Suggestion:** Restate: "By Theorem 3.3 applied to the Bernoulli C1 violation, the total system hazard is preserved exactly at all orders for any $m \geq 2$, not merely at first order; this generalization is a corollary of the profile-likelihood factorization argument rather than of the m=2 first-order expansion in (c1-delta)."

### M-3: The "first-order analysis for Weibull" claim is empirical, not analytical
- **Location:** sensitivity_framework.tex lines 612-628 (discussion of robustness hierarchy item 3 asymmetry).
- **Quoted text:** "Under C2 and C3 violations, the candidate-set bias remains balanced enough that it does not propagate to the time portion at first order. Under C1 violation, the candidate set sometimes omits the true cause entirely, producing more aggressive candidate-set bias that does propagate to the time-portion estimates at first order through the shared parameterization."
- **Problem:** This is a verbal explanation of the asymmetry. The actual evidence cited is the simulation slope tests in Table 2 (Section 4). The paper does not derive the Weibull first-order coefficients analytically. The "structural origin" prose suggests an analytical result, but the paper provides only an empirical one. A reader expecting a Weibull score-equation analysis will be disappointed.
- **Suggestion:** Either provide the Weibull-specific first-order analysis (the calculation is more involved than the exponential because of the shape parameter dependence in the hazard), or label the verbal explanation as a heuristic interpretation of the empirical pattern: "We conjecture that the structural reason for the C1 versus C2/C3 asymmetry is that..." rather than asserting it.

### M-4: Proposition 3.4 (System Hazard Robustness) overstates the "moderately violated" claim
- **Location:** sensitivity_framework.tex lines 164-177 (Proposition 3.4 and its proof).
- **Quoted text:** "When C2 is moderately violated, (eq:hazard-condition) holds approximately and the total hazard bias is small."
- **Problem:** The condition (eq:hazard-condition) is stated as an equality, but the proposition's prose calls "moderate" violation a case where the condition "holds approximately." No formal quantification of "approximately" or "moderate" is given. The supporting evidence is the simulation observation that "the total hazard bias remains below 3% across all severity levels," but that 3% bound is empirical for a specific 5-component Weibull system, not a theoretical bound. The Bakoyannis 2025 motivation mentioned in the intro emphasizes the need for *bounded* sensitivity intervals (E-values, robustness values), which this proposition does not provide.
- **Suggestion:** Either rephrase to "in the simulation studies of Section 4, the total hazard bias remains below 3% across all severity levels" (empirical observation, not a theoretical result), or strengthen the proposition by providing a bound in terms of $\|\pi_{k,c} - \bar\pi_c\|$. The exponential case is exactly handled by Theorem 3.3, so the "approximate" language is most pressing for the Weibull case, which the proposition does not explicitly cover.

### M-5: Theorem 3.4 (misspec-C2) does not actually establish the bias structure it claims
- **Location:** sensitivity_framework.tex lines 125-162.
- **Quoted text:** "the pseudo-true parameter differs from $\vtheta^*$ whenever the masking weights $\pi_{k,c}$ vary with $k$, i.e., whenever C2 is violated."
- **Problem:** The proof derives the misspecified score and notes that "the misspecified score replaces the weights $\pi_{k,c}$ with uniform weights. Setting $\E[\partial\ell^{wrong}/\partial\theta_j] = 0$ at $\vtheta^\dagger$ yields a system whose solution differs from $\vtheta^*$ whenever the masking weights vary with $k$." This is asserted, not proven. A trivial counterexample to the "differs" claim: if $|c| = 1$ for all observed candidate sets, then the score under both true and misspecified models is identical (since the sum collapses to a single term), and $\vtheta^\dagger = \vtheta^*$ regardless of $\pi_{k,c}$ variation. The "whenever" claim is too strong. The correct statement should be: "the pseudo-true parameter may differ from $\vtheta^*$ when the masking weights vary with $k$, and generically does."
- **Suggestion:** Strengthen the proof by either (a) showing that the score equation system is non-degenerate in a precise sense (e.g., when at least one observation has a multi-element candidate set with informative weights), or (b) softening the claim to "generically differs from $\vtheta^*$" with a footnote noting the singleton-candidate-set degenerate case.

## Minor Findings

### m-1: Theorem 3.3 proof uses retention condition but its statement is buried
- **Location:** sensitivity_framework.tex Theorem 3.3 statement (lines 218-231) vs Remark 3.4 scope (lines 264-276).
- **Problem:** The theorem statement says "suppose data are generated under any masking model such that retained system lifetimes satisfy $T \sim \mathrm{Exp}(S^*)$". This is a strong condition that the reader must verify in each application. The remark below the proof discusses scope, but a reader skimming the theorem statement might miss that this is the key requirement. The corollaries (3.5 for C2, 3.6 for C1) verify this condition for their specific Bernoulli models.
- **Suggestion:** Add an explicit name to this condition (e.g., "Independent Retention Condition") and reference it by name in the corollaries, so the structural argument is more visible.

### m-2: Definition of "robustness interval" is not formalized in the paper
- **Location:** introduction.tex contribution 3 (lines 89-96); discussion.tex (lines 184-191); conclusion.tex contribution 3 (lines 33-40).
- **Problem:** The paper repeatedly promises "robustness intervals" as a contribution, but the formal definition is never given. The discussion section references "the robustness intervals of Section 3.7" but Section 3.7 only defines ISNI (Definition 3.10) and states the robustness hierarchy corollary; no robustness interval construction is given. The conclusion claims the implementation is in the mdrelax R package, but the appendix software description does not mention robustness-interval functions.
- **Suggestion:** Either formally define a robustness interval (e.g., $RI_\alpha = [\hat\theta - \alpha \cdot |\nabla_\alpha \theta^\dagger|, \hat\theta + \alpha \cdot |\nabla_\alpha \theta^\dagger|]$ at calibrated severity $\alpha$) and demonstrate its construction in the application, or remove "robustness intervals" as a numbered contribution and reframe as "robustness hierarchy informed by ISNI". Currently the contribution is asserted but not delivered in the main text.

### m-3: The empty-candidate-set retention argument has a subtle issue
- **Location:** sensitivity_framework.tex lines 387-393 and Corollary 3.6 proof.
- **Quoted text:** "Failures with no diagnostic information are excluded from the parametric MLE."
- **Problem:** Excluding observations with empty candidate sets introduces a sample-size reduction that the C1-C2-C3 likelihood does not account for. The corollary proof argues that retention is independent of $T$, so the conditional distribution of $T$ given retention equals the marginal. But for the likelihood to be valid, the practitioner workflow of "drop the empty candidate sets and fit C1-C2-C3" treats the remaining observations as if they were a random sample from the standard model. Since the retention probability is independent of $T$ and $K$, this is asymptotically valid for parameter estimation, but the effective sample size is reduced. This is implicit in the analysis but not made explicit; a careful reviewer will note that the variance estimates from the standard FIM are computed on the retained sample, which is what one would do, but the connection to the original sampling design (where some observations have missing data) is not discussed.
- **Suggestion:** Add a sentence: "The retained sample is a random thinning of the original sample independent of $T$ and $K$, so the C1-C2-C3 MLE on the retained sample is consistent for $S^*$ with effective sample size $n \cdot [1 - \alpha_{C_1}(1-p_0)^{m-1}]$ at the exponential / total hazard level. Variance estimates from the standard FIM correspondingly reflect this reduced effective sample size."

### m-4: The non-identifiability theorem (3.8) proof is presented as constructive for m=2 only
- **Location:** sensitivity_framework.tex lines 297-344.
- **Quoted text:** "We prove the result constructively for $m = 2$ components; the general case follows analogously."
- **Problem:** "Follows analogously" is doing a lot of work. The constructed deformation in m=2 explicitly preserves three functionals (system hazard, weighted candidate-set hazard for $c=\{1,2\}$, single-component hazards for $c=\{1\}$, $c=\{2\}$). For general $m$, the number of possible candidate sets grows exponentially ($2^m - 1$), and the constraint structure is more involved. The claim of analogous proof is plausible but not trivial. The non-identifiability claim is also supported by the simulation in the appendix (lines 47-69 of appendix.tex), which is empirical, not analytical.
- **Suggestion:** Either provide the general-$m$ proof in the appendix (counting argument: the parameter space for $(\theta, P)$ has dimension $m + m(m-1)$, the candidate-set distribution has at most $2^m - 1$ degrees of freedom plus 1 for $S^*$, so identifiability requires $m + m(m-1) \leq 2^m - 1 + 1$ which fails for small $m$), or restrict the theorem to m=2 with a remark on the general case being open.

### m-5: Lemma 3.1 (set-valued vs point-cause non-nesting) proof direction (ii) requires m=4
- **Location:** appendix.tex Lemma 3.1 proof (lines 113-139).
- **Quoted text:** "For $m = 4$ and any $K$, take..."
- **Problem:** The proof of direction (ii) uses m=4 explicitly. The lemma statement claims non-nesting "for $m \geq 2$", but direction (ii) is only proven for $m \geq 4$. The two directions could be presented as separate sub-lemmas with their applicable $m$ ranges, or the proof could be generalized to $m = 2, 3$.
- **Suggestion:** Either generalize the direction (ii) proof to all $m \geq 2$ (likely possible: take $\pi_{2|1} = 1$ with $\rho(\{2, 3\}|2) = 1/2$, $\rho(\{2\}|2) = 1/2$ for $m = 3$), or restate the lemma as "for $m \geq 4$" with $m = 2$ and $m = 3$ noted as boundary cases.

## Suggestions

1. Make explicit in Section 3.1 that "first-order expansion" terminology is used in two distinct senses in this paper: (a) the analytic Taylor expansion of $\theta^\dagger$ in $\alpha$ at $\alpha = 0$ (the C1 m=2 exponential case in Theorem 3.5), and (b) the empirical slope-test verdict that the linear coefficient of bias-vs-alpha is statistically indistinguishable from zero (the Weibull C2 and C3 cases in Section 4). These are different things and conflating them weakens the contribution narrative.
2. Move the application section's m=3 turbine engine analysis to use the C1, C2, AND C3 sweeps (currently only C2 is swept), to demonstrate the full unified framework on real data. The current application is essentially a C2-only analysis.
3. The "robustness hierarchy" corollary (3.11) item 3 cites a specific empirical bound ("below 4% in the 5-component sweep"). This is appropriate but should be flagged as system-specific rather than universal: a 5-component system with shape and scale heterogeneity matching the chosen design may behave differently from a 10-component homogeneous system. The discussion section already lists this as a limitation; it could be cross-referenced from the corollary.
