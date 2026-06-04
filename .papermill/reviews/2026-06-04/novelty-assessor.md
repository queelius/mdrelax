# Novelty Assessor Report

**Date**: 2026-06-04
**Focus**: contribution clarity, differentiation, significance. Priority: is the contribution distinct from the foundational masked-causes paper and from prior masked-data sensitivity work?

## Contribution list status

The off-by-one from the 2026-05-27 review (MAJ-2) is RESOLVED. Introduction (five enumerated, lines 68-113) and Conclusion (five enumerated, lines 11-56) now agree, and the in-context state.md describes five contributions with the specification test correctly positioned as Future Work (discussion.tex lines 166-180). The five:
1. Calibrated severity scales for C1/C2/C3 as a multi-dimensional sensitivity geometry.
2. Closed-form local sensitivity indices specializing ISNI to set-valued masking (Definition 3.19).
3. First-order robustness intervals (Definition 3.20, Proposition 3.21) extending Bakoyannis 2025 from Cox-semiparametric to parametric series systems.
4. Robustness hierarchy classifying functionals (Corollary 3.17).
5. Reproducible software (mdrelax).

Contribution 3 was the substance of CRIT-2. It is now delivered in the text (def:robustness-interval + prop:robustness-coverage), so the contribution list and the manuscript body finally match. This is the single most important novelty fix since the last review.

## Differentiation from the foundational companion (towell2026masked)

The boundary is real and the contribution is distinct. The foundational paper establishes the C1-C2-C3 likelihood (it is cited for exactly that, background.tex Theorem 3.x and introduction.tex line 12). Nothing in the foundational paper does pseudo-true expansions, ISNI specialization, robustness intervals, the exact-preservation theorem, or the robustness hierarchy. So there is no double-counting of contributions across the two papers.

WEAKNESS (MINOR, framing): the manuscript never states the boundary in a sentence. towell2026masked appears only as a co-citation for the likelihood. For a reader who knows the series exists (and for the sibling coarsening papers that will cite this one), the relationship should be explicit. Add to the introduction, after the contribution list: "Whereas the foundational framework of [towell2026masked] establishes when the C1-C2-C3 likelihood is valid, the present paper quantifies what is lost when those conditions fail." One sentence; closes the only real novelty-framing gap.

## Differentiation from prior masked-data / cause-uncertainty sensitivity work

Strong and now well-defended. Three structural differentiators (introduction.tex lines 48-62):
1. **Set-valued vs single-label.** Every prior cause-uncertainty treatment (Ebrahimi, Van Rompaye, Ha-Tsodikov, Moreno-Betancur, Bakoyannis-Yiannoutsos, Bakoyannis 2025) treats the cause as observed/missing/misclassified to ONE label. The candidate set is a set. This is backed by a proved non-nesting result (Lemma C.2), not just asserted. This is the paper's strongest novelty claim and it holds.
2. **Three orthogonal coarsening conditions vs a scalar departure.** ISNI / E-value / robustness-value frameworks operate on one departure parameter; here C1, C2, C3 carry distinct semantics and are swept independently.
3. **Closed-form parametric bias expressions** where the biostatistics literature is semiparametric. The exponential exact-preservation theorem (3.6) is a genuinely new structural result: it upgrades the usual first-order robustness claim to all-orders for the system hazard, via the profile-likelihood factorization. This is the headline.

Against Bakoyannis 2025 specifically (the concurrent robustness-interval paper, the main scooping risk): the paper cites it as concurrent, extends its idea to a different model class (parametric set-valued vs Cox binary-missing), and even quotes Bakoyannis's own stated gap ("no formula for the E-value in competing risks with MNAR causes") as a target the closed-form expansions fill. This is the right way to handle concurrent work. Verified the arXiv ID and title resolve and match.

## Significance

For the reliability community this is a useful and citable result: the practitioner takeaway ("system-level estimands are robust to all three masking violations, exactly so for exponential; component-level estimands are not") is sharp, actionable, and backed by both theory and simulation. The exact-preservation theorem is the kind of clean structural fact that sibling papers and downstream reliability work can cite directly. The robustness hierarchy (Cor 3.17) is the natural citation anchor.

## Residual novelty-side issues

1. (MINOR, carryover MIN-9) "Multi-dimensional sensitivity geometry" (contribution 1, also conclusion item 1) overstates three INDEPENDENT univariate sweeps. There is no joint/interaction analysis across (alpha_C1, alpha_C2, alpha_C3). Rephrase as "calibrated severities along three independent violation axes," or add a genuine 2-D joint sweep if the geometry claim is to be kept. As written a referee will call this out.
2. (MINOR) Contribution 5 and the conclusion both assert "None of the closest related work releases code." This is a verifiable comparative claim; it is plausible (the cited biostatistics papers do not ship packages) but it is the kind of statement a referee may probe. Keep it, but it is the weakest of the five as a "contribution."
3. (MINOR) The robustness interval (contribution 3) is now defined and proved but never DEMONSTRATED on data (see methodology-auditor). A contribution that is defined but never exercised reads as thinner than one shown in action. Exhibiting one RI on the Guo data would raise the apparent significance of contribution 3 substantially.

## Cross-check requested by area chair (is unclear writing hiding a weak contribution?)

No. The contributions are real and the writing is clear about them. The opposite mild risk applies: contribution 1 ("geometry") is framed slightly more grandly than the actual independent-univariate execution, and contribution 3 is real but under-demonstrated. Neither is a case of weak substance hidden by prose; both are framing/demonstration calibration.

## Bottom line (novelty)

Contributions are distinct from both the foundational companion and the prior point-cause sensitivity literature, with the non-nesting lemma and the exact-preservation theorem as the load-bearing novelty. The contribution-count and robustness-interval-delivery problems from the last review are fixed. Remaining items are framing-grade: add the one-sentence foundational-paper boundary, soften "geometry," and demonstrate the robustness interval on data.
