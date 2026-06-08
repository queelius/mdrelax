# Citation Verifier Report

**Date**: 2026-06-08
**Focus**: citation accuracy, missing references, bibliography integrity.

## Bibliography integrity (re-checked this pass)

- refs.bib parses cleanly; bibtex log (`paper/main.blg`) has no errors or warnings. The clean build resolves every `\cite` (0 undefined citations in main.log excluding font-shape lines).
- **One dead entry (MINOR, carryover, still present):** `Stefanski2002` (Stefanski & Boos, "The Calculus of M-Estimation," Am. Statistician 2002) is defined in refs.bib but cited nowhere in the included sections (verified by grep across paper/sections/ and main.tex: no match). Either cite it (it is the natural reference for the M-estimation framing of the pseudo-true parameter and the sandwich variance named in Limitations) or remove it. Citing it in the Limitations sandwich-variance sentence is the cleaner fix.

## Citation accuracy (spot-checked against known metadata + live DOI)

- **Bakoyannis2025** (load-bearing concurrent work): arXiv:2511.20980; refs.bib title matches verbatim. The in-text quote (introduction.tex line 76-79: "there is no formula for computing the E-value in the context of competing risks analysis with MNAR causes of failure") is consistent with the paper's scope. Author list plausible.
- **towell2026masked** (foundational): DOI 10.5281/zenodo.18725577; title "Masked Causes of Failure in Series Systems: A Likelihood Framework" matches.
- **Huairu-2013** (application data, Guo/Niu/Szidarovszky RAMS 2013): DOI 10.1109/rams.2013.6517765, consistent. NOTE the bibkey is `Huairu-2013` but the first author surname is Guo (Huairui is the given name); the in-text prose correctly says "Guo et al." This is a slightly misleading bibkey but not an error in the rendered citation. Cosmetic.
- Foundational reliability chain (Usher-1988, Lin-1993, LinGuess1994, Guttman1995, Mukhopadhyay2006, CraiuReiser2006, craiu2004inference, flehinger2002parametric): DOIs well-formed and consistent.
- ISNI / sensitivity chain (TroxelMaHeitjan2004 [no DOI, correct: that Statistica Sinica volume predates DOI assignment], ZhangHeitjan2006, CopasEguchi2005, SiannisCopasLu2005, HeitjanRubin1991, JacobsenKeiding1995, Rubin1976, HarelSchafer2009, XieGaoHeitjan2018, CinelliHazlett2020): all consistent.
- Misclassified-cause chain (Ebrahimi1996, VanRompaye2012, HaTsodikov2015, MorenoBetancur2015, BakoyannisYiannoutsos2015): all consistent.

No hallucinated or misattributed reference found.

## Internal citation/metadata consistency note (shared with format + prose)

The refs.bib HEADER comment (lines 1-3) still titles the paper "Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption", the OLD framing. The on-page title, .zenodo.json, and CITATION.cff use the new "Coarsening Conditions..." framing. This is a comment, so it does not affect any rendered citation, but it should be updated when the title is reconciled (see review.md). The `towell2026mdrelax` citation key that sibling papers use is unaffected.

## Prior-art freshness / currency adds (carryover, optional)

- **Recommended (MINOR):** Hu, Huang, Shen (2023), "Maintenance optimization of a two-component series system considering masked causes of failure," QREI, 10.1002/qre.3423. Most recent reliability-side masked-cause paper; a maintenance/decision paper (does not threaten novelty). The reliability prior-work paragraph (background.tex lines 142-155) currently stops at 2006 + Flehinger 2002; one recent cite closes the gap.
- **Optional:** Lo, Ma, Manuguerra, Moreno-Betancur (2022), SMMR, 10.1177/09622802211070254 (same family as the cited MorenoBetancur2015; only if that paragraph is expanded).
- **Not needed:** 2025/2026 interval-censored-competing-risks-with-missing-types papers; all semiparametric, single-label, on the point-cause side of the paper's own non-nesting Lemma C.2.

## Bottom line (citations)

Bibliography is clean and accurate; the concurrent-work and foundational citations resolve and match. Items: cite-or-cut the dead Stefanski2002 entry, update the stale refs.bib header title when the title is reconciled, and (optional) add the qre.3423 currency cite. No missing must-cite threatens the contribution.
