# Citation Verifier Report

**Date**: 2026-06-04
**Focus**: citation accuracy, missing references, bibliography integrity. Prior-art freshness checked live via DOI resolution + Crossref/arXiv APIs.

## Bibliography integrity

- 35 entries defined in refs.bib; 34 cited in the included sections. No broken citations (every \cite key resolves to a bib entry). Bibtex log is clean (no errors/warnings).
- **One dead entry (MINOR):** `Stefanski2002` (Stefanski & Boos, "The Calculus of M-Estimation," Am. Statistician 2002) is defined but never cited. Either cite it (it is a natural reference for the M-estimation framing of the pseudo-true parameter and the sandwich variance mentioned in Limitations, sensitivity_framework.tex / discussion.tex) or remove it. Citing it in the Limitations sandwich-variance sentence is the better fix.

## Citation accuracy (spot-checked against resolved metadata)

- **Bakoyannis2025** (the load-bearing concurrent-work citation): arXiv:2511.20980 resolves; title matches the bib entry verbatim ("Robustness intervals for competing risks analysis with causes of failure missing not at random"). The in-text quote at introduction.tex line 62 ("there is no formula for computing the E-value in the context of competing risks analysis with MNAR causes of failure") is consistent with the paper's abstract/scope. Author list (Bakoyannis, Rontogiannis, Zhang, Tu, Mwangi, Yiannoutsos) is plausible and matches the arXiv record.
- **towell2026masked** (foundational companion): DOI 10.5281/zenodo.18725577 resolves; citation_title "Masked Causes of Failure in Series Systems: A Likelihood Framework" matches the bib entry.
- **CinelliHazlett2020**: DOI 10.1111/rssb.12348 resolves to JRSS-B 82(1):39. Matches.
- **Huairu-2013** (the application data, Guo/Niu/Szidarovszky RAMS 2013): DOI 10.1109/rams.2013.6517765, consistent.
- Foundational reliability chain (Usher-1988 10.1109/24.9880; Lin-1993 10.1109/24.273596; LinGuess1994; Guttman1995; Mukhopadhyay2006; CraiuReiser2006; craiu2004inference Biometrika 91(3):543; flehinger2002parametric LDA 2002): all DOIs well-formed and consistent with known records.
- ISNI / sensitivity chain (TroxelMaHeitjan2004 Statistica Sinica 14(4):1221 [no DOI, correct: that volume predates DOI assignment for the journal]; ZhangHeitjan2006 10.1111/j.1541-0420.2006.00582.x; CopasEguchi2005; SiannisCopasLu2005; HeitjanRubin1991; JacobsenKeiding1995; Rubin1976): all consistent.

No citation appears to be a hallucinated or misattributed reference.

## MIN-3 from prior review (add towell2026masked): RESOLVED

The 2026-05-27 review asked for the foundational paper to be cited (it was previously citing only the thesis and Lin-1993). towell2026masked is now cited in introduction.tex line 12 and background.tex Theorem (thm:like-C123). The foundational citation chain across the research series is now in place. (See novelty-assessor/prose-auditor for the separate point that the boundary should also be stated in prose, not just co-cited.)

## Prior-art freshness (live searches)

WebSearch was not available; Crossref + arXiv + DOI resolution via curl WERE, and were used.

**Recommended add (MINOR):** Hu, Huang, Shen (2023), "Maintenance optimization of a two-component series system considering masked causes of failure," Quality and Reliability Engineering International, 10.1002/qre.3423. This is the most recent reliability-side masked-cause paper and is currently absent. It is a maintenance/decision paper (does not threaten novelty), but its inclusion demonstrates the reliability lineage is current through 2023. The paper's reliability prior-work paragraph (background.tex lines 142-155) currently stops at 2006 (Mukhopadhyay/Craiu-Reiser) plus Flehinger 2002; one recent citation would close that ~20-year gap.

**Optional add:** Lo, Ma, Manuguerra, Moreno-Betancur (2022), SMMR, 10.1177/09622802211070254, MAR penalized cause-specific Cox. Same family as the already-cited MorenoBetancur2015; only worth adding if the missing-cause biostatistics paragraph is expanded.

**Not needed:** the 2025/2026 interval-censored-competing-risks-with-missing-types papers (CSDA 10.1016/j.csda.2025.108229; LDA 10.1007/s10985-026-09698-x; SMMR 10.1177/09622802261420820). All are semiparametric, single-label, interval-censored; they fall on the point-cause side of the paper's own non-nesting Lemma C.2 and do not need to be engaged individually.

## Bottom line (citations)

Bibliography is clean and accurate; the key concurrent-work and foundational citations resolve and match. One dead entry to cite-or-cut (Stefanski2002), one recommended currency add (qre.3423). No missing must-cite that threatens the contribution.
