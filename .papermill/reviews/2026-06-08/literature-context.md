# Literature Context Packet

**Date**: 2026-06-08
**Paper**: Sensitivity of Series System Reliability Estimation to the Coarsening Conditions on Masked Failure Data (towell2026mdrelax)
**Method note**: This is a near-final paper with an exhaustive prior-art reassessment (2026-05-04, three parallel scouts) and a current literature packet from the 2026-06-04 review. The present pass re-confirms freshness and identity rather than re-surveying from scratch; live DOI/Zenodo checks were run.

## Identity / DOI verification (live, 2026-06-08)

| Item | DOI / ID | Status (this pass) |
|---|---|---|
| This paper (concept) | 10.5281/zenodo.20414727 | Resolves -> latest version record 20468529 |
| This paper (latest published version) | 10.5281/zenodo.20468529 | Resolves; **published title is the OLD "...non-informative masking assumption"** |
| Foundational companion (towell2026masked) | 10.5281/zenodo.18725577 | Cited; in refs.bib |
| Bakoyannis 2025 (the work extended) | arXiv:2511.20980 | Cited; title matches refs.bib verbatim |

The live Zenodo title check is the load-bearing freshness fact this pass: the published DOI landing page still reads "Sensitivity of series system reliability estimation to the non-informative masking assumption", whereas the on-page PDF title, .zenodo.json, and CITATION.cff now read "...Coarsening Conditions on Masked Failure Data". See the format/prose findings in review.md.

## Three lineages the paper sits in (stable since 2026-06-04)

1. **Masked cause-of-failure reliability (home literature).** Usher-Hodgson 1988, Usher-Lin-Guess 1993, Lin-Guess 1994, Guttman et al. 1995, Flehinger et al. 2002, Craiu-Duchesne 2004, Mukhopadhyay-Basu 2006, Craiu-Reiser 2006, Guo/Niu/Szidarovszky 2013 (application data), towell2023reliability, towell2026masked (foundational), towell2026binary (k-out-of-m companion). All cited.

2. **Sensitivity to nonignorability / local misspecification (method ancestry).** Heitjan-Rubin 1991 (CAR), Jacobsen-Keiding 1995, Rubin 1976 (parameter distinctness), Troxel-Ma-Heitjan 2004 (ISNI, direct ancestor), Zhang-Heitjan 2006 (closest construction), Copas-Eguchi 2005, Siannis-Copas-Lu 2005, Cinelli-Hazlett 2020, VanderWeele-Ding 2017, White 1982 / Huber 1967. All cited. Stefanski-Boos 2002 (M-estimation calculus) is in refs.bib but uncited (see citation findings).

3. **Missing / misclassified cause of failure (the non-nested neighbor).** Ebrahimi 1996, Van Rompaye 2012, Ha-Tsodikov 2015, Moreno-Betancur 2015, Bakoyannis-Yiannoutsos 2015, Bakoyannis et al. 2025 (concurrent robustness intervals). The paper's Appendix C non-nesting lemma (lem:non-nesting) establishes this lineage is a neighbor, not a competitor.

## Freshness verdict (unchanged, re-confirmed)

The prior-art framing remains current and defensible. The single largest "are they scooped?" risk, the concurrent Bakoyannis 2025 robustness-interval paper, is cited, characterized as concurrent, and differentiated on two axes (semiparametric Cox binary-missing vs parametric set-valued masking). No 2024-2026 work supersedes the contribution.

Currency adds carried from prior passes (both optional, neither a must-cite, neither threatens novelty):
- Hu, Huang, Shen (2023), QREI, 10.1002/qre.3423 (most recent reliability-side masked-cause; closes a roughly two-decade gap in the reliability prior-work paragraph, which currently stops at 2006 + Flehinger 2002).
- Lo, Ma, Manuguerra, Moreno-Betancur (2022), SMMR, 10.1177/09622802211070254 (same family as the already-cited MorenoBetancur2015; only if that paragraph is expanded).

## Boundary against the foundational companion (towell2026masked)

Real and distinct: the foundational paper establishes the C1-C2-C3 likelihood; this paper quantifies the bias when those conditions fail (pseudo-true expansions, ISNI specialization, robustness intervals, exact-preservation theorem, robustness hierarchy, the non-nesting lemma). None of those appear in the foundational paper, so there is no double-counting. The manuscript still co-cites towell2026masked only for the likelihood and never states the boundary in one sentence (see novelty findings).
