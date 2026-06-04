# Literature Context Packet (merged scouts)

**Date**: 2026-06-04
**Paper**: Sensitivity of Series System Reliability Estimation to the Coarsening Conditions on Masked Failure Data (towell2026mdrelax)
**Tooling note**: WebSearch proper was not available, but `curl` + DOI resolution + the Crossref and arXiv REST APIs WERE reachable and were used for prior-art freshness. Publisher landing pages (Wiley, SAGE) are Cloudflare-gated; clean metadata was pulled from the Crossref API instead.

## DOI / identity verification (resolved live)

| Item | DOI / ID | Status |
|---|---|---|
| This paper (concept) | 10.5281/zenodo.20414727 | Resolves; citation_title matches manuscript title |
| Foundational companion (towell2026masked) | 10.5281/zenodo.18725577 | Resolves; "Masked Causes of Failure in Series Systems: A Likelihood Framework" |
| Bakoyannis 2025 (the work the paper extends) | arXiv:2511.20980 | Resolves; title matches refs.bib exactly: "Robustness intervals for competing risks analysis with causes of failure missing not at random" |
| CinelliHazlett2020 | 10.1111/rssb.12348 | Resolves (JRSS-B 82(1):39) |

## Lineage 1: Masked cause-of-failure reliability (the home literature)

Already engaged and current: Usher-Hodgson 1988, Usher-Lin-Guess 1993, Lin-Guess 1994, Guttman et al. 1995, Flehinger et al. 2002, Craiu-Duchesne 2004, Mukhopadhyay-Basu 2006, Craiu-Reiser 2006, Guo/Niu/Szidarovszky 2013 (the application data), towell2023reliability (thesis), towell2026masked (foundational), towell2026binary (k-out-of-m companion).

Recent reliability-side masked-cause work surfaced via Crossref (2023-2026) NOT cited:
- **Hu, Huang, Shen (2023), Quality and Reliability Engineering International, 10.1002/qre.3423** -- "Maintenance optimization of a two-component series system considering masked causes of failure." Reliability-domain, masked-cause, two-component series. It is a maintenance/decision paper, not a misspecification/sensitivity paper, so it does not threaten novelty, but it is the most recent reliability-side masked-cause article and a reasonable currency citation.
- A 2023 storage-reliability masked-data paper (10.17531/ein/172922) and a 2023 software-hardware embedded masked-data paper (10.1016/j.cie.2023.109746) exist; both are applied estimation, peripheral to a sensitivity contribution.
- CRAN siblings maskedcauses (2026) and maskedhaz (2026) are the author's own ecosystem.

## Lineage 2: Sensitivity to nonignorability / local misspecification (the method ancestry)

Already engaged and current: Heitjan-Rubin 1991 (CAR), Jacobsen-Keiding 1995 (CAR general spaces), Rubin 1976 (parameter distinctness), Troxel-Ma-Heitjan 2004 (ISNI), Zhang-Heitjan 2006 (ISNI for coarse data, the closest ancestor), Copas-Eguchi 2005, Siannis-Copas-Lu 2005, Cinelli-Hazlett 2020 (robustness value), VanderWeele-Ding 2017 (E-value), Stefanski-Boos 2002 (M-estimation calculus), White 1982 / Huber 1967 (misspecified MLE).

## Lineage 3: Missing / misclassified cause of failure (the non-nested neighbor)

Already engaged: Ebrahimi 1996, Van Rompaye et al. 2012, Ha-Tsodikov 2015, Moreno-Betancur et al. 2015, Bakoyannis-Yiannoutsos 2015, Bakoyannis et al. 2025 (concurrent robustness intervals). The paper argues (Appendix C, Lemma C.2 / lem:non-nesting) that set-valued candidate-set masking is non-nested with point-cause misclassification, so this lineage is a neighbor rather than a competitor.

Recent additions to this lineage surfaced via Crossref (2022-2026), all semiparametric biostatistics, NONE of which threaten the positioning because they are point-cause / single-label and interval-censored, not set-valued:
- Lo, Ma, Manuguerra, Moreno-Betancur (2022, SMMR, 10.1177/09622802211070254): MAR penalized cause-specific Cox. Same family as the already-cited MorenoBetancur2015; optional add.
- Lou, Ma, Xiang, Sun (2025, CSDA, 10.1016/j.csda.2025.108229): left-truncated interval-censored competing risks with missing event types.
- Two 2026 papers (LDA 10.1007/s10985-026-09698-x additive-hazards interval-censored missing types; SMMR 10.1177/09622802261420820 missing-cause sensitivity) confirm the topic is active but remain semiparametric/single-label.

## Freshness verdict

The prior-art framing is current and defensible. The concurrent Bakoyannis 2025 robustness-interval paper (the single most important "are they scooped?" risk) is correctly cited, correctly characterized as concurrent, and correctly differentiated (Cox semiparametric binary-missing vs parametric set-valued; the paper even quotes Bakoyannis's "no E-value formula" gap as a target). No 2024-2026 work was found that supersedes the contribution or that is a clear must-cite omission. The one reasonable currency add is Hu-Huang-Shen 2023 (reliability-side masked-cause), and optionally Lo et al. 2022.

## Boundary against the foundational companion (towell2026masked)

The foundational paper establishes the C1-C2-C3 likelihood framework; this paper is the sensitivity/robustness companion that asks what happens when C2 (and C1, C3) are violated. The boundary is real and the contributions are distinct (the foundational paper does not do pseudo-true expansions, ISNI specialization, robustness intervals, or the exact-preservation theorem). However, in the manuscript towell2026masked appears only as a co-citation for the C1-C2-C3 likelihood (background.tex line 89, introduction.tex line 12); there is no sentence explicitly stating "the foundational paper develops X; the present paper develops the sensitivity analysis of Y." Recommended to add one such sentence so the reader and the sibling-paper citation chain see the boundary cleanly.
