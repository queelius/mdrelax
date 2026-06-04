# Literature Context: Sensitivity Analysis for Masked Series-System Reliability

**Date:** 2026-05-27
**Paper:** "Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption"

## Direct Prior Art (Same Problem)

### Reliability literature: dependent masking
- **Usher and Hodgson 1988** (cited): Foundation of C1-C2-C3 framework. Original masked data MLE for series systems.
- **Usher, Lin, Guess 1993** (cited as Lin-1993): Exact MLE for masked exponential/Weibull. Closest predecessor for the standard estimator the paper studies.
- **Lin and Guess 1994** (cited): Proportional dependent masking, m=2. Closest in intent to the C2 relaxation approach: parameterizes departure from C2 with a fixed ratio.
- **Guttman, Lin, Reiser, Usher 1995** (cited): Bayesian dependent masking, m=2.
- **Mukhopadhyay and Basu 2006** (cited): General m Bayesian, source of the power-weight C3 violation form. The paper's C3 sweep parameterization comes from here.
- **Craiu and Reiser 2006** (cited): Conditional masking probability models; identifiability under structural assumptions.
- **Craiu and Duchesne 2004** (cited): EM-based inference for masked causes.
- **Flehinger, Reiser, Yashchin 2002** (cited): Parametric Weibull masked competing risks with two-stage diagnostic procedures.

**Verdict:** The paper engages with the direct reliability lineage well. No obviously missing references.

### Biostatistics: missing/misclassified cause of failure
- **Ebrahimi 1996** (cited): Earliest cause-of-failure misclassification likelihood. Point-cause (not set-valued); paper's appendix lemma demonstrates non-nesting.
- **Van Rompaye, Jaffar, Goetghebeur 2012** (cited): Cox-model misclassification sensitivity. Point-cause.
- **Ha and Tsodikov 2015** (cited): Semiparametric proportional hazards with misclassified cause. Point-cause.
- **Moreno-Betancur, Rey, Latouche 2015** (cited): Pattern-mixture sensitivity analysis for competing risks with missing causes. Self-described as ad hoc GOF.
- **Bakoyannis and Yiannoutsos 2015** (cited): Closed-form bias under cumulative-incidence misclassification.
- **Bakoyannis et al. 2025** arXiv:2511.20980 (cited as concurrent work): Robustness intervals for Cox cause-specific hazards under MNAR cause. Submitted Nov 2025; paper claims this is concurrent with theirs.

**Verdict:** This corpus is well-engaged with. The set-valued vs point-cause distinction is supported by a non-nesting lemma in the appendix.

## Methodological Foundations (Sensitivity Analysis)

### Coarsening-at-Random and Ignorability
- **Rubin 1976** (cited): Inference and missing data, parameter distinctness condition.
- **Heitjan and Rubin 1991** (cited): Coarsening at random. C1+C2+C3 in the paper's framing is equivalent to CAR for the candidate-set coarsening map.
- **Jacobsen and Keiding 1995** (cited): CAR in general sample spaces.

**Verdict:** Engaged correctly. The paper appropriately presents C1-C2-C3 as a specialization of CAR.

### Local Sensitivity Indices
- **Troxel, Ma, Heitjan 2004** (cited): ISNI. Direct methodological ancestor of paper's first-order expansions.
- **Zhang and Heitjan 2006** (cited): ISNI for general coarse-data model. Closest prior local-sensitivity construction; paper specializes from interval-censoring to set-valued candidate-set coarsening.
- **Copas and Eguchi 2005** (cited): Local model uncertainty and first-order bias under nuisance perturbation.
- **Siannis, Copas, Lu 2005** (cited): Sensitivity for informative censoring in parametric survival.
- **Siannis 2004** (cited): Applications of parametric model for informative censoring.

**Verdict:** Strong engagement. The ISNI tradition is fully credited.

### Calibrated Sensitivity Bounds
- **Cinelli and Hazlett 2020** (cited): Omitted-variable robustness value. Stylistic anchor for severity-scale framing.
- **VanderWeele and Ding 2017** (cited): E-value for observational research.
- **Rosenbaum 2002** (cited): Observational studies.

**Verdict:** Engaged. The "calibrated severity" framing is positioned as an extension to multi-dimensional sensitivity geometry.

### Misspecification Theory
- **White 1982**, **Huber 1967** (cited): Foundation for pseudo-true MLE convergence.
- **Stefanski and Boos 2002** (cited): M-estimation calculus; cited but not load-bearing.

**Verdict:** Foundational references present.

### Identifiability and Partial Identification
- **Tsiatis 1975** (cited): Non-identifiability of competing risks.
- **Crowder 2001** (cited): Classical competing risks book.
- **Manski 2003** (cited): Partial identification of probability distributions.

**Verdict:** Standard non-identifiability literature engaged.

## Potential Gaps and Missing References

### Possibly missing, Cox / proportional hazards extensions
- **Bordes et al.** literature on identifiability in mixture and competing risks settings could be added if positioning is desired.
- **Tsiatis (1975)** and **Heckman and Honore (1989)** on competing-risks identifiability. Only Tsiatis is cited.
- **Bickel et al.** efficient semiparametric estimation theory. Not cited, but plausibly not needed for parametric setting.

### Possibly missing, Information-theoretic / Fisher information for missing-data
- **Louis 1982** (observed-information formula for incomplete data). Not cited but standard reference; could strengthen the discussion of FIM/coverage.
- **Orchard and Woodbury 1972** (missing information principle). Not cited.

### Possibly missing, Reliability identifiability under masked data
- **Reiser et al.** (multiple papers). Guttman and Reiser are cited; others may be relevant.
- **Sarhan 2003, 2005** on masked failure cause Bayesian analysis. Not cited; possibly relevant.

### Possibly missing, Software/reproducibility comparison points
- The contribution "no related work releases code" should be checked against:
  - `riskRegression` R package (Gerds et al.). Covers competing risks regression but not the masked-set-valued case the paper claims novelty in.
  - `survival` package. Standard for survival analysis.
  - These are tangential; the contribution claim "none of the closest related work releases code" is defensible if scoped to direct sensitivity-analysis methodology papers.

### Concurrent work concern
- **Bakoyannis et al. 2025** arXiv:2511.20980. Paper labels this as concurrent work. The differentiation argument (paper covers set-valued candidate sets vs their binary missing cause; parametric series systems vs semiparametric Cox) is reasonable. The non-nesting lemma in the appendix supports the structural distinctness claim.

## Bibliography Adequacy Assessment

The paper's bibliography (33 entries, 14 added in the May 2026 reframing) is well-targeted to the reframed thesis. It engages with:
- The classical reliability dependent-masking literature (Lin and Guess, Guttman, Mukhopadhyay, Craiu, Reiser, Flehinger).
- The CAR/missing-data foundation (Rubin, Heitjan-Rubin, Jacobsen-Keiding).
- The ISNI methodological lineage (Troxel-Heitjan, Zhang-Heitjan).
- The local-bias / sensitivity bound corpus (Copas-Eguchi, Siannis, Cinelli-Hazlett, VanderWeele-Ding).
- The biostatistics missing-cause analogues (Ebrahimi, Van Rompaye, Ha-Tsodikov, Moreno-Betancur, Bakoyannis x2).
- The misspecification theory (White, Huber).
- Identifiability foundations (Tsiatis, Crowder, Manski).

No major gaps are evident relative to the paper's positioned scope. Minor possible additions for completeness (Louis 1982 for FIM, Sarhan for related Bayesian masked-data work) would strengthen but are not load-bearing.
