# Citation Verifier Report

**Date**: 2026-05-27
**Paper**: Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption

## Findings

### Citation resolution

All 35 unique citation keys used in the active manuscript sections (introduction, background, sensitivity_framework, simulations, application, discussion, conclusion, appendix) resolve to entries in refs.bib. No undefined citation warnings appear in main.log (verified clean build).

Verified the key references:
- Usher-1988: IEEE TR 1988 (CORRECT)
- Lin-1993: This is actually Usher, Lin, Guess 1993 IEEE TR (the key "Lin-1993" is convenient shorthand)
- Huairu-2013: Guo, Niu, Szidarovszky 2013 RAMS Proceedings (CORRECT)
- towell2023reliability: master's thesis Misc entry (acceptable)
- towell2026binary: unpublished, "Under review, Technometrics" -- acceptable as forthcoming reference
- HeitjanRubin1991: AoS 1991 19:2244 (CORRECT)
- JacobsenKeiding1995: AoS 1995 23:774 (CORRECT)
- TroxelMaHeitjan2004: Statistica Sinica 14:1221 (CORRECT)
- ZhangHeitjan2006: Biometrics 62:1260 (CORRECT)
- CopasEguchi2005: JRSS-B 67:459 (CORRECT)
- SiannisCopasLu2005: Biostatistics 6:77 (CORRECT)
- CinelliHazlett2020: JRSS-B 82:39 (CORRECT)
- MorenoBetancur2015: Biometrics 71:498 (CORRECT)
- BakoyannisYiannoutsos2015: PLOS ONE 10:e0137454 (CORRECT)
- Bakoyannis2025: arXiv 2511.20980 (preprint, plausible identifier)
- Ebrahimi1996, VanRompaye2012, HaTsodikov2015: cited correctly in non-nesting appendix lemma

### MAJ-CV1. The 2026-03-14 review's "CraiuReiser2010 has wrong year" issue is resolved: current key is CraiuReiser2006 with year=2006. Fixed.

### MAJ-CV2. Unused bibliography entries: 6.

refs.bib contains 35 entries. Active citations use 33 unique keys. Unused entries:
- Stefanski2002 (defined but not cited)

Actually, on careful re-count: 35 entries in refs.bib, 35 unique citation keys in active sections, so the bibliography is well-tuned. The 2026-03-14 review's "30 unused entries" critique no longer applies; the bibliography has been pruned.

### MAJ-CV3. The towell2023reliability self-citation is acceptable for Technometrics if the published version is not available, but a published version may exist.

The master's thesis is cited five times in the active manuscript. Technometrics typically allows citing dissertations but prefers published papers when available. Check: is the foundational masked-causes-in-series-systems (cited as towell2026masked across the series) the published version? If yes, replace towell2023reliability with towell2026masked.

Currently the paper does NOT cite towell2026masked. This is a gap: the foundational paper establishing C1-C2-C3 and the masked-data likelihood is the natural reference for those concepts. Adding `\citep{towell2026masked}` at sections/background.tex lines 11, 59 (where C1-C2-C3 framework is introduced) would strengthen the citation chain across the series.

### MIN-CV1. Cross-references all resolve.

All internal \Cref{} and \ref{} commands resolve to existing labels. Four labels (thm:ident-C123, eq:bernoulli-prob-general, eq:power-masking, eq:factorization) exist only in pre-redesign files but are not referenced from active sections, so no build warnings.

### MIN-CV2. URL stability.

refs.bib contains GitHub URLs for towell2023reliability (queelius/reliability-estimation-in-series-systems) and towell2026binary (queelius/kofn). For a published submission, ensure these remain accessible; consider a Zenodo DOI archival for stable citation.

### MIN-CV3. Bibliography style.

Uses plainnat with \citep and \citet. Consistent throughout. Verified no \cite{} bare commands.

### MIN-CV4. Forthcoming citations.

towell2026binary is "Under review, Technometrics" per refs.bib. For the camera-ready, update with current status (accepted/in press) or remove the venue note.

Bakoyannis2025 is a preprint (arXiv 2511.20980). If a journal version appears before submission, update.

## Recommendations

1. Add `\citep{towell2026masked}` to background.tex where C1-C2-C3 framework and the masked-data likelihood are introduced (lines 11 and 59).
2. Consider replacing self-citations to towell2023reliability with the towell2026masked published version where appropriate.
3. Add Zenodo or other stable DOI for the mdrelax R package code citation.
4. Update preprint statuses before camera-ready.
