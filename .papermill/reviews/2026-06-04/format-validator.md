# Format Validator Report

**Date**: 2026-06-04
**Focus**: build verification, label resolution, venue formatting.

## Build

- `make pdf` (pdflatex x3 + bibtex) completes with exit 0 on a forced clean rebuild (touched main.tex).
- No undefined references, no "multiply defined" labels, no LaTeX rerun warnings in the included build. All \Cref / \cite resolve.
- 41 pages, double-spaced.

## Verdict on prior format issues

- **ORCID missing (MAJ-9): FIXED.** main.tex lines 105-110 author block now includes "ORCID: 0000-0001-6443-9897" with an \href to orcid.org. (Implemented as a plain line in the author block rather than \thanks, which is fine for the article class.)
- **Figure path (MAJ-5): NOT FIXED.** main.tex line 26 still has \graphicspath{{../inst/simulations/figures/}{figures/}}. The first path reaches up out of the paper sub-repo into the R package's inst/ directory. The paper-local figures/ directory IS now on the path as a fallback (good intent), but I confirmed figures/ exists; check that the actual figure PDFs (fig_c1_bias_vs_alpha.pdf, fig_bias_vs_severity.pdf, fig_c3_bias_vs_alpha.pdf, and the rmse/coverage analogues) are present in paper/figures/ and not only in ../inst/. For a ManuscriptCentral/ASA tarball the figures must be co-located with the source; the ../inst path will not exist in the submission bundle. Recommend: copy the needed PDFs into paper/figures/ and drop the ../inst path (or keep both, but verify the build succeeds with ONLY paper/figures/ reachable).
- **Orphan pre-redesign section files (MAJ-7): NOT FIXED, and now demonstrably harmful.** Four files are not \input by main.tex: identifiability.tex, relaxed_models.tex, introduction.pre-unified.tex, specification_test_sketch.tex. They do not affect the build, BUT they re-define labels that collide with the live files: thm:misspec-C2, thm:misspec-C3, eq:P-matrix, eq:like-C1, eq:bernoulli-like, sec:introduction, sec:identifiability, sec:misspec. A referee or editor inspecting the source tarball will see duplicate \label definitions for core theorems. Move these four files to sections/_unused/ or delete them before submission.
- **Vanity count "1276 unit tests" (MIN-4): NOT FIXED.** appendix.tex line 174. Replace with a description of what the suite verifies (likelihood correctness, score-gradient agreement, FIM consistency, MLE convergence) and drop the integer.
- **Hyperref Unicode warnings (MIN-7): NOT FIXED.** Two "Token not allowed in a PDF string (Unicode)" warnings remain in main.log. Cosmetic (affects only PDF bookmark strings). Wrap the offending math/symbol in section or title strings with \texorpdfstring.
- **html_paper/ artifact (MIN-11):** not present in the current paper/ listing, so either already cleaned or never regenerated. The Makefile still has an html target writing to html_paper/; ensure that directory is gitignored so a stray `make html` does not pollute the tarball.

## NEW format/consistency issue

- **Title vs metadata mismatch (MINOR, also flagged by prose-auditor):** the on-page title (main.tex line 102) is "...to the Coarsening Conditions on Masked Failure Data" while refs.bib header, state.md, and the Zenodo concept-DOI record (verified live: citation_title = "Sensitivity of series system reliability estimation to the non-informative masking assumption") use the older "...Non-Informative Masking Assumption". The PDF title and the DOI landing-page title will not match at submission. Pick one and reconcile the Zenodo deposit metadata.

## Venue formatting (Technometrics)

- **Document class is `article`, not the ASA/Technometrics template (MIN-5):** retemplating to the Technometrics ASA style is needed before submission. Currently 11pt letterpaper, 1in margins, double-spaced.
- **Page length (MIN-6):** 41 pages double-spaced (up from 39 at the last review, because Section 3.8 robustness-interval material and the rewritten abstract were added). This is over the typical Technometrics regular-article length. Levers: move the Software appendix (Appendix D) and possibly the Score Functions appendix (Appendix A) to supplementary materials; tighten the Limitations bullet list; the Application section is short so it is not the place to cut. The robustness-interval material is core and should stay.
- backups (IEEE TR, Lifetime Data Analysis, JSS/R Journal): page limits are less binding for LDA and the R-package venues; for IEEE TR a two-column retemplate and significant compression would be required.

## Reproducibility hooks (production-positive)

- regenerate_tables.R regenerates Tables 1-3 from committed data and is referenced in the Table 1 caption. Strong.
- run_sensitivity_sweep*.R reproduce the figures. The figure-path issue above is the only thing standing between the current tree and a self-contained tarball.

## Bottom line (format)

Build is clean, ORCID fixed, no broken references. Remaining items, in submission-priority order: (1) co-locate figures into paper/figures/ and drop the ../inst path; (2) remove the four orphan section files (they create colliding labels in the tarball); (3) reconcile the title vs Zenodo metadata; (4) retemplate to ASA + trim to the Technometrics page limit (move appendices to supplement); (5) drop the 1276 vanity count and the two hyperref Unicode warnings. None block correctness; (1) and (2) block a clean tarball.
