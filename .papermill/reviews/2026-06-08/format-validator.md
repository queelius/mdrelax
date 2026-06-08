# Format Validator Report

**Date**: 2026-06-08
**Focus**: build verification, label resolution, venue formatting.

## Build (verified this pass)

- `make -C paper pdf` (pdflatex x3 + bibtex) completes with exit 0.
- **0 undefined references** (`LC_ALL=C grep -ai undefined main.log | grep -vi "Font shape" | wc -l` = 0).
- **0 multiply-defined labels** in the compiled document.
- A transient "Label(s) may have changed. Rerun" warning appeared after the first full pass but clears on a subsequent pdflatex run (cross-references settle); the final log has no rerun warning.
- 42 pages, double-spaced (up 1 page from the 41 reported 2026-06-04).
- Two harmless `hyperref Warning: Token not allowed in a PDF string (Unicode)` remain (the math in the title; affects only PDF bookmark strings).
- Per task instructions, the tracked `paper/main.pdf` was reverted after building (`git checkout -- paper/main.pdf`); the tree is clean.

## Verdict on prior format items

- **ORCID (MAJ-9): FIXED** (confirmed, main.tex lines 100-105).
- **Figure path (MAJ-5): STILL OPEN.** main.tex line 26 is `\graphicspath{{../inst/simulations/figures/}{figures/}}`. The first path reaches up out of the paper sub-repo into the R package's inst/ tree, which will not exist in a ManuscriptCentral/ASA submission tarball. paper/figures/ is on the path as a fallback; the build currently succeeds because the figures resolve via one of the two paths. Before submission: copy the needed figure PDFs (fig_c1_bias_vs_alpha.pdf, fig_bias_vs_severity.pdf, fig_c3_bias_vs_alpha.pdf and the rmse/coverage analogues) into paper/figures/ and verify the build with ONLY paper/figures/ reachable, then drop the ../inst path.
- **Orphan pre-redesign section files (MAJ-7): STILL OPEN and demonstrably harmful to the tarball.** Four files are NOT `\input` by main.tex (confirmed: main.tex inputs only introduction, background, sensitivity_framework, simulations, application, discussion, conclusion, appendix):
  - `sections/relaxed_models.tex` -- redefines eq:P-matrix, eq:bernoulli-like, eq:like-C1, and many others.
  - `sections/identifiability.tex` -- redefines thm:misspec-C2, thm:misspec-C3, sec:identifiability, sec:misspec.
  - `sections/introduction.pre-unified.tex` -- redefines sec:introduction.
  - `sections/specification_test_sketch.tex` -- a working draft (defines its own spec-test labels).
  They do not affect the build (not included), but a referee/editor inspecting the source tarball will see duplicate `\label` definitions for core theorems and sections. Move to `sections/_unused/` or delete before submission.
- **Vanity count "1276 unit tests" (MIN-4): STILL OPEN** (appendix.tex line 174).
- **Hyperref Unicode warnings (MIN-7): STILL OPEN** (two, cosmetic; wrap the math in the title with `\texorpdfstring`).

## NEW format/consistency issue: title vs published metadata

The on-page title (main.tex) and .zenodo.json and CITATION.cff all now use "Sensitivity of Series System Reliability Estimation to the Coarsening Conditions on Masked Failure Data". Verified live this pass: the **published Zenodo record (latest version 10.5281/zenodo.20468529 under concept 10.5281/zenodo.20414727) still shows the OLD title** "Sensitivity of series system reliability estimation to the non-informative masking assumption". So the deposited .zenodo.json (which would govern the NEXT version) and the live published landing page disagree, and the PDF title differs from the current DOI landing page. The refs.bib header comment and state.md `title:` also still carry the old title. Reconcile: confirm the on-page "Coarsening Conditions" title as canonical, push a new Zenodo version (or edit metadata) so the published landing page matches, and update refs.bib header + state.md.

## Venue formatting (Technometrics)

- **Document class is `article`, not the ASA/Technometrics template (carryover):** retemplate to the Technometrics ASA style before submission (currently 11pt letterpaper, 1in margins, double-spaced).
- **Page length (carryover):** 42 pages double-spaced, over typical Technometrics regular-article length. Levers: move the Software appendix (app:software) and possibly the Score Functions appendix (app:scores) to supplementary materials; tighten the Limitations bullet list. The robustness-interval material in 3.8 is core and should stay; the Application is short and is not the place to cut.

## Reproducibility hooks (production-positive)

- `regenerate_tables.R` regenerates Tables 1-3 from committed data and is named in the Table 1 caption. Verified to reproduce all three tables exactly this pass. Strong.
- `run_sensitivity_sweep*.R` reproduce the figures.
- GAP (cross-ref methodology lens): the Application table (Table 4, tab:guo-sensitivity) is NOT regenerated from a committed path and its s=0 baseline scales drift from the package's canonical Guo MLE. Bringing it under a committed regenerator (like the simulation tables) would close the last reproducibility soft spot.

## Bottom line (format)

Build is clean (0 undefined, 0 multiply-defined, rerun warning clears, 42 pp). Submission-priority order: (1) co-locate figures into paper/figures/ and drop the ../inst path; (2) remove the four orphan section files (colliding labels in the tarball); (3) reconcile the title vs published Zenodo metadata; (4) retemplate to ASA + trim to the page limit (appendices to supplement); (5) drop the 1276 vanity count and the two hyperref Unicode warnings. None blocks correctness; (1) and (2) block a clean tarball.
