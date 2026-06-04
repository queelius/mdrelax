# Format Validator Report

**Date**: 2026-05-27
**Paper**: Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption

## Build Status

The current main.pdf (39 pages, 541 KB, last build 2026-05-27 20:47) builds cleanly:
- Zero LaTeX warnings (verified via grep on main.log)
- Zero undefined references
- Zero multiply-defined labels
- All figures resolve via the graphicspath `../inst/simulations/figures/`
- Bibliography compiles via natbib + plainnat (main.bbl present)

## Findings

### MAJ-F1. ORCID missing from main.tex title block.

main.tex lines 97-101 contain only:
```
\author{
    Alex Towell\\
    Southern Illinois University Edwardsville\\
    \texttt{lex@metafunctor.com}
}
```

DESCRIPTION file lists ORCID 0000-0001-6443-9897. State file polish_checklist flags this as pending. For Technometrics, ORCID is recommended (not strictly required but expected for new submissions). Add:

```
\author{
    Alex Towell\thanks{Department of Computer Science, Southern Illinois University Edwardsville. ORCID 0000-0001-6443-9897.}\\
    \texttt{lex@metafunctor.com}
}
```

The state file also notes the full name should be Alexander Richard Towell or Alexander Towell (formal); current "Alex Towell" is casual.

### MAJ-F2. Cross-package graphics path is awkward for Technometrics camera-ready.

main.tex line 18: `\graphicspath{{../inst/simulations/figures/}}`

This pulls figures from outside the paper/ directory (from the R package's inst/simulations/figures/ directory). For Technometrics camera-ready submission via tarball, this creates a directory dependency that will likely break ManuscriptCentral / ASA submission systems.

Options:
(a) Symlink or copy figures into paper/figures/ before tarballing.
(b) Use a Makefile target `make submission-tarball` that flattens the figure tree.
(c) Modify the R simulation scripts to write figures to paper/figures/ instead of inst/simulations/figures/.

State file already flags this as pending; recommend (c) to permanently resolve.

### MAJ-F3. Pre-redesign remnant section files in paper/sections/.

Four files exist in paper/sections/ but are not included in main.tex:
- relaxed_models.tex (502 lines)
- identifiability.tex (382 lines)
- introduction.pre-unified.tex (74 lines)
- specification_test_sketch.tex (207 lines)

These do not affect the build, but they contain orphan labels (thm:ident-C123, eq:bernoulli-prob-general, eq:power-masking, eq:factorization). For a Technometrics submission tarball, including these creates reviewer confusion. State file polish_checklist flags removing them as pending.

Recommend: move them to paper/sections/_unused/ or delete them before submission. The introduction.pre-unified.tex backup is safe in git history; no need to keep in the source tree.

### MAJ-F4. Document class is article, not Technometrics ASA template.

main.tex line 5: `\documentclass[11pt,letterpaper]{article}`

State file polish_checklist explicitly flags: "Reformat to Technometrics ASA template before submission." The ASA template (asaproc or specific Technometrics class) is preferred. For now this is acceptable for editorial review; for submission, retemplate.

### MAJ-F5. Page count.

Current PDF is 39 pages (double-spaced 11pt). Technometrics page limits vary by article type (regular article approximately 30 pages double-spaced including all). At 39 pages, the paper is over typical limit. Suggestion: when retemplated to Technometrics, consolidate the Limitations section (currently 6 bulleted items spanning ~60 lines) and consider moving the appendix software section to a supplementary materials document.

### MIN-F1. Vanity count flagged.

appendix.tex line 174: "1276 unit tests verifying likelihood correctness, score--gradient agreement, FIM consistency, and MLE convergence properties."

The "1276" is a vanity count that the soul plugin convention flags. Replace with: "The package's test suite verifies likelihood correctness, score-gradient agreement, FIM consistency, and MLE convergence properties." (Drop the integer.)

### MIN-F2. \todo and \placeholder macros.

Defined in main.tex lines 42-43 but not used in any active section file (verified via grep). Safe to leave defined or remove.

### MIN-F3. html_paper/ directory.

Contains LaTeXML output. Not part of the deliverable for a journal submission. Recommend gitignoring or moving to a build artifacts directory.

### MIN-F4. Hyperref Unicode warnings.

main_build.log shows two "Package hyperref Warning: Token not allowed in a PDF string (Unicode)" warnings. These are caused by the equation-laden title or section titles being interpreted as PDF metadata. Resolvable with \texorpdfstring around offending math in section headings, or by switching to a hyperref-aware title. Cosmetic; not a blocker.

### MIN-F5. R package version consistency.

DESCRIPTION says Version 1.0.0; appendix.tex line 166 says "mdrelax R package (version 1.0.0)". Consistent. (The 2026-03-14 review's 1.0.0 vs 1.1.0 inconsistency is resolved.)

## Cross-Verification

All findings verified by:
- Reading main.tex and counting included sections.
- Running `pdflatex -interaction=nonstopmode main.tex` and inspecting the log.
- Listing paper/sections/ to find remnant files.
- Listing the inst/simulations/figures/ directory to confirm the cross-package path resolves.
- Grepping for \todo, \placeholder, and version strings.
