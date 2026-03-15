# Paper Completion Plan: Sensitivity Analysis of C2 Violations

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete all remaining work to bring the redesigned sensitivity-analysis paper to submission readiness.

**Architecture:** Five independent workstreams: (1) refresh project state file to match the redesigned paper, (2) rebuild the PDF to verify the manuscript compiles, (3) conduct a prior-art survey for masked data / sensitivity analysis literature, (4) select a target venue, (5) run a multi-agent editorial review. Streams 1-2 are prerequisites for stream 5; streams 3-4 can run in parallel.

**Tech Stack:** LaTeX (pdflatex + bibtex), R (mdrelax package), papermill skills

---

## Chunk 1: Refresh State & Rebuild PDF

### Task 1: Update `.papermill.md` to match redesigned paper

**Files:**
- Modify: `.papermill.md`

The paper was completely rewritten in commit `109fd8a`. The current state file reflects the OLD framing (general relaxed-model framework). The NEW framing is: "Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption."

- [ ] **Step 1: Rewrite `.papermill.md` with correct metadata**

Update the following fields to match the redesigned paper:
- `title`: "Sensitivity of Series System Reliability Estimation to the Non-Informative Masking Assumption"
- `thesis.claim`: When C2 is violated but assumed to hold, the MLE converges to a pseudo-true parameter. System hazard remains consistent; individual components absorb masking asymmetry. Joint estimation fails due to non-identifiability. Sensitivity analysis is recommended over model elaboration.
- `thesis.novelty`: First systematic characterization of C2 misspecification bias in masked series systems, proving system-hazard robustness and component-level degradation. Bernoulli perturbation model for controlled C2 violations. Non-identifiability proof for joint (θ, P) estimation.
- `experiments`: Update to reflect the sensitivity sweep design (11 severity levels, 5-component Weibull, B=200, n=500)
- `structure`: Update section table to match current `main.tex` includes (no more `relaxed_models.tex` or `identifiability.tex` as standalone — they're merged into `sensitivity_framework.tex`)

- [ ] **Step 2: Verify updated state file is consistent with paper**

Cross-check: main.tex includes match section table, figure references match existing files, experiment descriptions match simulation scripts.

### Task 2: Rebuild PDF

**Files:**
- Build: `paper/main.pdf`

- [ ] **Step 1: Compile LaTeX**

```bash
cd paper && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

- [ ] **Step 2: Check for warnings**

Review build output for: undefined references, missing citations, missing figures, overfull hboxes.

- [ ] **Step 3: Fix any issues found**

Address any compilation errors or missing references.

## Chunk 2: Prior-Art Survey

### Task 3: Conduct systematic literature search

This is the biggest gap. No prior-art survey has been done (`last_survey: null`).

The paper's positioning: "sensitivity of C1-C2-C3 estimator to C2 violations" — complementary to the dependent masking literature (Lin & Guess 1994, Guttman 1995, Mukhopadhyay 2006, Craiu & Reiser 2006). Need to find:

1. **Sensitivity analysis in reliability** — papers studying robustness of reliability estimators to assumption violations
2. **Masked data / competing risks** — any masked data papers not already cited
3. **Model misspecification in survival analysis** — White (1982) pseudo-true parameter theory, Huber sandwich estimator
4. **Non-identifiability in mixture/masked models** — results on when masking + component parameters are jointly non-identifiable

- [ ] **Step 1: Run papermill prior-art survey**

Use `/papermill:prior-art` to conduct systematic search. Focus areas:
- "masked failure data" + "dependent masking" + sensitivity
- "series system reliability" + "model misspecification"
- "competing risks" + "non-informative censoring" + sensitivity
- "pseudo-true parameter" + misspecification (White 1982 lineage)

- [ ] **Step 2: Cross-reference findings with existing refs.bib**

Identify which found papers are already cited vs. need adding.

- [ ] **Step 3: Update refs.bib with new references**

Add BibTeX entries for important missing references.

- [ ] **Step 4: Update paper sections with new citations**

Primary targets: `background.tex` (prior work subsection), `sensitivity_framework.tex` (misspecification theorem context), `discussion.tex` (comparison to alternative approaches).

## Chunk 3: Venue Selection

### Task 4: Evaluate and select target venue

**Files:**
- Modify: `.papermill.md` (venue section)

Current candidates: JSS, Technometrics, Lifetime Data Analysis, RESS.

- [ ] **Step 1: Run papermill venue evaluation**

Use `/papermill:venue` to evaluate candidates against paper characteristics:
- Theory-heavy (~60%): 2 theorems with proofs, non-identifiability result
- Simulation study: systematic sensitivity sweep
- Software: R package with 1276 tests
- ~20 pages
- Single author

- [ ] **Step 2: Update `.papermill.md` with selected venue**

Set `venue.target` and add notes on formatting requirements.

## Chunk 4: Editorial Review

### Task 5: Run multi-agent review

**Prerequisite:** Tasks 1-2 complete (state file updated, PDF compiles).

- [ ] **Step 1: Run papermill review**

Use `/papermill:review` for comprehensive editorial feedback.

- [ ] **Step 2: Triage review findings**

Categorize as: must-fix (blocking), should-fix (improves quality), nice-to-have.

- [ ] **Step 3: Address must-fix issues**

Implement fixes in the manuscript.

- [ ] **Step 4: Rebuild PDF after fixes**

```bash
cd paper && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```
