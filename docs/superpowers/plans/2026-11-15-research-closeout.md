# Research Closeout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Finish one existing masked-data result with reproducible evidence and discoverable publication records.

**Architecture:** Keep the mathematical note and executable evidence in mdrelax. Link a structured central record and the R registry to the owning source and existing paper archive.

**Tech Stack:** R, testthat, covr, Python standard-library JSON/HTTP, Git, Zenodo public API.

**Spec:** ../specs/2026-11-15-research-closeout.md

## Global Constraints

- No new application, package API, or general paper correctness claim.
- Exact fitted identity is separate from population consistency and unbiasedness.
- Preserve existing archive families, unrelated edits, and predecessor attribution.
- Structured evidence and dated source pins are authoritative.

## Review Focus

- Censoring: exposure includes all retained observations; failures alone contribute log-rate.
- Full masking: identified total must not be mistaken for identified components.
- Selection: dropping early failures changes the inferential target.
- Boundaries: zero failures and binding rate bounds are outside the interior identity.
- Publication: historical PDF must not be claimed to contain newly added checks.

### Task 1: Verify the selected result

Files: tests/testthat/test-total-rate-factorization.R; paper/reproduction/exponential-total-rate/verify.R and reproduce.sh.

- [x] Run existing package suite and inspect coverage.
- [x] Add independent factorization, derivative and numerical fit checks, including censoring and full masking.
- [x] Add executable finite examples and rejecting controls for changed component rankings, rate bounds, zero failures, and selective retention.
- [x] Run against a clean installed source package, with no user R profile or preloaded workspace.

### Task 2: Explain and preserve

Files: paper/reproduction/exponential-total-rate/README.md and evidence.json; README.md; MAINTENANCE.md.

- [x] Explain claim, proof, assumptions, limitations, reproduction, citation and credit.
- [x] Verify existing archive metadata and PDF checksum; record exact DOI and version.
- [x] Correct stale export status and distinguish AI review from external peer review.
- [x] Obtain one independent whole-change review, resolve important findings, rerun relevant checks.
- [ ] Commit and publish scoped source changes; verify public source.

### Task 3: Connect the ecosystem

Files: research-records/records/masked-exponential-total-rate/{record.json,README.md}; catalog.json; registry README.md.

- [ ] Add one stable central record with exact source pins, evidence and archive scope.
- [ ] Make the verified research result the ecosystem's entry point; retain the latency example as an example.
- [ ] Validate JSON, links, hashes and diff, publish, and verify public content.

## Execution ledger

2026-11-15: User authorized bounded completion and publication under standing research instructions. Execute inline; no further design/permission round is needed. Isolated source worktree /tmp/mdrelax-closeout-20261115 preserves the original checkout. Existing archive verified via public API; v0.2.0 contains the exact factorization. New checks will be identified separately from that archive.

2026-11-15: Tasks 1 and source documentation complete. Independent AI review of 2932b16 reran the reproduction, suite and coverage and found no critical, important or reportable minor issues. Ruling: broader Weibull/relaxed-model and historical manuscript correctness, cross-platform verification and legacy simulation-loader repairs remain outside this bounded result; those limitations are public. Raw sessionInfo whitespace is retained as evidence. Publication will reuse the paper family for a supplemental reproducibility edition, preserving its PDF byte-for-byte; no new concept DOI.
