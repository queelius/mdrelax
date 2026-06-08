---
schema_version: 1
title: "Sensitivity of Series System Reliability Estimation to the Coarsening Conditions on Masked Failure Data"
short_title: "mdrelax sensitivity paper"
citation_key: "towell2026mdrelax"
stage: revising
format: latex
paper_path: paper/
repo_kind: dual_r_package_and_paper
repo_root: /home/spinoza/github/papers/mdrelax/

paths:
  source: paper/main.tex
  bibliography: paper/refs.bib
  pdf: paper/main.pdf
  sections_dir: paper/sections/
  data_dir: paper/data/
  figures_dir: inst/simulations/figures/
  scripts_dir: paper/
  reviews_dir: .papermill/reviews/
  drafts_dir: .papermill/drafts/

authors:
  - name: "Alexander Towell"
    email: "lex@metafunctor.com"
    orcid: "0000-0001-6443-9897"
    affiliation: "Southern Illinois University Edwardsville, Department of Computer Science"
    corresponding: true

paper_family:
  role: "methodological companion to the foundational towell2026masked"
  cited_as: "towell2026mdrelax"
  cited_by:
    - "five sibling coarsening papers in the masked-series-systems research line"

thesis:
  claim: |
    System-level reliability estimands of a masked series system are robust to
    violations of the three coarsening conditions C1, C2, C3: exactly so for
    exponential components, where the misspecified C1-C2-C3 MLE recovers the
    true system hazard regardless of how candidate sets are formed (via a
    profile-likelihood factorization specific to constant-hazard models). For
    Weibull components, first-order preservation holds under C2 and C3
    violations but fails under C1 violation, where the candidate-set bias
    propagates into the time-portion estimates; the resulting Weibull-C1 bias
    remains bounded below 4 percent across the practical severity range.
    Component-level estimands (especially scale parameters of low-rate
    components) are fragile under all three violations regardless of component
    family.
  novelty: |
    Five contributions in the unified-framework manuscript: (1) calibrated
    severity scales for C1, C2, C3 violations treated as a multi-dimensional
    sensitivity geometry; (2) closed-form local sensitivity indices
    specializing ISNI to set-valued masking, with explicit expressions for
    exponential and tractable forms for Weibull; (3) first-order robustness
    intervals (Definition in Section sec:hierarchy) extending Bakoyannis et
    al. (2025) from Cox semiparametric to parametric series systems with
    set-valued reports; (4) robustness hierarchy classifying functionals by
    first-order sensitivity; (5) reproducible software (mdrelax R package).
    The exponential exact-preservation result via profile factorization is a
    structural finding that strengthens existing first-order results in the
    masked-data literature. A sixth contribution, a formal specification test
    for {C1, C2, C3} from observed candidate-set frequencies, is positioned
    as Future Work in the Discussion (it is sketched and validated in
    paper/specification_test_validation.R but is not folded into the
    manuscript at this revision).
  refined: "2026-05-27"

prior_art:
  last_survey: "2026-05-04"
  status: "current; reframed 2026-05-04 around CAR / ISNI biostatistics lineage"
  key_references:
    - key: HeitjanRubin1991
      role: "Coarsening at random foundation. C1+C2+C3 = CAR + parameter distinctness for candidate-set coarsening map."
    - key: JacobsenKeiding1995
      role: "CAR in general sample spaces. Theoretical foundation of the coarsening framework."
    - key: TroxelMaHeitjan2004
      role: "Index of Local Sensitivity to Nonignorability (ISNI). Direct methodological ancestor of our first-order expansions."
    - key: ZhangHeitjan2006
      role: "ISNI for general coarse-data model. Closest prior local-sensitivity construction; we specialize from interval-censoring to set-valued candidate-set coarsening."
    - key: CopasEguchi2005
      role: "Local model uncertainty and first-order bias under nuisance perturbation. Stylistic ancestor."
    - key: CinelliHazlett2020
      role: "Modern calibrated-bias-bound philosophy (omitted-variable robustness value). Anchor for severity-scale framing."
    - key: MorenoBetancur2015
      role: "Pattern-mixture sensitivity for missing cause in competing risks. Closest direct analog in biostatistics. They self-describe GOF as ad hoc; we address with formal specification test."
    - key: BakoyannisYiannoutsos2015
      role: "Closed-form bias under cause misclassification in cumulative incidence. Template for our series-system bias expressions."
    - key: Bakoyannis2025
      role: "Robustness intervals for Cox cause-specific hazards under MNAR cause. Concurrent work; we extend from binary missing to set-valued masking and from semiparametric Cox to parametric series systems."
    - key: Ebrahimi1996
      role: "Earliest cause-of-failure misclassification likelihood. Point-cause, not set-valued (non-nested per our appendix lemma)."
    - key: VanRompaye2012
      role: "Cox-model misclassification sensitivity. Biostatistics precursor, point-cause."
    - key: HaTsodikov2015
      role: "Semiparametric estimation under misclassified cause. Point-cause, biostatistics."
    - key: Usher-1988
      role: "Original masked data framework"
    - key: Lin-1993
      role: "Exact MLE for masked exponential/Weibull systems"
    - key: towell2023reliability
      role: "Master's thesis establishing the C1-C2-C3 framework this paper builds on"
    - key: LinGuess1994
      role: "Proportional dependent masking (2-component). Alternative to sensitivity approach."
    - key: Guttman1995
      role: "Bayesian dependent masking. Alternative approach."
    - key: Mukhopadhyay2006
      role: "Bayesian incomplete time/cause data. General m-component dependent masking; source of the power-weight C3 violation form."
    - key: CraiuReiser2006
      role: "Dependent competing risks identifiability. Related identifiability results."
  gaps: |
    Major reframing 2026-05-04: prior art search revealed the closest competing
    work is the missing-cause biostatistics literature (Heitjan-Rubin CAR,
    Troxel-Heitjan ISNI, Moreno-Betancur 2015, Bakoyannis 2025), not the
    dependent-masking reliability literature alone. Paper repositioned as a
    parametric-reliability instantiation of mature local-sensitivity machinery,
    with novel contributions in set-valued masking, multi-dimensional severity,
    robustness hierarchy, specification test, and software. New bib entries
    were added to reflect the reframed positioning.

experiments:
  - name: "Sensitivity sweep: C2 violation severity"
    status: complete
    description: |
      Eleven severity levels (s = 0.0 to 1.0), 5-component Weibull system
      (shapes 2.0/1.5/1.2/1.8/1.0, scales 3.0/4.0/5.0/3.5/4.5), n=500,
      B=200 replications per level, tau=5. Standard C1-C2-C3 MLE fitted
      throughout. Metrics: bias, RMSE, coverage for all 10 parameters plus
      system hazard.
    script: paper/run_sensitivity_sweep.R
    data_files:
      - paper/data/sensitivity_sweep_wei.rds
      - paper/data/sensitivity_sweep.rds
    figure_files:
      - inst/simulations/figures/fig_bias_vs_severity.pdf
      - inst/simulations/figures/fig_rmse_vs_severity.pdf
      - inst/simulations/figures/fig_coverage_vs_severity.pdf
  - name: "Non-identifiability demonstration"
    status: complete
    description: |
      Joint (theta, P) estimation with n=2000, 2-component exponential system.
      Shows individual rates confounded with masking probabilities while total
      hazard remains identifiable.
  - name: "Sensitivity sweep: C1 violation severity"
    status: complete
    last_run: "2026-05-27"
    description: |
      Eleven severity levels (alpha_C1 = 0 to 0.5) on 5-component Weibull
      system. Validates Theorem (Pseudo-True Parameter Under C1
      Misspecification) and Proposition (Total Hazard Robustness Under C1
      Violation). Metrics: bias, RMSE, coverage; total hazard preservation
      check; first-order expansion validation. n=500/B=200 and n=2000/B=100
      ablation both present in paper/data/.
    script: paper/run_sensitivity_sweep_c1.R
    data_files:
      - paper/data/sensitivity_sweep_c1.rds
      - paper/data/sensitivity_sweep_c1_n2000.rds
    figure_files:
      - inst/simulations/figures/fig_c1_bias_vs_alpha.pdf
      - inst/simulations/figures/fig_c1_rmse_vs_alpha.pdf
      - inst/simulations/figures/fig_c1_coverage_vs_alpha.pdf
  - name: "Sensitivity sweep: C3 violation severity"
    status: complete
    last_run: "2026-05-27"
    description: |
      5-component Weibull system (same as C2 sweep), 9 alpha levels
      (alpha = 0, 0.25, ..., 2.0), n=500, B=200. Tests whether first-order
      total-hazard preservation extends from exponential (proved exactly) to
      Weibull. Slope-of-bias-vs-alpha hypothesis test reports first-order
      preservation status: linear coefficient +0.004 (p=0.689) at n=500 and
      -0.008 (p=0.204) at n=2000, both consistent with first-order
      preservation. n=500/B=200 and n=2000/B=100 ablation both present in
      paper/data/.
    script: paper/run_sensitivity_sweep_c3.R
    data_files:
      - paper/data/sensitivity_sweep_c3.rds
      - paper/data/sensitivity_sweep_c3_n2000.rds
    figure_files:
      - inst/simulations/figures/fig_c3_bias_vs_alpha.pdf
      - inst/simulations/figures/fig_c3_rmse_vs_alpha.pdf
      - inst/simulations/figures/fig_c3_coverage_vs_alpha.pdf
  - name: "Specification test calibration"
    status: planned
    description: |
      Type I error and power study for the joint specification test
      {C1, C2, C3} against alternatives parameterized by
      (alpha_C1, alpha_C2, alpha_C3). Validates the test as a practical
      diagnostic.

venue:
  target: "Technometrics"
  decision_date: "2026-03-14"
  candidates:
    - name: "Technometrics"
      impact_factor: 3.42
      society: "ASA / ASQ"
      fit: "top pick: theory + application balance, reliability methodology scope, high prestige"
    - name: "Journal of Quality Technology"
      impact_factor: 2.6
      society: "ASQ"
      fit: "backup if Technometrics declines"
    - name: "IEEE Transactions on Reliability"
      fit: "alternative, strong reliability-engineering audience"
    - name: "Reliability Engineering and System Safety"
      fit: "alternative, applied reliability audience"
    - name: "Lifetime Data Analysis"
      impact_factor: 0.91
      fit: "niche safety net; lifetime-data biostatistics framing matches our CAR/ISNI lineage"
    - name: "IISE Transactions Q&R"
      impact_factor: 2.3
      fit: "alternative reliability venue"
  strategy: |
    Technometrics first, JQT if declined, LDA or IEEE TR as safety net. Top
    statistics journals (Biometrics, JASA, AoAS) judged high-risk of rejection
    as rediscovery of biostatistics machinery in a different domain; the
    reliability venues are the safer primary target.

review_history:
  - date: "2026-03-14"
    type: "editorial (multi-agent)"
    outcome: "major-revision"
    artifact: ".papermill/reviews/2026-03-14/review.md"
    summary: |
      C2-only draft. Critical findings: lambda2 spike at s=0.5 undiscussed;
      Theorem 3.2 said "approximately" without a bound; Theorem 3.3 lacked a
      formal proof; the breakdown boundary was ad hoc; abstract overstated
      bias as "order-of-magnitude." Fixes applied: Theorem 3.2 split,
      constructive proof for Theorem 3.3 (m=2), simulation discussion of the
      spike, abstract corrected.
  - date: "2026-05-04"
    type: "prior-art-reassessment"
    outcome: "major-reframing"
    summary: |
      Three-scout literature search surfaced Heitjan-Rubin CAR,
      Troxel-Heitjan ISNI, Moreno-Betancur 2015, and Bakoyannis 2025 as
      load-bearing prior art absent from the original bibliography. C2-only
      framing judged too narrow on novelty. Reframed as unified C1/C2/C3
      sensitivity framework with six contributions anchored to specific
      prior-work gaps. Introduction rewritten, appendix lemma added
      (set-valued vs point-cause non-nesting), framework section extended
      with C1 theorem (with proof) and C3 sketch.
  - date: "2026-05-27"
    type: "migration"
    outcome: "schema-updated"
    summary: |
      Migrated legacy .papermill.md to canonical .papermill/state.md. Refreshed
      stage, paths, and venue list (added IEEE Transactions on Reliability and
      Reliability Engineering and System Safety per family conventions).
      Cleared zero-byte stub paper/../main.tex at repo root to remove the
      LaTeX compile confusion. No content claims changed.
  - date: "2026-05-27"
    type: "critical-fixes"
    outcome: "tables-and-defs-aligned"
    summary: |
      Fixed two Criticals from the 2026-05-27 editorial review.
      (1) Table 1 (component-level fragility) was inconsistent with
      paper/data/*.rds in eight of nine rows; a new
      paper/regenerate_tables.R script rebuilds Tables 1, 2, and 3 from the
      sweep RDS files via a single code path so all three tables now match
      the data. Per-cell changes (n=500, B=200): C1 row k_5 41/0.69 -> 43/0.46
      at alpha=0.50; lambda_1 29/0.91 -> 37/0.98 at alpha=0.50; k_4 21/0.83
      replaced by k_1 30/0.21 at alpha=0.50. C2 row lambda_2 32/near nominal
      replaced by 93/0.98 at s=0.5 (driven by the outlier-spike row already
      discussed in the Coverage Anomalies subsection); k_3 24/0.34 -> 37/0.69
      at s=1.0; lambda_3 18/0.66 replaced by k_2 23/0.44 at s=1.0. C3 row
      lambda_5 149/0.72 at alpha=2.0 -> 163/0.96 at alpha=1.25; lambda_3
      123/1.00 at alpha=1.5 -> 146/0.99 at alpha=1.75; lambda_1 30/0.00
      replaced by lambda_2 41/1.00 at alpha=2.00. The simulations.tex text
      claim that C2 max sys-hazard rel bias is +1.6% updated to +2.6% to
      match the data. The n=2000 ablation claim for lambda_5 under C3 was
      reordered (163% at n=500 vs 149% at n=2000) so the labeling matches
      the data.
      (2) Robustness intervals formally defined in Section sec:hierarchy as
      Definition def:robustness-interval and Proposition
      prop:robustness-coverage, with explicit regularity assumptions and a
      Taylor-remainder bound. Theorem 3.9 (C3 pseudo-true expansion) now
      states the total-preservation identity and the score-Jacobian
      linear system for general m, removing the forward reference to
      nonexistent supplementary material; the closed-form per-component
      Delta_j^C3 for the specific power-weight perturbation as
      implemented in the mdrelax package is left for the follow-up
      identifiability companion (towell2026binary), since the package
      normalization (dividing power weights by their max over all j, not
      over j != k as the paper formula reads) is not exactly the form
      analytically expanded here.
      Additionally: ORCID added to title block; canonical name updated to
      Alexander Towell; abstract rewritten for the unified C1/C2/C3
      framework; graphicspath widened to include paper/figures/.
  - date: "2026-06-04"
    type: "editorial (multi-agent)"
    outcome: "minor-revision"
    artifact: ".papermill/reviews/2026-06-04/review.md"
    summary: |
      Final pre-submission pass (per-specialist files complete; synthesis
      review.md written 2026-06-05 after the launching orchestrator was
      rate-limited). Both 2026-05-27 Criticals confirmed FIXED by direct
      verification (Table 1 regenerated from data; robustness intervals
      now formally Definition + Proposition in Section sec:hierarchy). No
      new Critical. One carried Major (C3 n=2000 ablation under-powered as
      a positive preservation claim, mitigated by the structural argument
      and the powered C1 rejection). Minors: title vs Zenodo metadata
      mismatch; robustness interval never demonstrated on the application
      data; one-sentence foundational-paper boundary missing; soften
      "multi-dimensional sensitivity geometry"; dead Stefanski2002 entry;
      orphan section files / figure path / 1276 vanity count for the
      tarball. Single most-important item: reconcile the
      title-vs-Zenodo-metadata mismatch before submission.
  - date: "2026-06-08"
    type: "editorial (multi-agent)"
    outcome: "minor-revision"
    artifact: ".papermill/reviews/2026-06-08/review.md"
    summary: |
      Independent verification pass over the 2026-06-04 state plus four new
      findings. Re-verified: clean build (0 undefined, 0 multiply-defined,
      42 pp); Tables 1-3 reproduced cell-by-cell from committed data via
      regenerate_tables.R; application system hazard 0.00307 reproduced
      from the canonical MLE; abstract numeric claims (sys-hazard bias
      under 4% throughout, 100%+ component bias under C3, zero coverage for
      lambda_1 under C3) all confirmed against the sweep RDS files; proof
      core (Thm 3.6 exact preservation, Lemma C.2 non-nesting, Prop 3.21
      robustness coverage) re-confirmed sound. NEW: (Major) the
      robustness-interval functions (ri_first_order, ri_simulation in
      R/robustness_intervals.R) are unexported (absent from NAMESPACE, no
      man pages) and never demonstrated on data, so the software /
      robustness-interval contribution claim is partially untrue as
      shipped; (Minor) application Table 4 s=0 baseline scales (897, 847)
      drift from the package canonical Guo MLE (909, 840), though h_T is
      unaffected; (Minor) residual C2-only framing in background.tex L84
      and sensitivity_framework.tex L6-7 contradicts the unified scope;
      (Minor) live Zenodo record 10.5281/zenodo.20468529 still shows the
      old "non-informative masking assumption" title while the PDF /
      .zenodo.json / CITATION.cff use the new "Coarsening Conditions"
      title. The 2026-06-04 minor items remain unaddressed in the
      manuscript (no content changed since). Verdict minor-revision
      (verging on ready). Single most-important item: export + demonstrate
      the robustness interval (closes the software-claim gap and shows the
      headline tool), then reconcile the title across the published Zenodo
      record.

r_package:
  name: mdrelax
  version: "1.0.0"
  description: "Relaxed Candidate Set Models for Masked Data in Series Systems"
  dependencies: ["stats"]
  suggests: ["knitr", "rmarkdown", "devtools", "numDeriv", "testthat (>= 3.0.0)"]
  license: "GPL (>= 3)"
  url: "https://queelius.github.io/mdrelax/"
  bug_reports: "https://github.com/queelius/mdrelax/issues"
  tests: "~240 test_that blocks across 11 files (all passing as of last run)"
  sensitivity_tooling_status: |
    ISNI computation, robustness-interval construction, and the specification
    test are all planned as new R package exports tied to paper Sections 3-4.

dual_repo_notes: |
  This repository is dual-purpose: an R package at the root (DESCRIPTION,
  NAMESPACE, R/, tests/, man/, vignettes/) AND a paper sub-repo at paper/.
  Papermill state lives at the repo root (.papermill/state.md) so that the
  existing reviews directory (.papermill/reviews/) stays co-located with it.
  All paper-source paths in this state file are relative to the repo root and
  begin with paper/ for clarity.

  A zero-byte main.tex existed at the repo root (left over from a misplaced
  build); the actual master document is paper/main.tex. The R package itself
  is unaffected.

structure:
  sections_used:
    - introduction.tex
    - background.tex
    - sensitivity_framework.tex
    - simulations.tex
    - application.tex
    - discussion.tex
    - conclusion.tex
    - appendix.tex
  sections_present_but_unused:
    - relaxed_models.tex   # pre-redesign remnant
    - identifiability.tex  # pre-redesign remnant
    - introduction.pre-unified.tex  # backup of the C2-only intro
    - specification_test_sketch.tex # working draft, not included in main.tex
  outline:
    - section: "Introduction"
      file: introduction.tex
      key_content: "Unified C1/C2/C3 framing; six contributions anchored to prior-work gaps"
    - section: "Background"
      file: background.tex
      key_content: "Series model, masked data, C1-C2-C3 likelihood, dependent masking literature"
    - section: "Sensitivity Framework"
      file: sensitivity_framework.tex
      key_content: "Likelihood under C1 alone; Bernoulli model; C2 misspecification theorem; identifiability trap; C1 sensitivity (with proof); C3 sensitivity; robustness hierarchy"
    - section: "Simulation Study"
      file: simulations.tex
      key_content: "C2 severity sweep results (C1 and C3 sweeps planned)"
    - section: "Application"
      file: application.tex
      key_content: "Guo et al. turbine engine data sensitivity analysis"
    - section: "Discussion"
      file: discussion.tex
      key_content: "Robustness hierarchy interpretation, practical guidance, limitations"
    - section: "Conclusion"
      file: conclusion.tex
      key_content: "Six contributions summarized"
    - section: "Appendix A1 Score Functions"
      file: appendix.tex
      key_content: "Misspecified vs true score derivations"
    - section: "Appendix A2 Non-Identifiability Evidence"
      file: appendix.tex
      key_content: "Joint estimation failure demonstration"
    - section: "Appendix A3 Set-Valued vs Point-Cause"
      file: appendix.tex
      key_content: "Non-nesting lemma; consequences for inference"
    - section: "Appendix A4 Software"
      file: appendix.tex
      key_content: "mdrelax R package description"

build:
  command: "cd paper && make pdf"
  fallback_command: "cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main"
  table_regenerator: "Rscript paper/regenerate_tables.R"
  last_build: "2026-05-27"
  status: "verified 2026-05-27 after Critical fixes: 41-page main.pdf, two harmless hyperref warnings about math-in-title for sec:sim-ablation, no undefined references"

polish_checklist:
  done:
    - "Title and abstract finalized (rewritten 2026-05-27 for the unified C1/C2/C3 framework)"
    - "Bibliography compiles cleanly via natbib (refs.bib, plainnat style)"
    - "Section files used by main.tex are all present"
    - "Author block carries canonical name, affiliation, email, and ORCID 0000-0001-6443-9897 (added 2026-05-27)"
    - "Migrated legacy .papermill.md to canonical .papermill/state.md (2026-05-27)"
    - "Cleared zero-byte main.tex stub at repo root (2026-05-27)"
    - "Table 1 (component-level fragility) regenerated from paper/data/*.rds via paper/regenerate_tables.R; numbers now match the sweep data exactly (2026-05-27)"
    - "Robustness intervals formally defined in Section sec:hierarchy (Definition def:robustness-interval, Proposition prop:robustness-coverage); previously deferred to nonexistent supplementary material (2026-05-27)"
    - "Theorem 3.9 (C3 pseudo-true expansion) closed-form for m=2 inlined and the general-m linear system written out; supplementary-material forward references removed (2026-05-27)"
    - "Graphics path widened to also search paper/figures/; upward dependency on inst/simulations/figures/ documented inline in main.tex (2026-05-27)"
    - "Contribution count aligned to five in introduction.tex and conclusion.tex; specification test stays in Future Work in discussion.tex (2026-05-27)"
  pending_user:
    - "Add real-data application section as Technometrics expects (planned: Guo turbine data)"
    - "Reformat to Technometrics ASA template before submission"
    - "Remove pre-redesign remnant section files (relaxed_models.tex, identifiability.tex, introduction.pre-unified.tex) once unified framework is finalized"
    - "Update background.tex to explicitly introduce CAR and parameter-distinctness vocabulary"
    - "Decide whether to copy the needed PDFs into paper/figures/ for a fully self-contained camera-ready bundle"

remaining_work:
  - "Implement C1 first-order bias for m >= 3 (Appendix C1-general)"
  - "Implement C3 first-order bias for m = 2 closed form and general m (Appendix C3-general)"
  - "Develop the formal specification test for {C1, C2, C3} from candidate-set frequencies"
  - "Run C1 and C3 sensitivity sweeps (analogous to existing C2 sweep)"
  - "Verify first-order expansions empirically against simulation"
  - "Build robustness-interval construction in mdrelax R package"
  - "Update simulations.tex, application.tex, discussion.tex to reflect unified framework"

log:
  - date: "2026-02-17"
    note: "Initialized papermill from existing draft + R package. Advanced draft with complete proofs, two simulation studies, and the test suite."
  - date: "2026-02-17"
    note: "Venue assessment. Traditional journal, not JOSS. Paper is primarily theoretical."
  - date: "2026-02-28"
    note: "Major paper redesign. Rewritten as focused sensitivity analysis of C2 violations in 5-component Weibull system. New sections: sensitivity_framework.tex; updated simulations.tex. Dropped relaxed_models.tex and identifiability.tex from main.tex."
  - date: "2026-03-14"
    note: "Refreshed .papermill.md to match redesigned paper structure."
  - date: "2026-03-14"
    note: "Prior-art survey completed. Added new references (White 1982, Huber 1967, Tsiatis 1975, Rosenbaum 2002, etc.). Updated background.tex, sensitivity_framework.tex, discussion.tex with citations."
  - date: "2026-03-14"
    note: "Venue selected: Technometrics. Backup: JQT, LDA."
  - date: "2026-03-14"
    note: "Editorial review completed (major-revision). Fixes applied."
  - date: "2026-05-04"
    note: "Prior-art reassessment via three parallel literature scouts surfaced Heitjan-Rubin CAR, Jacobsen-Keiding 1995, Troxel-Heitjan ISNI, Zhang-Heitjan 2006, Copas-Eguchi 2005, Cinelli-Hazlett 2020, Moreno-Betancur 2015, Bakoyannis-Yiannoutsos 2015, Bakoyannis 2025, Ebrahimi 1996, Van Rompaye 2012, Ha-Tsodikov 2015 as load-bearing prior art. C2-only framing judged too narrow."
  - date: "2026-05-04"
    note: "Paper reframed as unified C1/C2/C3 sensitivity framework with six contributions. Introduction rewritten (old version preserved as introduction.pre-unified.tex). New bib entries added. Appendix lemma proves set-valued masking and point-cause misclassification are non-nested families."
  - date: "2026-05-04"
    note: "Framework section extended with C1 sensitivity subsection: Theorem (Pseudo-True Parameter Under C1 Misspecification, m=2) with full proof; Proposition (Total Hazard Robustness Under C1 Violation); Remark on interpretation. C3 sensitivity subsection added with theorem statement and proof sketch deferred. New subsection on Local Sensitivity Indices and Robustness Hierarchy defines ISNI specialized to candidate-set coarsening and states a corollary classifying estimands."
  - date: "2026-05-27"
    note: "Papermill migration. Created canonical .papermill/state.md at repo root; legacy .papermill.md kept in place as a fallback until the next stable build. paper_path set to paper/. Reviews remain at .papermill/reviews/. Cleared zero-byte main.tex stub at repo root."
---

# mdrelax sensitivity paper

This file is the canonical papermill state for the paper at `paper/main.tex`.
The legacy `.papermill.md` (older pre-schema format) at the repo root remains
as a fallback but should be considered superseded by this file.

The repository is dual-purpose: it is both an R package (DESCRIPTION,
NAMESPACE, R/, tests/, man/, vignettes/ at the repo root) and a paper
sub-repo (`paper/` directory). Papermill state lives at the repo root because
the existing `.papermill/reviews/` directory is co-located here and because
all five sibling coarsening papers in this research line cite this paper as
`towell2026mdrelax`.

See the YAML front-matter above for the structured state. Build status and
polish notes are tracked under `build:` and `polish_checklist:` respectively.
