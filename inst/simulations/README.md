# Simulation Studies for Relaxed Candidate Set Models

This directory contains the complete simulation framework for the paper
"Relaxed Candidate Set Models for Masked Data in Series Systems."

## Overview

The simulation studies examine three key research questions:

1. **KL-Divergence Efficiency**: How does informativeness of masking (measured
   by KL-divergence from non-informative baseline) affect MLE efficiency?

2. **Misspecification Bias**: What bias arises when condition C2 (non-informative
   masking) is violated but assumed to hold in the analysis?

3. **Identifiability Analysis**: How does candidate set structure affect
   parameter identifiability, and when does the FIM become singular?

## File Structure

```
inst/simulations/
|-- README.md                  # This file
|-- sim_utils.R               # Core simulation infrastructure
|-- sim_kl_efficiency.R       # Study 1: KL-divergence and efficiency
|-- sim_misspecification.R    # Study 2: Misspecification bias
|-- sim_identifiability.R     # Study 3: Identifiability via FIM
|-- sim_plots.R               # Publication-quality visualization
|-- run_all.R                 # Master script to run all studies
`-- results/                  # Output directory (created on run)
```

## Quick Start

### Running All Studies

From the command line:

```bash
# Full study (may take several hours)
Rscript run_all.R

# Quick test run (reduced replications)
Rscript run_all.R --quick

# Run specific study
Rscript run_all.R --study=kl --cores=4

# See all options
Rscript run_all.R --help
```

From R:

```r
# Source the framework
source("inst/simulations/sim_utils.R")
source("inst/simulations/sim_kl_efficiency.R")

# Run KL-divergence study
results <- run_kl_efficiency_study(
  rates = c(1, 1.5, 2),
  sample_sizes = c(50, 100, 200),
  kl_levels = c(0, 0.5, 1.0),
  B = 500,
  seed = 42
)

# View results
print_kl_efficiency_results(results)

# Generate plots
source("inst/simulations/sim_plots.R")
plot_rmse_by_kl(results)
```

## Simulation Studies

### Study 1: KL-Divergence Efficiency

**File**: `sim_kl_efficiency.R`

**Research Question**: How does the KL-divergence from non-informative masking
affect MLE efficiency?

**Design**:
- True parameters: lambda = (1, 1.5, 2) (3-component exponential)
- Sample sizes: n in {50, 100, 200, 500}
- KL-divergence levels: d in {0, 0.1, 0.5, 1.0, 2.0}
- Replications: B = 500 per scenario

**Metrics**:
- Bias of each lambda_hat_j
- Variance and MSE
- 95% CI coverage probability
- Relative efficiency vs d=0 baseline

**Key Functions**:
- `run_kl_efficiency_study()`: Main driver
- `print_kl_efficiency_results()`: Console output
- `latex_kl_efficiency_table()`: LaTeX table generation

### Study 2: Misspecification Bias

**File**: `sim_misspecification.R`

**Research Question**: What is the bias when C2 is violated but assumed to hold?

**Design**:
- Generate data with informative masking (alpha, beta parameters)
- Estimate using:
  1. Correct model (oracle with known masking weights)
  2. Misspecified model (assumes C2 holds)
- Compare bias and coverage

**Scenarios**:
- alpha in {0, 1, 5, 10} (uninformative to highly informative)
- beta = 0.3 (fixed)
- n = 200, B = 500

**Interpretation**:
- alpha = 0: Non-informative masking (C2 holds)
- alpha > 0: Informative masking (C2 violated)
- RMSE ratio > 1 indicates misspecification penalty

**Key Functions**:
- `run_misspecification_study()`: Main driver
- `print_misspec_results()`: Console output
- `latex_misspec_table()`: LaTeX table generation

### Study 3: Identifiability Analysis

**File**: `sim_identifiability.R`

**Research Question**: How does candidate set structure affect identifiability?

**Approach**:
- Vary correlation (rho) of component co-occurrence in candidate sets
- Monitor Fisher Information Matrix (FIM) eigenvalues
- Detect when smallest eigenvalue approaches 0 (non-identifiability)

**Design**:
- rho in {0, 0.1, ..., 0.9, 0.99}
- Monte Carlo estimation of expected FIM
- MLE performance analysis at each rho level

**Key Insight**: When components always co-occur (high rho), only their sum
is identifiable, not individual rates. This is demonstrated by the block
model analysis.

**Key Functions**:
- `run_identifiability_study()`: Main driver
- `analyze_block_model()`: Demonstrates non-identifiability
- `fim_eigenvalue_analysis()`: FIM diagnostic
- `print_identifiability_results()`: Console output

## Visualization

**File**: `sim_plots.R`

All plots use a publication-quality theme compatible with academic journals.
Output formats: PDF (vector) and PNG (raster).

**KL-Divergence Plots**:
- `plot_bias_by_kl()`: Bias vs KL-divergence
- `plot_rmse_by_kl()`: RMSE by KL level and sample size
- `plot_efficiency_by_kl()`: Relative efficiency
- `plot_coverage_by_kl()`: CI coverage probability

**Misspecification Plots**:
- `plot_misspec_bias_comparison()`: Correct vs misspecified model
- `plot_misspec_rmse_comparison()`: RMSE by model and alpha
- `plot_misspec_rmse_ratio()`: Misspecification penalty

**Identifiability Plots**:
- `plot_fim_eigenvalues()`: Smallest eigenvalue vs correlation
- `plot_condition_number()`: Condition number diagnostic
- `plot_rmse_by_correlation()`: RMSE degradation with correlation

**Utility Plots**:
- `plot_mse_heatmap()`: MSE across scenarios
- `plot_sampling_distribution()`: Distribution of MLE estimates
- `plot_qq_normal()`: Normality assessment

## Output Structure

After running the full study, results are organized as:

```
results/
|-- all_results.rds              # Combined results object
|-- kl_efficiency/
|   |-- kl_efficiency_raw.rds    # Raw simulation results
|   |-- kl_efficiency_summary.csv
|   |-- kl_efficiency_relative.csv
|   |-- table_rmse.tex           # LaTeX table
|   |-- table_efficiency.tex
|   `-- figures/
|       |-- kl_bias.pdf
|       |-- kl_rmse.pdf
|       |-- kl_efficiency.pdf
|       `-- kl_coverage.pdf
|-- misspecification/
|   |-- misspec_raw.rds
|   |-- misspec_summary.csv
|   |-- misspec_comparison.csv
|   |-- table.tex
|   `-- figures/
|       |-- misspec_bias.pdf
|       |-- misspec_rmse.pdf
|       `-- misspec_ratio.pdf
`-- identifiability/
    |-- identifiability_raw.rds
    |-- identifiability_summary.csv
    |-- block_model.rds
    |-- table.tex
    `-- figures/
        |-- ident_eigenvalues.pdf
        |-- ident_condition.pdf
        |-- ident_rmse.pdf
        `-- ident_combined.pdf
```

## Dependencies

**Required**:
- R >= 4.0
- stats, MASS (base R)
- tibble, dplyr (data manipulation)

**Optional** (for full functionality):
- ggplot2 (visualization)
- gridExtra (combined plots)
- parallel (parallel execution)

**Package Dependencies** (if available):
- md.tools (masked data utilities) - falls back to stubs if not installed
- algebraic.mle (MLE framework) - not required for simulations

## Reproducibility

All simulations are designed for full reproducibility:

1. **Random Seeds**: Every function accepts a `seed` parameter. The master
   script uses deterministic seed derivation for each study.

2. **Intermediate Results**: Results are saved to disk at each stage,
   allowing restart from checkpoints.

3. **Environment Recording**: Consider using `sessionInfo()` to record
   the R environment when running production simulations.

4. **Version Control**: All simulation code is version controlled with
   the package source.

## Performance Notes

- **Memory**: Each scenario stores B estimate vectors. For B=500, m=3,
  this is minimal (~12KB per scenario).

- **Time Estimates** (single core):
  - KL study (full): ~2-4 hours
  - Misspecification study: ~1-2 hours
  - Identifiability study: ~1-2 hours
  - Quick mode: ~10-20 minutes total

- **Parallel Execution**: Use `--cores=N` to parallelize across scenarios.
  Recommended: N = physical cores (not hyperthreads).

## Extending the Framework

To add a new simulation study:

1. Create `sim_newstudy.R` following the template of existing studies
2. Implement:
   - `run_newstudy()`: Main driver function
   - `print_newstudy_results()`: Console output
   - `latex_newstudy_table()`: LaTeX generation (optional)
3. Add plotting functions to `sim_plots.R`
4. Add to `run_all.R` master script

## Citation

If using this simulation framework, please cite:

```
@article{towell2024relaxed,
  title={Relaxed Candidate Set Models for Masked Data in Series Systems},
  author={Towell, Alexander},
  year={2024}
}
```

## Contact

For questions or issues with the simulation framework:
- Alexander Towell <lex@metafunctor.com>
- GitHub: https://github.com/queelius/md_series_system_relaxed_candidate_set_models
