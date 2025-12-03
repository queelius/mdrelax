#!/usr/bin/env Rscript
# =============================================================================
# Master Script: Run All Simulation Studies
# =============================================================================
#
# This script executes all simulation studies for the paper:
# "Relaxed Candidate Set Models for Masked Data in Series Systems"
#
# Studies:
# 1. KL-Divergence Efficiency: Effect of informativeness on MLE efficiency
# 2. Misspecification Bias: Bias when C2 is violated but assumed to hold
# 3. Identifiability Analysis: FIM eigenvalues vs candidate set correlation
#
# Usage:
#   Rscript run_all.R [--quick] [--study=NAME] [--cores=N] [--seed=S]
#
# Options:
#   --quick       Run with reduced replications (for testing)
#   --study=NAME  Run specific study only (kl, misspec, ident, all)
#   --cores=N     Number of cores for parallel execution
#   --seed=S      Random seed for reproducibility
#   --output=DIR  Output directory for results
#
# Author: Alexander Towell
# Package: md_series_system_relaxed_candidate_set_models
# =============================================================================

# -----------------------------------------------------------------------------
# Parse Command Line Arguments
# -----------------------------------------------------------------------------

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    # Defaults
    opts <- list(
        quick = FALSE,
        study = "all",
        cores = 1,
        seed = 42,
        output = NULL
    )

    for (arg in args) {
        if (arg == "--quick") {
            opts$quick <- TRUE
        } else if (grepl("^--study=", arg)) {
            opts$study <- sub("^--study=", "", arg)
        } else if (grepl("^--cores=", arg)) {
            opts$cores <- as.integer(sub("^--cores=", "", arg))
        } else if (grepl("^--seed=", arg)) {
            opts$seed <- as.integer(sub("^--seed=", "", arg))
        } else if (grepl("^--output=", arg)) {
            opts$output <- sub("^--output=", "", arg)
        } else if (arg == "--help" || arg == "-h") {
            cat("Usage: Rscript run_all.R [OPTIONS]\n\n")
            cat("Options:\n")
            cat("  --quick       Run with reduced replications\n")
            cat("  --study=NAME  Run specific study (kl, misspec, ident, all)\n")
            cat("  --cores=N     Number of cores\n")
            cat("  --seed=S      Random seed\n")
            cat("  --output=DIR  Output directory\n")
            quit(save = "no")
        }
    }

    opts
}

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

# Get script directory
script_dir <- tryCatch({
    dirname(sys.frame(1)$ofile)
}, error = function(e) {
    "."
})

# Source simulation files
cat("Loading simulation framework...\n")
source(file.path(script_dir, "sim_utils.R"))
source(file.path(script_dir, "sim_kl_efficiency.R"))
source(file.path(script_dir, "sim_misspecification.R"))
source(file.path(script_dir, "sim_identifiability.R"))
source(file.path(script_dir, "sim_plots.R"))

# Parse arguments
opts <- parse_args()

# Set output directory
if (is.null(opts$output)) {
    opts$output <- file.path(script_dir, "results")
}
dir.create(opts$output, recursive = TRUE, showWarnings = FALSE)

# Print configuration
cat("\n")
cat("=======================================================\n")
cat("Simulation Study Configuration\n")
cat("=======================================================\n")
cat(sprintf("Quick mode: %s\n", opts$quick))
cat(sprintf("Study: %s\n", opts$study))
cat(sprintf("Cores: %d\n", opts$cores))
cat(sprintf("Seed: %d\n", opts$seed))
cat(sprintf("Output: %s\n", opts$output))
cat("=======================================================\n\n")

# Set simulation parameters based on mode
if (opts$quick) {
    # Quick testing parameters
    B_kl <- 100
    B_misspec <- 100
    B_ident <- 100
    n_mc <- 50
    sample_sizes <- c(50, 100)
    kl_levels <- c(0, 0.5, 1.0)
    alpha_levels <- c(0, 5)
    rho_levels <- c(0, 0.5, 0.9)
    cat("Running in QUICK mode (reduced replications)\n\n")
} else {
    # Full study parameters
    B_kl <- 500
    B_misspec <- 500
    B_ident <- 200
    n_mc <- 100
    sample_sizes <- c(50, 100, 200, 500)
    kl_levels <- c(0, 0.1, 0.5, 1.0, 2.0)
    alpha_levels <- c(0, 1, 5, 10)
    rho_levels <- seq(0, 0.99, by = 0.1)
    cat("Running FULL simulation study\n\n")
}

# Common parameters
rates <- c(1, 1.5, 2)  # 3-component exponential series
p_base <- 0.3           # Base masking probability
target_cens <- 0.2      # Target 20% censoring

# -----------------------------------------------------------------------------
# Study 1: KL-Divergence Efficiency
# -----------------------------------------------------------------------------

run_kl_study <- function() {
    cat("\n")
    cat("=======================================================\n")
    cat("STUDY 1: KL-Divergence Efficiency\n")
    cat("=======================================================\n\n")

    start_time <- Sys.time()

    results <- run_kl_efficiency_study(
        rates = rates,
        sample_sizes = sample_sizes,
        kl_levels = kl_levels,
        p = p_base,
        target_cens_prop = target_cens,
        B = B_kl,
        n_cores = opts$cores,
        seed = opts$seed,
        output_dir = file.path(opts$output, "kl_efficiency"),
        progress = TRUE
    )

    # Print summary
    print_kl_efficiency_results(results)

    # Generate plots
    cat("\nGenerating plots...\n")
    plots <- plot_kl_efficiency_all(
        results,
        output_dir = file.path(opts$output, "kl_efficiency", "figures"),
        format = "both"
    )

    # Generate LaTeX tables
    cat("\nGenerating LaTeX tables...\n")
    latex_rmse <- latex_kl_efficiency_table(results, metric = "rmse")
    latex_eff <- latex_kl_efficiency_table(results, metric = "efficiency")

    writeLines(latex_rmse, file.path(opts$output, "kl_efficiency", "table_rmse.tex"))
    writeLines(latex_eff, file.path(opts$output, "kl_efficiency", "table_efficiency.tex"))

    end_time <- Sys.time()
    cat(sprintf("\nStudy 1 completed in %.1f minutes\n",
                difftime(end_time, start_time, units = "mins")))

    results
}

# -----------------------------------------------------------------------------
# Study 2: Misspecification Bias
# -----------------------------------------------------------------------------

run_misspec_study <- function() {
    cat("\n")
    cat("=======================================================\n")
    cat("STUDY 2: Misspecification Bias\n")
    cat("=======================================================\n\n")

    start_time <- Sys.time()

    results <- run_misspecification_study(
        rates = rates,
        n = 200,
        alpha_levels = alpha_levels,
        beta = p_base,
        target_cens_prop = target_cens,
        B = B_misspec,
        n_cores = opts$cores,
        seed = opts$seed + 1000,
        output_dir = file.path(opts$output, "misspecification"),
        progress = TRUE
    )

    # Print summary
    print_misspec_results(results)

    # Generate plots
    cat("\nGenerating plots...\n")
    plots <- plot_misspec_all(
        results,
        output_dir = file.path(opts$output, "misspecification", "figures"),
        format = "both"
    )

    # Generate LaTeX table
    cat("\nGenerating LaTeX table...\n")
    latex_table <- latex_misspec_table(results)
    writeLines(latex_table, file.path(opts$output, "misspecification", "table.tex"))

    end_time <- Sys.time()
    cat(sprintf("\nStudy 2 completed in %.1f minutes\n",
                difftime(end_time, start_time, units = "mins")))

    results
}

# -----------------------------------------------------------------------------
# Study 3: Identifiability Analysis
# -----------------------------------------------------------------------------

run_ident_study <- function() {
    cat("\n")
    cat("=======================================================\n")
    cat("STUDY 3: Identifiability Analysis\n")
    cat("=======================================================\n\n")

    start_time <- Sys.time()

    results <- run_identifiability_study(
        rates = rates,
        n = 200,
        rho_levels = rho_levels,
        p_base = p_base,
        target_cens_prop = target_cens,
        n_mc = n_mc,
        B = B_ident,
        seed = opts$seed + 2000,
        output_dir = file.path(opts$output, "identifiability"),
        progress = TRUE
    )

    # Print summary
    print_identifiability_results(results)

    # Generate plots
    cat("\nGenerating plots...\n")
    plots <- plot_identifiability_all(
        results,
        output_dir = file.path(opts$output, "identifiability", "figures"),
        format = "both"
    )

    # Generate LaTeX table
    cat("\nGenerating LaTeX table...\n")
    latex_table <- latex_identifiability_table(results)
    writeLines(latex_table, file.path(opts$output, "identifiability", "table.tex"))

    # Block model analysis
    cat("\n--- Block Model Analysis ---\n")
    block_results <- analyze_block_model(
        rates = rates,
        n = 200,
        tau = compute_tau_for_censoring(rates, target_cens),
        B = min(B_ident, 200),
        seed = opts$seed + 3000,
        progress = TRUE
    )

    save_results(block_results, file.path(opts$output, "identifiability", "block_model.rds"))

    end_time <- Sys.time()
    cat(sprintf("\nStudy 3 completed in %.1f minutes\n",
                difftime(end_time, start_time, units = "mins")))

    results
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

all_results <- list()
total_start <- Sys.time()

if (opts$study %in% c("all", "kl")) {
    all_results$kl <- run_kl_study()
}

if (opts$study %in% c("all", "misspec")) {
    all_results$misspec <- run_misspec_study()
}

if (opts$study %in% c("all", "ident")) {
    all_results$ident <- run_ident_study()
}

# -----------------------------------------------------------------------------
# Summary Report
# -----------------------------------------------------------------------------

total_end <- Sys.time()

cat("\n")
cat("=======================================================\n")
cat("SIMULATION STUDY COMPLETE\n")
cat("=======================================================\n")
cat(sprintf("Total time: %.1f minutes\n",
            difftime(total_end, total_start, units = "mins")))
cat(sprintf("Results saved to: %s\n", opts$output))
cat("\nOutput files:\n")

# List output files
output_files <- list.files(opts$output, recursive = TRUE, full.names = FALSE)
for (f in head(output_files, 30)) {
    cat(sprintf("  %s\n", f))
}
if (length(output_files) > 30) {
    cat(sprintf("  ... and %d more files\n", length(output_files) - 30))
}

cat("\n")
cat("To view results in R:\n")
cat(sprintf("  results <- readRDS('%s/kl_efficiency/kl_efficiency_raw.rds')\n",
            opts$output))
cat("  print_kl_efficiency_results(results)\n")
cat("\n")

# Save combined results
save_results(all_results, file.path(opts$output, "all_results.rds"))

cat("Done.\n")
