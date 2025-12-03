# =============================================================================
# Simulation Study 1: KL-Divergence and MLE Efficiency
# =============================================================================
#
# Research Question:
# How does the KL-divergence from non-informative masking affect MLE efficiency?
#
# This study examines the relationship between:
# - The degree of informativeness in masking (measured by KL-divergence from
#   the baseline Bernoulli model)
# - The statistical efficiency of maximum likelihood estimators
#
# Under conditions C1, C2, C3, the Bernoulli candidate set model provides
# non-informative masking. When these conditions are relaxed (particularly C2),
# the masking becomes informative and may provide additional information about
# component reliabilities.
#
# Hypothesis:
# - d=0 (non-informative): Baseline efficiency
# - As d increases: Initially may improve efficiency (more information),
#   but extreme informativeness may introduce bias if model is misspecified
#
# Author: Alexander Towell
# Package: md_series_system_relaxed_candidate_set_models
# =============================================================================

# Load utilities
source(file.path(dirname(sys.frame(1)$ofile %||% "."), "sim_utils.R"))

#' Run KL-Divergence Efficiency Study
#'
#' Examines how KL-divergence from non-informative masking affects MLE efficiency.
#'
#' @param rates True rate parameters (default: c(1, 1.5, 2))
#' @param sample_sizes Sample sizes to test (default: c(50, 100, 200, 500))
#' @param kl_levels KL-divergence levels (default: c(0, 0.1, 0.5, 1.0, 2.0))
#' @param p Base masking probability (default: 0.3)
#' @param target_cens_prop Target censoring proportion (default: 0.2)
#' @param B Number of replications per scenario (default: 500)
#' @param n_cores Number of cores for parallel execution (default: 1)
#' @param seed Random seed for reproducibility
#' @param output_dir Directory to save results
#' @param progress Show progress messages
#'
#' @return List with:
#'   - results: Raw simulation results
#'   - summary: Summary data frame
#'   - efficiency: Relative efficiency table
#'
#' @export
run_kl_efficiency_study <- function(
    rates = c(1, 1.5, 2),
    sample_sizes = c(50, 100, 200, 500),
    kl_levels = c(0, 0.1, 0.5, 1.0, 2.0),
    p = 0.3,
    target_cens_prop = 0.2,
    B = 500,
    n_cores = 1,
    seed = 42,
    output_dir = NULL,
    progress = TRUE
) {
    m <- length(rates)

    # Compute censoring time
    tau <- compute_tau_for_censoring(rates, target_cens_prop)

    if (progress) {
        cat("==============================================\n")
        cat("KL-Divergence Efficiency Study\n")
        cat("==============================================\n")
        cat(sprintf("True rates: (%s)\n", paste(rates, collapse = ", ")))
        cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse = ", ")))
        cat(sprintf("KL-divergence levels: %s\n", paste(kl_levels, collapse = ", ")))
        cat(sprintf("Base masking probability p: %.2f\n", p))
        cat(sprintf("Target censoring proportion: %.2f\n", target_cens_prop))
        cat(sprintf("Right-censoring time tau: %.4f\n", tau))
        cat(sprintf("Replications: %d\n", B))
        cat(sprintf("Number of cores: %d\n", n_cores))
        cat("==============================================\n\n")
    }

    # Build scenario list
    scenarios <- list()
    scenario_id <- 0

    for (n in sample_sizes) {
        for (d in kl_levels) {
            scenario_id <- scenario_id + 1

            if (d == 0) {
                # Use standard Bernoulli model
                scenarios[[scenario_id]] <- list(
                    n = n,
                    rates = rates,
                    tau = tau,
                    masking_model = "bernoulli",
                    masking_params = list(p = p),
                    kl_target = d
                )
            } else {
                # Use KL-divergence constrained model
                scenarios[[scenario_id]] <- list(
                    n = n,
                    rates = rates,
                    tau = tau,
                    masking_model = "kl_divergence",
                    masking_params = list(p = p, d = d),
                    kl_target = d
                )
            }
        }
    }

    n_scenarios <- length(scenarios)
    if (progress) {
        cat(sprintf("Total scenarios: %d\n\n", n_scenarios))
    }

    # Run simulations
    results <- run_simulation(
        scenarios = scenarios,
        B = B,
        n_cores = n_cores,
        seed = seed,
        progress = progress
    )

    # Attach scenario metadata to results
    for (i in seq_along(results)) {
        results[[i]]$kl_target <- scenarios[[i]]$kl_target
    }

    # Create summary data frame
    summary_df <- create_kl_efficiency_summary(results, rates)

    # Compute relative efficiencies
    efficiency_df <- compute_relative_efficiencies(summary_df)

    # Save results if output directory specified
    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        save_results(results, file.path(output_dir, "kl_efficiency_raw.rds"))
        utils::write.csv(summary_df, file.path(output_dir, "kl_efficiency_summary.csv"),
                         row.names = FALSE)
        utils::write.csv(efficiency_df, file.path(output_dir, "kl_efficiency_relative.csv"),
                         row.names = FALSE)

        if (progress) {
            cat(sprintf("\nResults saved to: %s\n", output_dir))
        }
    }

    list(
        results = results,
        summary = summary_df,
        efficiency = efficiency_df,
        params = list(
            rates = rates,
            sample_sizes = sample_sizes,
            kl_levels = kl_levels,
            p = p,
            tau = tau,
            B = B
        )
    )
}

#' Create summary data frame for KL efficiency study
#' @keywords internal
create_kl_efficiency_summary <- function(results, rates) {
    m <- length(rates)

    rows <- lapply(seq_along(results), function(i) {
        r <- results[[i]]
        params <- r$params
        stats <- r$stats
        kl_target <- r$kl_target

        # Build row for each component
        rows_i <- lapply(seq_len(m), function(j) {
            data.frame(
                n = params$n,
                kl_target = kl_target,
                component = j,
                theta_true = rates[j],
                bias = stats$bias[j],
                variance = stats$variance[j],
                mse = stats$mse[j],
                rmse = stats$rmse[j],
                rel_bias = stats$bias[j] / rates[j],
                coverage = if (!is.null(stats$coverage)) stats$coverage[j] else NA,
                mean_ci_width = if (!is.null(stats$mean_ci_width))
                    stats$mean_ci_width[j] else NA,
                n_valid = stats$n_valid,
                B = params$B,
                stringsAsFactors = FALSE
            )
        })
        do.call(rbind, rows_i)
    })

    do.call(rbind, rows)
}

#' Compute relative efficiencies vs d=0 baseline
#' @keywords internal
compute_relative_efficiencies <- function(summary_df) {
    # Get baseline (d=0) MSEs
    baseline <- summary_df[summary_df$kl_target == 0, ]

    # Merge with all results
    efficiency_df <- merge(
        summary_df,
        baseline[, c("n", "component", "mse")],
        by = c("n", "component"),
        suffixes = c("", "_baseline")
    )

    # Compute relative efficiency
    efficiency_df$rel_efficiency <- efficiency_df$mse_baseline / efficiency_df$mse

    # Clean up
    efficiency_df$mse_baseline <- NULL

    efficiency_df[order(efficiency_df$n, efficiency_df$kl_target,
                        efficiency_df$component), ]
}

#' Print KL efficiency study results
#'
#' @param study_results Results from run_kl_efficiency_study
#' @export
print_kl_efficiency_results <- function(study_results) {
    cat("\n")
    cat("==============================================\n")
    cat("KL-Divergence Efficiency Study Results\n")
    cat("==============================================\n\n")

    summary_df <- study_results$summary
    efficiency_df <- study_results$efficiency

    # Print by sample size
    for (n in unique(summary_df$n)) {
        cat(sprintf("--- Sample Size n = %d ---\n\n", n))

        cat("Bias by KL-Divergence Level:\n")
        cat("KL-div | ")
        for (j in unique(summary_df$component)) {
            cat(sprintf("Comp %d   | ", j))
        }
        cat("\n")
        cat(rep("-", 50), "\n", sep = "")

        for (d in unique(summary_df$kl_target)) {
            sub <- summary_df[summary_df$n == n & summary_df$kl_target == d, ]
            cat(sprintf("%6.2f | ", d))
            for (j in unique(sub$component)) {
                row <- sub[sub$component == j, ]
                cat(sprintf("%8.5f | ", row$bias))
            }
            cat("\n")
        }
        cat("\n")

        cat("RMSE by KL-Divergence Level:\n")
        cat("KL-div | ")
        for (j in unique(summary_df$component)) {
            cat(sprintf("Comp %d   | ", j))
        }
        cat("\n")
        cat(rep("-", 50), "\n", sep = "")

        for (d in unique(summary_df$kl_target)) {
            sub <- summary_df[summary_df$n == n & summary_df$kl_target == d, ]
            cat(sprintf("%6.2f | ", d))
            for (j in unique(sub$component)) {
                row <- sub[sub$component == j, ]
                cat(sprintf("%8.5f | ", row$rmse))
            }
            cat("\n")
        }
        cat("\n")

        cat("Relative Efficiency (vs d=0 baseline):\n")
        cat("KL-div | ")
        for (j in unique(efficiency_df$component)) {
            cat(sprintf("Comp %d   | ", j))
        }
        cat("\n")
        cat(rep("-", 50), "\n", sep = "")

        for (d in unique(efficiency_df$kl_target)) {
            sub <- efficiency_df[efficiency_df$n == n & efficiency_df$kl_target == d, ]
            cat(sprintf("%6.2f | ", d))
            for (j in unique(sub$component)) {
                row <- sub[sub$component == j, ]
                cat(sprintf("%8.3f | ", row$rel_efficiency))
            }
            cat("\n")
        }
        cat("\n")

        if (!is.null(summary_df$coverage) && !all(is.na(summary_df$coverage))) {
            cat("95% CI Coverage Probability:\n")
            cat("KL-div | ")
            for (j in unique(summary_df$component)) {
                cat(sprintf("Comp %d   | ", j))
            }
            cat("\n")
            cat(rep("-", 50), "\n", sep = "")

            for (d in unique(summary_df$kl_target)) {
                sub <- summary_df[summary_df$n == n & summary_df$kl_target == d, ]
                cat(sprintf("%6.2f | ", d))
                for (j in unique(sub$component)) {
                    row <- sub[sub$component == j, ]
                    cat(sprintf("%8.3f | ", row$coverage))
                }
                cat("\n")
            }
            cat("\n")
        }
    }
}

#' Generate LaTeX table for KL efficiency results
#'
#' @param study_results Results from run_kl_efficiency_study
#' @param metric One of "bias", "rmse", "coverage", "efficiency"
#' @param caption Table caption
#'
#' @return Character string with LaTeX table
#' @export
latex_kl_efficiency_table <- function(study_results,
                                       metric = "rmse",
                                       caption = NULL) {
    summary_df <- study_results$summary
    m <- max(summary_df$component)

    if (is.null(caption)) {
        caption <- switch(metric,
            bias = "Bias of MLE by KL-divergence level and sample size",
            rmse = "RMSE of MLE by KL-divergence level and sample size",
            coverage = "95\\% CI coverage probability by KL-divergence level",
            efficiency = "Relative efficiency vs.~non-informative masking (d=0)"
        )
    }

    # Header
    latex <- c(
        "\\begin{table}[htbp]",
        "\\centering",
        sprintf("\\caption{%s}", caption),
        sprintf("\\begin{tabular}{cc%s}", paste(rep("c", m), collapse = "")),
        "\\toprule",
        sprintf("$n$ & $d$ & %s \\\\",
                paste(sprintf("$\\hat{\\lambda}_%d$", seq_len(m)), collapse = " & ")),
        "\\midrule"
    )

    # Data rows
    for (n in sort(unique(summary_df$n))) {
        first_n <- TRUE
        for (d in sort(unique(summary_df$kl_target))) {
            sub <- summary_df[summary_df$n == n & summary_df$kl_target == d, ]
            sub <- sub[order(sub$component), ]

            if (metric == "efficiency" && d == 0) {
                # Skip baseline for efficiency
                vals <- rep("1.000", m)
            } else if (metric == "efficiency") {
                eff_df <- study_results$efficiency
                eff_sub <- eff_df[eff_df$n == n & eff_df$kl_target == d, ]
                eff_sub <- eff_sub[order(eff_sub$component), ]
                vals <- sprintf("%.3f", eff_sub$rel_efficiency)
            } else {
                vals <- switch(metric,
                    bias = sprintf("%.4f", sub$bias),
                    rmse = sprintf("%.4f", sub$rmse),
                    coverage = sprintf("%.3f", sub$coverage)
                )
            }

            n_str <- if (first_n) as.character(n) else ""
            latex <- c(latex, sprintf("%s & %.1f & %s \\\\",
                                       n_str, d, paste(vals, collapse = " & ")))
            first_n <- FALSE
        }
        latex <- c(latex, "\\midrule")
    }

    # Footer
    latex <- c(latex,
        "\\bottomrule",
        "\\end{tabular}",
        "\\end{table}"
    )

    paste(latex, collapse = "\n")
}

# =============================================================================
# Main execution when run as script
# =============================================================================

if (sys.nframe() == 0) {
    # Running as standalone script

    cat("Starting KL-Divergence Efficiency Study...\n\n")

    # Quick test run with reduced parameters
    results <- run_kl_efficiency_study(
        rates = c(1, 1.5, 2),
        sample_sizes = c(50, 100),
        kl_levels = c(0, 0.5, 1.0),
        p = 0.3,
        target_cens_prop = 0.2,
        B = 100,  # Reduced for testing
        n_cores = 1,
        seed = 42,
        output_dir = file.path(dirname(sys.frame(1)$ofile %||% "."),
                               "results", "kl_efficiency"),
        progress = TRUE
    )

    print_kl_efficiency_results(results)

    # Generate LaTeX table
    cat("\n\nLaTeX Table (RMSE):\n")
    cat(latex_kl_efficiency_table(results, metric = "rmse"))
    cat("\n")
}
