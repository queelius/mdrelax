# =============================================================================
# Simulation Study 3: Identifiability Analysis via FIM Eigenvalues
# =============================================================================
#
# Research Question:
# How does candidate set structure affect parameter identifiability?
#
# Background:
# Parameter identifiability in masked data analysis depends critically on
# the structure of candidate sets. When components always co-occur in
# candidate sets (high correlation), it becomes impossible to distinguish
# their individual failure rates -- only their sum is identifiable.
#
# This is demonstrated by the block candidate model (md_block_candidate_m3)
# where components 1 and 2 always appear together: only lambda_1 + lambda_2
# is identifiable, not the individual rates.
#
# Approach:
# 1. Vary the correlation of component co-occurrence in candidate sets
# 2. Monitor Fisher Information Matrix (FIM) eigenvalues
# 3. Detect when smallest eigenvalue approaches 0 (non-identifiability)
# 4. Characterize the "identifiability boundary"
#
# Key Insight:
# The FIM is positive semi-definite. When the smallest eigenvalue is 0,
# the model has a flat direction in parameter space -- infinitely many
# parameter combinations give the same likelihood.
#
# Author: Alexander Towell
# Package: md_series_system_relaxed_candidate_set_models
# =============================================================================

# Load utilities
source(file.path(dirname(sys.frame(1)$ofile %||% "."), "sim_utils.R"))

# =============================================================================
# Candidate Set Structure Models
# =============================================================================

#' Generate candidate sets with controlled co-occurrence correlation
#'
#' Creates candidate sets where the correlation between component j and k
#' appearing together is controlled by parameter rho.
#'
#' @param n Number of observations
#' @param m Number of components
#' @param k_failed Vector of failed component indices
#' @param p_base Base probability of inclusion
#' @param rho Correlation parameter for component co-occurrence (0 to 1)
#'            rho = 0: independent inclusion
#'            rho = 1: components always co-occur (block structure)
#' @param blocks Optional list of component groups that co-occur
#'
#' @return n x m Boolean matrix of candidate sets
#' @keywords internal
generate_correlated_candidates <- function(n, m, k_failed, p_base, rho,
                                            blocks = NULL) {
    C <- matrix(FALSE, nrow = n, ncol = m)

    # Failed component always included
    for (i in seq_len(n)) {
        C[i, k_failed[i]] <- TRUE
    }

    if (is.null(blocks)) {
        # Default: all non-failed components form one correlated group
        for (i in seq_len(n)) {
            # Draw a common random variable for correlation
            u_common <- stats::runif(1)

            for (j in seq_len(m)) {
                if (j == k_failed[i]) next  # Already included

                # Mix common and independent random variables
                u_indep <- stats::runif(1)
                u_mix <- rho * u_common + (1 - rho) * u_indep

                C[i, j] <- (u_mix <= p_base)
            }
        }
    } else {
        # Block structure: components in same block have same inclusion
        for (i in seq_len(n)) {
            # Process each block
            for (block in blocks) {
                block_comps <- block[block != k_failed[i]]
                if (length(block_comps) == 0) next

                # Single decision for the block
                u_block <- stats::runif(1)
                include_block <- (u_block <= p_base)

                # Apply correlation within block
                u_common <- stats::runif(1)
                for (j in block_comps) {
                    u_indep <- stats::runif(1)
                    u_mix <- rho * u_common + (1 - rho) * u_indep

                    if (rho >= 0.99) {
                        C[i, j] <- include_block
                    } else {
                        C[i, j] <- (u_mix <= p_base)
                    }
                }
            }
        }
    }

    C
}

#' Generate masked data with controlled candidate set correlation
#'
#' @param n Sample size
#' @param rates True rate parameters
#' @param tau Right-censoring time
#' @param p_base Base probability
#' @param rho Correlation parameter
#' @param blocks Component blocks (optional)
#' @param seed Random seed
#'
#' @return Masked data frame
#' @keywords internal
generate_correlated_masked_data <- function(n, rates, tau, p_base, rho,
                                             blocks = NULL, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    m <- length(rates)

    # Generate component times
    Tm <- generate_component_times(n, rates)

    # System lifetime and failed component
    sys_time <- apply(Tm, 1, min)
    failed_comp <- apply(Tm, 1, which.min)

    # Apply censoring
    tau_vec <- rep(tau, length.out = n)
    delta <- sys_time > tau_vec
    obs_time <- pmin(sys_time, tau_vec)

    # Generate correlated candidate sets
    C <- generate_correlated_candidates(n, m, failed_comp, p_base, rho, blocks)

    # Set to FALSE for censored observations
    C[delta, ] <- FALSE

    # Build data frame
    md <- tibble::tibble(
        t = obs_time,
        k = failed_comp,
        delta = delta
    )
    md <- dplyr::bind_cols(md, md_encode_matrix(Tm, "t"))
    md <- dplyr::bind_cols(md, md_encode_matrix(C, "x"))
    md <- md_mark_latent(md, c(paste0("t", seq_len(m)), "k"))

    md
}

# =============================================================================
# FIM Analysis Functions
# =============================================================================

#' Compute expected FIM via Monte Carlo averaging
#'
#' @param n Sample size
#' @param rates True parameters
#' @param tau Censoring time
#' @param p_base Base masking probability
#' @param rho Correlation parameter
#' @param blocks Component blocks
#' @param n_mc Number of MC samples for averaging
#' @param seed Random seed
#'
#' @return List with expected FIM and eigenvalue analysis
#' @keywords internal
compute_expected_fim_correlated <- function(n, rates, tau, p_base, rho,
                                             blocks = NULL,
                                             n_mc = 50,
                                             seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    m <- length(rates)
    fim_sum <- matrix(0, m, m)
    n_valid <- 0

    for (i in seq_len(n_mc)) {
        md <- generate_correlated_masked_data(
            n = n, rates = rates, tau = tau,
            p_base = p_base, rho = rho, blocks = blocks
        )

        # Compute FIM at true parameters
        C <- md_decode_matrix(md, "x")
        s <- md$t
        delta <- md$delta

        fim <- matrix(0, m, m)
        for (k in seq_len(n)) {
            if (!delta[k]) {
                C_k <- C[k, ]
                hazard_sum <- sum(rates[C_k])
                if (hazard_sum > 0) {
                    inv_sq <- 1 / hazard_sum^2
                    for (j1 in seq_len(m)) {
                        if (C_k[j1]) {
                            for (j2 in seq_len(m)) {
                                if (C_k[j2]) {
                                    fim[j1, j2] <- fim[j1, j2] + inv_sq
                                }
                            }
                        }
                    }
                }
            }
        }

        if (all(is.finite(fim))) {
            fim_sum <- fim_sum + fim
            n_valid <- n_valid + 1
        }
    }

    if (n_valid == 0) {
        return(list(
            fim = matrix(NA, m, m),
            eigenvalues = rep(NA, m),
            condition_number = NA,
            smallest_eigenvalue = NA,
            near_singular = TRUE
        ))
    }

    fim_avg <- fim_sum / n_valid
    eig_analysis <- fim_eigenvalue_analysis(fim_avg)

    list(
        fim = fim_avg,
        eigenvalues = eig_analysis$eigenvalues,
        condition_number = eig_analysis$condition_number,
        smallest_eigenvalue = eig_analysis$smallest,
        near_singular = eig_analysis$near_singular,
        n_valid = n_valid
    )
}

# =============================================================================
# Main Simulation Functions
# =============================================================================

#' Run identifiability analysis study
#'
#' @param rates True rate parameters (default: c(1, 1.5, 2))
#' @param n Sample size (default: 200)
#' @param rho_levels Correlation levels to test (default: seq(0, 0.99, 0.1))
#' @param p_base Base masking probability (default: 0.3)
#' @param target_cens_prop Target censoring proportion (default: 0.2)
#' @param n_mc Number of MC samples for FIM estimation (default: 100)
#' @param B Number of replications for MLE analysis (default: 200)
#' @param seed Random seed
#' @param output_dir Directory to save results
#' @param progress Show progress
#'
#' @return List with:
#'   - fim_results: FIM eigenvalue analysis by rho
#'   - mle_results: MLE performance by rho
#'   - summary: Summary data frame
#'
#' @export
run_identifiability_study <- function(
    rates = c(1, 1.5, 2),
    n = 200,
    rho_levels = seq(0, 0.99, by = 0.1),
    p_base = 0.3,
    target_cens_prop = 0.2,
    n_mc = 100,
    B = 200,
    seed = 42,
    output_dir = NULL,
    progress = TRUE
) {
    m <- length(rates)
    tau <- compute_tau_for_censoring(rates, target_cens_prop)

    if (progress) {
        cat("==============================================\n")
        cat("Identifiability Analysis Study\n")
        cat("==============================================\n")
        cat(sprintf("True rates: (%s)\n", paste(rates, collapse = ", ")))
        cat(sprintf("Sample size: %d\n", n))
        cat(sprintf("Rho levels: %s\n",
                    paste(sprintf("%.2f", rho_levels), collapse = ", ")))
        cat(sprintf("Base masking probability: %.2f\n", p_base))
        cat(sprintf("Target censoring: %.2f\n", target_cens_prop))
        cat(sprintf("MC samples for FIM: %d\n", n_mc))
        cat(sprintf("MLE replications: %d\n", B))
        cat("==============================================\n\n")
    }

    if (!is.null(seed)) set.seed(seed)

    # Storage
    fim_results <- list()
    mle_results <- list()

    for (rho in rho_levels) {
        if (progress) {
            cat(sprintf("\n--- Rho = %.2f ---\n", rho))
        }

        # Part 1: FIM eigenvalue analysis
        if (progress) cat("  Computing expected FIM...\n")

        fim_analysis <- compute_expected_fim_correlated(
            n = n, rates = rates, tau = tau,
            p_base = p_base, rho = rho,
            n_mc = n_mc,
            seed = NULL  # Use current RNG state
        )

        fim_results[[as.character(rho)]] <- fim_analysis

        if (progress) {
            cat(sprintf("  Smallest eigenvalue: %.6f\n", fim_analysis$smallest_eigenvalue))
            cat(sprintf("  Condition number: %.2f\n", fim_analysis$condition_number))
        }

        # Part 2: MLE performance analysis
        if (progress) cat(sprintf("  Running %d MLE replications...\n", B))

        estimates <- matrix(NA_real_, B, m)
        vcov_list <- vector("list", B)
        converged <- logical(B)

        for (b in seq_len(B)) {
            if (progress && b %% 50 == 0) {
                cat(sprintf("    Replication %d/%d\n", b, B))
            }

            tryCatch({
                md <- generate_correlated_masked_data(
                    n = n, rates = rates, tau = tau,
                    p_base = p_base, rho = rho
                )

                mle_result <- compute_mle_exp_series(md, theta0 = rates)

                estimates[b, ] <- mle_result$theta_hat
                vcov_list[[b]] <- mle_result$vcov
                converged[b] <- mle_result$converged
            }, error = function(e) {
                # Leave as NA
            })
        }

        stats <- compute_mle_stats(estimates, rates, vcov_list)

        mle_results[[as.character(rho)]] <- list(
            estimates = estimates,
            vcov_list = vcov_list,
            converged = converged,
            stats = stats
        )

        if (progress) {
            cat(sprintf("  Valid estimates: %d/%d\n", stats$n_valid, B))
            cat(sprintf("  RMSE: (%s)\n",
                        paste(sprintf("%.4f", stats$rmse), collapse = ", ")))
        }
    }

    # Create summary data frame
    summary_df <- create_identifiability_summary(
        fim_results, mle_results, rates, rho_levels, n, B
    )

    # Save results
    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        save_results(list(fim = fim_results, mle = mle_results),
                     file.path(output_dir, "identifiability_raw.rds"))
        utils::write.csv(summary_df, file.path(output_dir, "identifiability_summary.csv"),
                         row.names = FALSE)

        if (progress) {
            cat(sprintf("\nResults saved to: %s\n", output_dir))
        }
    }

    list(
        fim_results = fim_results,
        mle_results = mle_results,
        summary = summary_df,
        params = list(
            rates = rates,
            n = n,
            rho_levels = rho_levels,
            p_base = p_base,
            tau = tau,
            n_mc = n_mc,
            B = B
        )
    )
}

#' Create summary data frame for identifiability study
#' @keywords internal
create_identifiability_summary <- function(fim_results, mle_results,
                                            rates, rho_levels, n, B) {
    m <- length(rates)

    rows <- lapply(as.character(rho_levels), function(rho_str) {
        rho <- as.numeric(rho_str)
        fim <- fim_results[[rho_str]]
        mle <- mle_results[[rho_str]]

        # FIM-level row
        fim_row <- data.frame(
            rho = rho,
            smallest_eigenvalue = fim$smallest_eigenvalue,
            condition_number = fim$condition_number,
            near_singular = fim$near_singular,
            n_valid_mle = mle$stats$n_valid,
            stringsAsFactors = FALSE
        )

        # Component-level rows
        comp_rows <- lapply(seq_len(m), function(j) {
            data.frame(
                rho = rho,
                component = j,
                theta_true = rates[j],
                bias = mle$stats$bias[j],
                variance = mle$stats$variance[j],
                mse = mle$stats$mse[j],
                rmse = mle$stats$rmse[j],
                coverage = if (!is.null(mle$stats$coverage))
                    mle$stats$coverage[j] else NA,
                stringsAsFactors = FALSE
            )
        })

        # Combine
        comp_df <- do.call(rbind, comp_rows)
        comp_df$smallest_eigenvalue <- fim$smallest_eigenvalue
        comp_df$condition_number <- fim$condition_number
        comp_df$near_singular <- fim$near_singular
        comp_df$n <- n
        comp_df$B <- B

        comp_df
    })

    do.call(rbind, rows)
}

#' Print identifiability study results
#'
#' @param study_results Results from run_identifiability_study
#' @export
print_identifiability_results <- function(study_results) {
    cat("\n")
    cat("==============================================\n")
    cat("Identifiability Analysis Results\n")
    cat("==============================================\n\n")

    summary_df <- study_results$summary
    m <- max(summary_df$component)
    rho_levels <- sort(unique(summary_df$rho))

    # Print FIM eigenvalue analysis
    cat("=== FIM Eigenvalue Analysis ===\n\n")
    cat("  Rho  | Smallest EV | Cond. Number | Near Singular\n")
    cat(rep("-", 55), "\n", sep = "")

    for (rho in rho_levels) {
        sub <- summary_df[summary_df$rho == rho & summary_df$component == 1, ]
        near_sing <- if (sub$near_singular) "YES" else "no"
        cat(sprintf(" %5.2f | %11.6f | %12.2f | %s\n",
                    rho, sub$smallest_eigenvalue, sub$condition_number, near_sing))
    }
    cat("\n")

    # Print RMSE by correlation
    cat("\n=== RMSE by Correlation Level ===\n\n")
    cat("  Rho  |")
    for (j in seq_len(m)) cat(sprintf(" Comp %d  |", j))
    cat("\n")
    cat(rep("-", 10 + 10*m), "\n", sep = "")

    for (rho in rho_levels) {
        sub <- summary_df[summary_df$rho == rho, ]
        sub <- sub[order(sub$component), ]
        cat(sprintf(" %5.2f |", rho))
        for (j in seq_len(m)) {
            cat(sprintf(" %7.4f |", sub$rmse[sub$component == j]))
        }
        cat("\n")
    }
    cat("\n")

    # Identify critical correlation threshold
    threshold_idx <- which(summary_df$near_singular & summary_df$component == 1)[1]
    if (!is.na(threshold_idx)) {
        threshold_rho <- summary_df$rho[threshold_idx]
        cat(sprintf("Critical correlation threshold: rho >= %.2f\n", threshold_rho))
        cat("(FIM becomes near-singular, parameters not fully identifiable)\n\n")
    }
}

#' Generate LaTeX table for identifiability results
#'
#' @param study_results Results from run_identifiability_study
#' @param caption Table caption
#'
#' @return Character string with LaTeX table
#' @export
latex_identifiability_table <- function(study_results, caption = NULL) {
    summary_df <- study_results$summary
    m <- max(summary_df$component)

    if (is.null(caption)) {
        caption <- "FIM eigenvalue analysis and MLE performance by correlation level"
    }

    latex <- c(
        "\\begin{table}[htbp]",
        "\\centering",
        sprintf("\\caption{%s}", caption),
        sprintf("\\begin{tabular}{c|cc|%s}", paste(rep("c", m), collapse = "")),
        "\\toprule",
        sprintf("$\\rho$ & $\\lambda_{\\min}$ & Cond. & %s \\\\",
                paste(sprintf("RMSE($\\hat{\\lambda}_%d$)", seq_len(m)),
                      collapse = " & ")),
        "\\midrule"
    )

    rho_levels <- sort(unique(summary_df$rho))

    for (rho in rho_levels) {
        sub <- summary_df[summary_df$rho == rho, ]
        sub <- sub[order(sub$component), ]

        eig <- sub$smallest_eigenvalue[1]
        cond <- sub$condition_number[1]
        rmse_vals <- sprintf("%.4f", sub$rmse)

        # Format eigenvalue with warning for near-singular
        eig_str <- sprintf("%.4f", eig)
        if (sub$near_singular[1]) {
            eig_str <- paste0("\\textbf{", eig_str, "}")
        }

        latex <- c(latex, sprintf("%.2f & %s & %.1f & %s \\\\",
                                   rho, eig_str, cond,
                                   paste(rmse_vals, collapse = " & ")))
    }

    latex <- c(latex,
        "\\bottomrule",
        "\\end{tabular}",
        "\\label{tab:identifiability}",
        "\\end{table}"
    )

    paste(latex, collapse = "\n")
}

# =============================================================================
# Block Model Analysis
# =============================================================================

#' Analyze block candidate model identifiability
#'
#' Demonstrates non-identifiability using the md_block_candidate_m3 model
#' where components 1 and 2 always co-occur.
#'
#' @param rates True rate parameters (must have length 3)
#' @param n Sample size
#' @param tau Censoring time
#' @param B Number of replications
#' @param seed Random seed
#' @param progress Show progress
#'
#' @return List with estimates showing only lambda_1 + lambda_2 is identifiable
#' @export
analyze_block_model <- function(rates = c(1, 1.5, 2),
                                 n = 200,
                                 tau = 0.5,
                                 B = 500,
                                 seed = 42,
                                 progress = TRUE) {
    stopifnot(length(rates) == 3)

    if (progress) {
        cat("==============================================\n")
        cat("Block Model Identifiability Analysis\n")
        cat("==============================================\n")
        cat("In this model, components 1 and 2 always co-occur.\n")
        cat("Only lambda_1 + lambda_2 is identifiable.\n\n")
    }

    if (!is.null(seed)) set.seed(seed)

    estimates <- matrix(NA_real_, B, 3)

    for (b in seq_len(B)) {
        if (progress && b %% 100 == 0) {
            cat(sprintf("Replication %d/%d\n", b, B))
        }

        tryCatch({
            # Generate component times
            Tm <- generate_component_times(n, rates)

            # System lifetime
            sys_time <- apply(Tm, 1, min)
            failed_comp <- apply(Tm, 1, which.min)

            # Censoring
            delta <- sys_time > tau
            obs_time <- pmin(sys_time, tau)

            # Block candidate model
            C <- matrix(FALSE, n, 3)
            for (i in seq_len(n)) {
                if (delta[i]) next  # No candidates for censored

                k <- failed_comp[i]
                if (k == 1 || k == 2) {
                    C[i, 1] <- TRUE
                    C[i, 2] <- TRUE
                } else {  # k == 3
                    if (stats::runif(1) < 0.1) {
                        C[i, ] <- TRUE  # All three
                    } else {
                        C[i, 3] <- TRUE  # Only 3
                    }
                }
            }

            # Build data frame
            md <- tibble::tibble(
                t = obs_time,
                k = failed_comp,
                delta = delta
            )
            md <- dplyr::bind_cols(md, md_encode_matrix(Tm, "t"))
            md <- dplyr::bind_cols(md, md_encode_matrix(C, "x"))

            # Fit MLE
            mle_result <- compute_mle_exp_series(md, theta0 = rates)
            estimates[b, ] <- mle_result$theta_hat

        }, error = function(e) {
            # Leave as NA
        })
    }

    # Compute sum of lambda_1 + lambda_2
    sum_12 <- estimates[, 1] + estimates[, 2]
    valid <- !is.na(sum_12)

    results <- list(
        estimates = estimates[valid, ],
        sum_12 = sum_12[valid],
        true_sum_12 = rates[1] + rates[2],
        true_lambda_3 = rates[3],
        mean_sum_12 = mean(sum_12, na.rm = TRUE),
        sd_sum_12 = stats::sd(sum_12, na.rm = TRUE),
        mean_lambda_3 = mean(estimates[valid, 3]),
        sd_lambda_3 = stats::sd(estimates[valid, 3]),
        mean_lambda_1 = mean(estimates[valid, 1]),
        sd_lambda_1 = stats::sd(estimates[valid, 1]),
        mean_lambda_2 = mean(estimates[valid, 2]),
        sd_lambda_2 = stats::sd(estimates[valid, 2])
    )

    if (progress) {
        cat("\n=== Results ===\n")
        cat(sprintf("True lambda_1 + lambda_2: %.2f\n", results$true_sum_12))
        cat(sprintf("Mean estimated sum: %.4f (SD: %.4f)\n",
                    results$mean_sum_12, results$sd_sum_12))
        cat(sprintf("\nTrue lambda_3: %.2f\n", results$true_lambda_3))
        cat(sprintf("Mean estimated lambda_3: %.4f (SD: %.4f)\n",
                    results$mean_lambda_3, results$sd_lambda_3))
        cat(sprintf("\nIndividual estimates (NOT identifiable):\n"))
        cat(sprintf("Mean lambda_1: %.4f (SD: %.4f)\n",
                    results$mean_lambda_1, results$sd_lambda_1))
        cat(sprintf("Mean lambda_2: %.4f (SD: %.4f)\n",
                    results$mean_lambda_2, results$sd_lambda_2))
        cat("\nNote: Large SD for individual lambdas confirms non-identifiability.\n")
        cat("The sum lambda_1 + lambda_2 has much smaller relative SD.\n")
    }

    results
}

# =============================================================================
# Main execution when run as script
# =============================================================================

if (sys.nframe() == 0) {
    cat("Starting Identifiability Analysis Study...\n\n")

    # Quick test run
    results <- run_identifiability_study(
        rates = c(1, 1.5, 2),
        n = 100,
        rho_levels = c(0, 0.3, 0.6, 0.9, 0.99),
        p_base = 0.3,
        target_cens_prop = 0.2,
        n_mc = 50,
        B = 100,  # Reduced for testing
        seed = 42,
        output_dir = file.path(dirname(sys.frame(1)$ofile %||% "."),
                               "results", "identifiability"),
        progress = TRUE
    )

    print_identifiability_results(results)

    # Block model analysis
    cat("\n\n")
    block_results <- analyze_block_model(
        rates = c(1, 1.5, 2),
        n = 100,
        tau = 0.5,
        B = 100,
        progress = TRUE
    )

    # Generate LaTeX table
    cat("\n\nLaTeX Table:\n")
    cat(latex_identifiability_table(results))
    cat("\n")
}
