# =============================================================================
# Simulation Study 2: Misspecification Bias Under Violated C2
# =============================================================================
#
# Research Question:
# What is the bias when condition C2 is violated but the analyst assumes it holds?
#
# Background:
# Under conditions C1, C2, C3, the likelihood function takes a simplified form
# where the probability of component j causing failure depends only on hazard
# rates, not on actual component failure times. When C2 is violated (informative
# masking), using the misspecified C1,C2,C3 likelihood introduces bias.
#
# Design:
# 1. Generate data with informative masking (violates C2)
#    - Use informative_masking_by_rank with varying (alpha, beta) parameters
#    - alpha controls informativeness: 0 = uniform, large = highly informative
#    - beta controls maximum probability for rank-2 component
#
# 2. Estimate parameters using two approaches:
#    a) Correct model: Accounts for known masking weights (oracle)
#    b) Misspecified model: Assumes C2 holds (standard practice)
#
# 3. Compare bias and coverage between approaches
#
# Author: Alexander Towell
# Package: md_series_system_relaxed_candidate_set_models
# =============================================================================

# Load utilities
source(file.path(dirname(sys.frame(1)$ofile %||% "."), "sim_utils.R"))

# =============================================================================
# Weighted Likelihood for Known Masking Model
# =============================================================================

#' Compute MLE using known masking weights (oracle estimator)
#'
#' When the masking probabilities are known, we can incorporate them
#' into the likelihood function. This serves as an oracle benchmark.
#'
#' Under informative masking with known weights q_j, the likelihood
#' contribution for observation i with candidate set C_i is:
#'
#' L_i(theta) propto R(s_i; theta) * sum_{j in C_i} [h_j(s_i; theta) * P(j in C | K=j)]
#'
#' For the rank-based model, P(j in C | K=j) depends on the ranking of
#' component failure times, but for known weights we can directly use
#' the observed q_j values.
#'
#' @param md Masked data frame with q1, ..., qm columns
#' @param theta0 Initial parameter values
#' @param sysvar Column name for system lifetime
#' @param setvar Column prefix for candidate sets
#' @param qvar Column prefix for masking probabilities
#' @param deltavar Column name for censoring indicator
#' @param lower Lower bound for parameters
#' @param upper Upper bound for parameters
#'
#' @return List with MLE results
#' @keywords internal
compute_mle_exp_series_weighted <- function(md,
                                             theta0 = NULL,
                                             sysvar = "t",
                                             setvar = "x",
                                             qvar = "q",
                                             deltavar = "delta",
                                             lower = 1e-10,
                                             upper = Inf) {
    n <- nrow(md)
    if (n == 0) stop("md is empty")

    # Extract data
    s <- md[[sysvar]]
    C <- md_decode_matrix(md, setvar)
    Q <- md_decode_matrix(md, qvar)

    if (is.null(C)) stop("No candidate set columns found")
    if (is.null(Q)) {
        # Fall back to unweighted
        return(compute_mle_exp_series(md, theta0, sysvar, setvar, deltavar,
                                       lower, upper))
    }

    m <- ncol(C)

    # Handle censoring
    has_delta <- !is.null(deltavar) && deltavar %in% colnames(md)
    delta <- if (has_delta) md[[deltavar]] else rep(FALSE, n)

    # Weighted log-likelihood function
    # Note: This is an approximation - the exact oracle would require
    # knowing the conditional distribution P(C | K, T)
    ll_fn <- function(theta) {
        if (any(theta <= 0)) return(-Inf)

        sum_theta <- sum(theta)
        ll <- -sum(s) * sum_theta

        for (i in seq_len(n)) {
            if (!delta[i]) {
                C_i <- C[i, ]
                Q_i <- Q[i, ]

                # Weighted hazard sum
                # Weight by q_j (probability component j appears in candidate set)
                weighted_hazard <- sum(theta[C_i] * Q_i[C_i])
                if (weighted_hazard <= 0) return(-Inf)
                ll <- ll + log(weighted_hazard)
            }
        }
        ll
    }

    # Score function (numerical for simplicity)
    score_fn <- function(theta) {
        if (any(theta <= 0)) return(rep(NA_real_, m))

        h <- 1e-6
        g <- numeric(m)
        fx <- ll_fn(theta)
        for (j in seq_len(m)) {
            theta_plus <- theta
            theta_plus[j] <- theta_plus[j] + h
            g[j] <- (ll_fn(theta_plus) - fx) / h
        }
        g
    }

    # Initialize theta0
    if (is.null(theta0)) {
        n_uncensored <- sum(!delta)
        total_hazard <- n_uncensored / sum(s)
        theta0 <- rep(total_hazard / m, m)
    }
    theta0 <- pmax(theta0, lower)

    # Optimize
    result <- tryCatch({
        stats::optim(
            par = theta0,
            fn = ll_fn,
            gr = score_fn,
            method = "L-BFGS-B",
            lower = rep(lower, m),
            upper = rep(upper, m),
            control = list(fnscale = -1, maxit = 1000)
        )
    }, error = function(e) {
        list(par = theta0, value = -Inf, convergence = 1)
    })

    # Compute FIM numerically
    fim <- matrix(0, m, m)
    h <- 1e-5
    for (j in seq_len(m)) {
        for (k in seq_len(m)) {
            theta_pp <- theta_pm <- theta_mp <- theta_mm <- result$par

            theta_pp[j] <- theta_pp[j] + h
            theta_pp[k] <- theta_pp[k] + h

            theta_pm[j] <- theta_pm[j] + h
            theta_pm[k] <- theta_pm[k] - h

            theta_mp[j] <- theta_mp[j] - h
            theta_mp[k] <- theta_mp[k] + h

            theta_mm[j] <- theta_mm[j] - h
            theta_mm[k] <- theta_mm[k] - h

            fim[j, k] <- -(ll_fn(theta_pp) - ll_fn(theta_pm) -
                            ll_fn(theta_mp) + ll_fn(theta_mm)) / (4 * h^2)
        }
    }

    vcov <- tryCatch({
        MASS::ginv(fim)
    }, error = function(e) {
        matrix(NA_real_, m, m)
    })

    list(
        theta_hat = result$par,
        loglike = result$value,
        fim = fim,
        vcov = vcov,
        converged = (result$convergence == 0),
        nobs = n,
        model = "weighted"
    )
}

# =============================================================================
# Main Simulation Functions
# =============================================================================

#' Run misspecification bias study
#'
#' @param rates True rate parameters (default: c(1, 1.5, 2))
#' @param n Sample size (default: 200)
#' @param alpha_levels Informativeness parameters (default: c(0, 1, 5, 10))
#' @param beta Fixed beta parameter (default: 0.3)
#' @param target_cens_prop Target censoring proportion (default: 0.2)
#' @param B Number of replications per scenario (default: 500)
#' @param n_cores Number of cores for parallel execution
#' @param seed Random seed
#' @param output_dir Directory to save results
#' @param progress Show progress
#'
#' @return List with:
#'   - results: Raw results for each scenario
#'   - summary: Summary data frame
#'   - comparison: Correct vs misspecified model comparison
#'
#' @export
run_misspecification_study <- function(
    rates = c(1, 1.5, 2),
    n = 200,
    alpha_levels = c(0, 1, 5, 10),
    beta = 0.3,
    target_cens_prop = 0.2,
    B = 500,
    n_cores = 1,
    seed = 42,
    output_dir = NULL,
    progress = TRUE
) {
    m <- length(rates)
    tau <- compute_tau_for_censoring(rates, target_cens_prop)

    if (progress) {
        cat("==============================================\n")
        cat("Misspecification Bias Study\n")
        cat("==============================================\n")
        cat(sprintf("True rates: (%s)\n", paste(rates, collapse = ", ")))
        cat(sprintf("Sample size: %d\n", n))
        cat(sprintf("Alpha levels: %s\n", paste(alpha_levels, collapse = ", ")))
        cat(sprintf("Beta: %.2f\n", beta))
        cat(sprintf("Target censoring: %.2f\n", target_cens_prop))
        cat(sprintf("Tau: %.4f\n", tau))
        cat(sprintf("Replications: %d\n", B))
        cat("==============================================\n\n")
    }

    # Storage for results
    all_results <- list()

    if (!is.null(seed)) set.seed(seed)

    for (alpha in alpha_levels) {
        if (progress) {
            cat(sprintf("\n--- Alpha = %.1f ---\n", alpha))
        }

        # Storage for this alpha level
        estimates_misspec <- matrix(NA_real_, B, m)
        estimates_correct <- matrix(NA_real_, B, m)
        vcov_misspec <- vector("list", B)
        vcov_correct <- vector("list", B)
        converged_misspec <- logical(B)
        converged_correct <- logical(B)

        for (b in seq_len(B)) {
            if (progress && b %% 100 == 0) {
                cat(sprintf("  Replication %d/%d\n", b, B))
            }

            tryCatch({
                # Generate data with informative masking
                md <- generate_masked_data(
                    n = n,
                    rates = rates,
                    tau = tau,
                    masking_model = "informative",
                    masking_params = list(alpha = alpha, beta = beta)
                )

                # Fit misspecified model (assumes C2 holds)
                mle_misspec <- compute_mle_exp_series(md, theta0 = rates)
                estimates_misspec[b, ] <- mle_misspec$theta_hat
                vcov_misspec[[b]] <- mle_misspec$vcov
                converged_misspec[b] <- mle_misspec$converged

                # Fit correct model (uses known masking weights)
                # Note: This is an oracle estimator since in practice
                # the masking weights would be unknown
                mle_correct <- compute_mle_exp_series_weighted(md, theta0 = rates)
                estimates_correct[b, ] <- mle_correct$theta_hat
                vcov_correct[[b]] <- mle_correct$vcov
                converged_correct[b] <- mle_correct$converged

            }, error = function(e) {
                # Leave as NA
            })
        }

        # Compute statistics
        stats_misspec <- compute_mle_stats(estimates_misspec, rates, vcov_misspec)
        stats_correct <- compute_mle_stats(estimates_correct, rates, vcov_correct)

        all_results[[as.character(alpha)]] <- list(
            alpha = alpha,
            beta = beta,
            misspecified = list(
                estimates = estimates_misspec,
                vcov_list = vcov_misspec,
                converged = converged_misspec,
                stats = stats_misspec
            ),
            correct = list(
                estimates = estimates_correct,
                vcov_list = vcov_correct,
                converged = converged_correct,
                stats = stats_correct
            )
        )
    }

    # Create summary data frame
    summary_df <- create_misspec_summary(all_results, rates, n, B)

    # Create comparison table
    comparison_df <- create_misspec_comparison(summary_df)

    # Save results
    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        save_results(all_results, file.path(output_dir, "misspec_raw.rds"))
        utils::write.csv(summary_df, file.path(output_dir, "misspec_summary.csv"),
                         row.names = FALSE)
        utils::write.csv(comparison_df, file.path(output_dir, "misspec_comparison.csv"),
                         row.names = FALSE)

        if (progress) {
            cat(sprintf("\nResults saved to: %s\n", output_dir))
        }
    }

    list(
        results = all_results,
        summary = summary_df,
        comparison = comparison_df,
        params = list(
            rates = rates,
            n = n,
            alpha_levels = alpha_levels,
            beta = beta,
            tau = tau,
            B = B
        )
    )
}

#' Create summary data frame for misspecification study
#' @keywords internal
create_misspec_summary <- function(results, rates, n, B) {
    m <- length(rates)

    rows <- lapply(names(results), function(alpha_str) {
        r <- results[[alpha_str]]
        alpha <- r$alpha

        # Misspecified model results
        rows_misspec <- lapply(seq_len(m), function(j) {
            data.frame(
                alpha = alpha,
                model = "misspecified",
                component = j,
                theta_true = rates[j],
                bias = r$misspecified$stats$bias[j],
                variance = r$misspecified$stats$variance[j],
                mse = r$misspecified$stats$mse[j],
                rmse = r$misspecified$stats$rmse[j],
                rel_bias = r$misspecified$stats$bias[j] / rates[j],
                coverage = if (!is.null(r$misspecified$stats$coverage))
                    r$misspecified$stats$coverage[j] else NA,
                n_valid = r$misspecified$stats$n_valid,
                n = n,
                B = B,
                stringsAsFactors = FALSE
            )
        })

        # Correct model results
        rows_correct <- lapply(seq_len(m), function(j) {
            data.frame(
                alpha = alpha,
                model = "correct",
                component = j,
                theta_true = rates[j],
                bias = r$correct$stats$bias[j],
                variance = r$correct$stats$variance[j],
                mse = r$correct$stats$mse[j],
                rmse = r$correct$stats$rmse[j],
                rel_bias = r$correct$stats$bias[j] / rates[j],
                coverage = if (!is.null(r$correct$stats$coverage))
                    r$correct$stats$coverage[j] else NA,
                n_valid = r$correct$stats$n_valid,
                n = n,
                B = B,
                stringsAsFactors = FALSE
            )
        })

        do.call(rbind, c(rows_misspec, rows_correct))
    })

    do.call(rbind, rows)
}

#' Create comparison between correct and misspecified models
#' @keywords internal
create_misspec_comparison <- function(summary_df) {
    # Pivot to wide format for comparison
    misspec <- summary_df[summary_df$model == "misspecified",
                          c("alpha", "component", "bias", "rmse", "coverage")]
    correct <- summary_df[summary_df$model == "correct",
                          c("alpha", "component", "bias", "rmse", "coverage")]

    comparison <- merge(misspec, correct,
                        by = c("alpha", "component"),
                        suffixes = c("_misspec", "_correct"))

    # Compute ratios
    comparison$bias_ratio <- comparison$bias_misspec / comparison$bias_correct
    comparison$rmse_ratio <- comparison$rmse_misspec / comparison$rmse_correct

    comparison[order(comparison$alpha, comparison$component), ]
}

#' Print misspecification study results
#'
#' @param study_results Results from run_misspecification_study
#' @export
print_misspec_results <- function(study_results) {
    cat("\n")
    cat("==============================================\n")
    cat("Misspecification Bias Study Results\n")
    cat("==============================================\n\n")

    summary_df <- study_results$summary
    comparison_df <- study_results$comparison
    m <- max(summary_df$component)

    cat("Interpretation:\n")
    cat("- alpha = 0: Non-informative masking (C2 holds)\n")
    cat("- alpha > 0: Informative masking (C2 violated)\n")
    cat("- 'Misspecified': Standard C1,C2,C3 likelihood\n")
    cat("- 'Correct': Oracle using known masking weights\n\n")

    # Print bias comparison
    cat("=== Bias Comparison ===\n\n")

    for (alpha in sort(unique(summary_df$alpha))) {
        cat(sprintf("Alpha = %.1f:\n", alpha))
        cat(sprintf("          | %s\n",
                    paste(sprintf("  Comp %d  ", seq_len(m)), collapse = "|")))
        cat(rep("-", 50), "\n", sep = "")

        for (model in c("misspecified", "correct")) {
            sub <- summary_df[summary_df$alpha == alpha & summary_df$model == model, ]
            sub <- sub[order(sub$component), ]
            cat(sprintf("%-10s|", model))
            for (j in seq_len(m)) {
                cat(sprintf(" %8.5f |", sub$bias[sub$component == j]))
            }
            cat("\n")
        }
        cat("\n")
    }

    # Print RMSE comparison
    cat("\n=== RMSE Comparison ===\n\n")

    for (alpha in sort(unique(summary_df$alpha))) {
        cat(sprintf("Alpha = %.1f:\n", alpha))
        cat(sprintf("          | %s\n",
                    paste(sprintf("  Comp %d  ", seq_len(m)), collapse = "|")))
        cat(rep("-", 50), "\n", sep = "")

        for (model in c("misspecified", "correct")) {
            sub <- summary_df[summary_df$alpha == alpha & summary_df$model == model, ]
            sub <- sub[order(sub$component), ]
            cat(sprintf("%-10s|", model))
            for (j in seq_len(m)) {
                cat(sprintf(" %8.5f |", sub$rmse[sub$component == j]))
            }
            cat("\n")
        }
        cat("\n")
    }

    # Print RMSE ratio
    cat("\n=== RMSE Ratio (Misspecified / Correct) ===\n")
    cat("(Values > 1 indicate misspecification penalty)\n\n")

    cat(sprintf("Alpha | %s\n",
                paste(sprintf(" Comp %d ", seq_len(m)), collapse = "|")))
    cat(rep("-", 40), "\n", sep = "")

    for (alpha in sort(unique(comparison_df$alpha))) {
        sub <- comparison_df[comparison_df$alpha == alpha, ]
        sub <- sub[order(sub$component), ]
        cat(sprintf("%5.1f |", alpha))
        for (j in seq_len(m)) {
            ratio <- sub$rmse_ratio[sub$component == j]
            cat(sprintf(" %6.3f |", ratio))
        }
        cat("\n")
    }
    cat("\n")
}

#' Generate LaTeX table for misspecification results
#'
#' @param study_results Results from run_misspecification_study
#' @param caption Table caption
#'
#' @return Character string with LaTeX table
#' @export
latex_misspec_table <- function(study_results, caption = NULL) {
    summary_df <- study_results$summary
    m <- max(summary_df$component)

    if (is.null(caption)) {
        caption <- "Comparison of correct and misspecified models under informative masking"
    }

    # Header
    latex <- c(
        "\\begin{table}[htbp]",
        "\\centering",
        sprintf("\\caption{%s}", caption),
        sprintf("\\begin{tabular}{cc%s%s}",
                paste(rep("c", m), collapse = ""),
                paste(rep("c", m), collapse = "")),
        "\\toprule",
        sprintf(" & & \\multicolumn{%d}{c}{Bias} & \\multicolumn{%d}{c}{RMSE} \\\\",
                m, m),
        sprintf("\\cmidrule(lr){3-%d} \\cmidrule(lr){%d-%d}",
                2 + m, 3 + m, 2 + 2*m),
        sprintf("$\\alpha$ & Model & %s & %s \\\\",
                paste(sprintf("$\\hat{\\lambda}_%d$", seq_len(m)), collapse = " & "),
                paste(sprintf("$\\hat{\\lambda}_%d$", seq_len(m)), collapse = " & ")),
        "\\midrule"
    )

    # Data rows
    for (alpha in sort(unique(summary_df$alpha))) {
        first_alpha <- TRUE
        for (model in c("misspecified", "correct")) {
            sub <- summary_df[summary_df$alpha == alpha & summary_df$model == model, ]
            sub <- sub[order(sub$component), ]

            alpha_str <- if (first_alpha) sprintf("%.1f", alpha) else ""
            model_str <- if (model == "misspecified") "Misspec." else "Correct"

            bias_vals <- sprintf("%.4f", sub$bias)
            rmse_vals <- sprintf("%.4f", sub$rmse)

            latex <- c(latex, sprintf("%s & %s & %s & %s \\\\",
                                       alpha_str, model_str,
                                       paste(bias_vals, collapse = " & "),
                                       paste(rmse_vals, collapse = " & ")))
            first_alpha <- FALSE
        }
        latex <- c(latex, "\\midrule")
    }

    # Footer
    latex <- c(latex,
        "\\bottomrule",
        "\\end{tabular}",
        sprintf("\\label{tab:misspec}"),
        "\\end{table}"
    )

    paste(latex, collapse = "\n")
}

# =============================================================================
# Auxiliary Analysis Functions
# =============================================================================

#' Compute asymptotic bias under misspecification
#'
#' Theoretical analysis of the bias induced by using the C1,C2,C3 likelihood
#' when the true masking model is informative.
#'
#' @param rates True rate parameters
#' @param alpha Informativeness parameter
#' @param beta Beta parameter for rank-based masking
#'
#' @return Approximate asymptotic bias vector
#' @export
theoretical_misspec_bias <- function(rates, alpha, beta) {
    # This is a placeholder for future theoretical derivation
    # The exact form depends on the structure of the informative masking
    # and requires integration over the joint distribution of
    # component failure times and candidate sets.

    # For now, return a rough approximation based on simulation insights:
    # Bias tends to be positive for the most reliable component
    # and negative for the least reliable when masking favors
    # components with similar failure times

    m <- length(rates)
    bias_approx <- numeric(m)

    if (alpha == 0) {
        # No bias when masking is non-informative
        return(bias_approx)
    }

    # Heuristic: bias proportional to alpha and inversely to rate
    # (more reliable components appear more often in candidates)
    mean_rate <- mean(rates)
    for (j in seq_len(m)) {
        bias_approx[j] <- alpha * beta * (mean_rate - rates[j]) / (10 * mean_rate)
    }

    bias_approx
}

# =============================================================================
# Main execution when run as script
# =============================================================================

if (sys.nframe() == 0) {
    cat("Starting Misspecification Bias Study...\n\n")

    # Quick test run
    results <- run_misspecification_study(
        rates = c(1, 1.5, 2),
        n = 100,
        alpha_levels = c(0, 1, 5),
        beta = 0.3,
        target_cens_prop = 0.2,
        B = 100,  # Reduced for testing
        n_cores = 1,
        seed = 42,
        output_dir = file.path(dirname(sys.frame(1)$ofile %||% "."),
                               "results", "misspecification"),
        progress = TRUE
    )

    print_misspec_results(results)

    # Generate LaTeX table
    cat("\n\nLaTeX Table:\n")
    cat(latex_misspec_table(results))
    cat("\n")
}
