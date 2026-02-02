# =============================================================================
# Simulation Study: Comparing Relaxed Candidate Set Models
# =============================================================================
#
# This file provides a framework for simulation studies comparing:
#   - C1-C2-C3 model (standard)
#   - Relaxed C2 model (general Bernoulli with P matrix)
#   - Relaxed C3 model (parameter-dependent masking)
#
# Key questions:
#   1. What is the cost of fitting relaxed C2 when C2 actually holds?
#   2. What is the bias when fitting C1-C2-C3 but C2 is violated?
#   3. What happens when we misspecify the C3 relaxation?
#
# =============================================================================

# -----------------------------------------------------------------------------
# Simulation Configuration
# -----------------------------------------------------------------------------

#' Create a simulation configuration
#'
#' @param n_sim Number of simulation replications
#' @param n_obs Number of observations per replication
#' @param theta True rate parameters
#' @param tau Right-censoring time
#' @param seed Random seed for reproducibility
#' @return Configuration list
#' @export
sim_config <- function(n_sim = 100, n_obs = 200, theta = c(1, 2),
                       tau = 5, seed = 42) {
    list(
        n_sim = n_sim,
        n_obs = n_obs,
        theta = theta,
        m = length(theta),
        tau = tau,
        seed = seed
    )
}

# -----------------------------------------------------------------------------
# Data Generation Scenarios
# -----------------------------------------------------------------------------

#' Generate data under C1-C2-C3 (uniform Bernoulli)
#'
#' @param config Simulation configuration
#' @param p Inclusion probability for non-failed components
#' @return List with t, C, delta, k, and true parameters
#' @keywords internal
generate_c1_c2_c3 <- function(config, p = 0.3) {
    sim <- rexp_series_md(
        n = config$n_obs,
        theta = config$theta,
        p = p,
        tau = config$tau
    )
    sim$true_theta <- config$theta
    sim$true_p <- p
    sim$scenario <- "C1-C2-C3"
    sim
}

#' Generate data under relaxed C2 (general Bernoulli)
#'
#' @param config Simulation configuration
#' @param P Inclusion probability matrix (m x m)
#' @return List with t, C, delta, k, and true parameters
#' @keywords internal
generate_relaxed_c2 <- function(config, P) {
    sim <- rexp_series_md_c1_c3(
        n = config$n_obs,
        theta = config$theta,
        P = P,
        tau = config$tau
    )
    sim$true_theta <- config$theta
    sim$true_P <- P
    sim$scenario <- "Relaxed-C2"
    sim
}

#' Generate data under relaxed C3 (power-weighted masking)
#'
#' This generates data where inclusion probability depends on hazard rates.
#' For exponential components: p_j ∝ theta_j^alpha
#'
#' @param config Simulation configuration
#' @param alpha Power parameter (0 = uniform, higher = more informative)
#' @param base_p Baseline inclusion probability
#' @return List with t, C, delta, k, and true parameters
#' @keywords internal
generate_relaxed_c3 <- function(config, alpha = 1, base_p = 0.5) {
    n <- config$n_obs
    m <- config$m
    theta <- config$theta

    # Generate component lifetimes
    Tm <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
        Tm[, j] <- rexp(n, rate = theta[j])
    }

    # System lifetime and failed component
    sys_time <- apply(Tm, 1, min)
    k <- apply(Tm, 1, which.min)

    # Right censoring
    tau <- rep(config$tau, length.out = n)
    delta <- as.numeric(sys_time <= tau)
    obs_time <- pmin(sys_time, tau)

    # Generate candidate sets using power-weighted hazards
    # p_j = base_p * (theta_j^alpha / max(theta^alpha))
    weights <- theta^alpha
    weights <- weights / max(weights)  # Normalize to [0, 1]

    C <- matrix(FALSE, nrow = n, ncol = m)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            ki <- k[i]
            for (j in seq_len(m)) {
                if (j == ki) {
                    C[i, j] <- TRUE  # C1: failed component always in
                } else {
                    # Probability proportional to hazard weight
                    C[i, j] <- runif(1) <= base_p * weights[j]
                }
            }
        }
    }

    list(
        t = obs_time,
        delta = delta,
        C = C,
        k = k,
        true_theta = theta,
        true_alpha = alpha,
        true_base_p = base_p,
        scenario = "Relaxed-C3"
    )
}

# -----------------------------------------------------------------------------
# Model Fitting Functions
# -----------------------------------------------------------------------------

#' Fit C1-C2-C3 model
#'
#' @param sim Simulated data
#' @return Fitted model results
#' @keywords internal
fit_c1_c2_c3 <- function(sim) {
    fit <- tryCatch(
        mle_exp_series(sim$t, sim$C, sim$delta),
        error = function(e) list(theta = rep(NA, ncol(sim$C)), converged = FALSE,
                                  loglik = NA_real_)
    )
    fit$model <- "C1-C2-C3"
    if (is.null(fit$loglik)) fit$loglik <- NA_real_
    fit
}

#' Fit relaxed C2 model with known P
#'
#' @param sim Simulated data
#' @param P Known P matrix (if NULL, uses true_P from sim)
#' @return Fitted model results
#' @keywords internal
fit_relaxed_c2_known_P <- function(sim, P = NULL) {
    if (is.null(P)) P <- sim$true_P
    fit <- tryCatch(
        mle_exp_series_c1_c3(sim$t, sim$C, sim$delta, fixed_P = P),
        error = function(e) list(theta = rep(NA, ncol(sim$C)), converged = FALSE,
                                  loglik = NA_real_)
    )
    fit$model <- "Relaxed-C2-known-P"
    if (is.null(fit$loglik)) fit$loglik <- NA_real_
    fit
}

#' Fit relaxed C2 model with unknown P (joint estimation)
#'
#' @param sim Simulated data
#' @param theta0 Initial theta (if NULL, uses method of moments)
#' @param P0 Initial P matrix (if NULL, uses uniform 0.5)
#' @return Fitted model results
#' @keywords internal
fit_relaxed_c2_unknown_P <- function(sim, theta0 = NULL, P0 = NULL) {
    fit <- tryCatch(
        mle_exp_series_c1_c3(sim$t, sim$C, sim$delta, theta0 = theta0, P0 = P0),
        error = function(e) list(theta = rep(NA, ncol(sim$C)), converged = FALSE,
                                  P = matrix(NA, ncol(sim$C), ncol(sim$C)),
                                  loglik = NA_real_)
    )
    fit$model <- "Relaxed-C2-unknown-P"
    if (is.null(fit$loglik)) fit$loglik <- NA_real_
    fit
}

# -----------------------------------------------------------------------------
# Single Replication Functions
# -----------------------------------------------------------------------------

#' Run a single simulation replication
#'
#' @param config Simulation configuration
#' @param scenario Which scenario to run (1-5)
#' @param ... Additional parameters for the scenario
#' @return Data frame with results from this replication
#' @keywords internal
run_single_replication <- function(config, scenario, ...) {
    args <- list(...)

    results <- list()

    if (scenario == 1) {
        # Scenario 1: C1-C2-C3 data, C1-C2-C3 model (baseline)
        p <- args$p %||% 0.3
        sim <- generate_c1_c2_c3(config, p = p)
        fit <- fit_c1_c2_c3(sim)

        results <- data.frame(
            scenario = 1,
            scenario_name = "C1-C2-C3 data, C1-C2-C3 model",
            model = fit$model,
            converged = fit$converged,
            theta_true = I(list(sim$true_theta)),
            theta_est = I(list(fit$theta)),
            loglik = fit$loglik
        )

    } else if (scenario == 2) {
        # Scenario 2: C1-C2-C3 data, Relaxed C2 model
        p <- args$p %||% 0.3
        sim <- generate_c1_c2_c3(config, p = p)

        # Create uniform P matrix matching the true model
        P_uniform <- make_P_matrix(config$m, "uniform", p = p)
        sim$true_P <- P_uniform

        fit <- fit_relaxed_c2_unknown_P(sim)

        results <- data.frame(
            scenario = 2,
            scenario_name = "C1-C2-C3 data, Relaxed-C2 model",
            model = fit$model,
            converged = fit$converged,
            theta_true = I(list(sim$true_theta)),
            theta_est = I(list(fit$theta)),
            loglik = fit$loglik
        )

    } else if (scenario == "2b") {
        # Scenario 2b: C1-C2-C3 data, Relaxed C2 model with KNOWN P
        # Tests if the model works when P is correctly specified
        p <- args$p %||% 0.3
        sim <- generate_c1_c2_c3(config, p = p)

        # Create uniform P matrix matching the true model
        P_uniform <- make_P_matrix(config$m, "uniform", p = p)
        sim$true_P <- P_uniform

        # Fit with KNOWN P (not estimating P)
        fit <- fit_relaxed_c2_known_P(sim, P = P_uniform)

        results <- data.frame(
            scenario = "2b",
            scenario_name = "C1-C2-C3 data, Relaxed-C2 (known P)",
            model = fit$model,
            converged = fit$converged,
            theta_true = I(list(sim$true_theta)),
            theta_est = I(list(fit$theta)),
            loglik = fit$loglik
        )

    } else if (scenario == "4b") {
        # Scenario 4b: Relaxed C2 data, Relaxed C2 model with KNOWN P
        # Tests if the model recovers parameters when P is correctly specified
        P <- args$P
        if (is.null(P)) {
            P <- make_P_matrix(config$m, "full",
                               values = seq(0.2, 0.6, length.out = config$m * (config$m - 1)))
        }
        sim <- generate_relaxed_c2(config, P = P)

        # Fit with KNOWN P (not estimating P)
        fit <- fit_relaxed_c2_known_P(sim, P = P)

        results <- data.frame(
            scenario = "4b",
            scenario_name = "Relaxed-C2 data, Relaxed-C2 (known P)",
            model = fit$model,
            converged = fit$converged,
            theta_true = I(list(sim$true_theta)),
            theta_est = I(list(fit$theta)),
            loglik = fit$loglik
        )

    } else if (scenario == 3) {
        # Scenario 3: Relaxed C2 data, C1-C2-C3 model (misspecified)
        P <- args$P
        if (is.null(P)) {
            # Default: asymmetric P that violates C2
            P <- make_P_matrix(config$m, "full",
                               values = seq(0.2, 0.6, length.out = config$m * (config$m - 1)))
        }
        sim <- generate_relaxed_c2(config, P = P)
        fit <- fit_c1_c2_c3(sim)

        results <- data.frame(
            scenario = 3,
            scenario_name = "Relaxed-C2 data, C1-C2-C3 model",
            model = fit$model,
            converged = fit$converged,
            theta_true = I(list(sim$true_theta)),
            theta_est = I(list(fit$theta)),
            loglik = fit$loglik
        )

    } else if (scenario == 4) {
        # Scenario 4: Relaxed C2 data, Relaxed C2 model (correct)
        P <- args$P
        if (is.null(P)) {
            P <- make_P_matrix(config$m, "full",
                               values = seq(0.2, 0.6, length.out = config$m * (config$m - 1)))
        }
        sim <- generate_relaxed_c2(config, P = P)
        fit <- fit_relaxed_c2_unknown_P(sim)

        results <- data.frame(
            scenario = 4,
            scenario_name = "Relaxed-C2 data, Relaxed-C2 model",
            model = fit$model,
            converged = fit$converged,
            theta_true = I(list(sim$true_theta)),
            theta_est = I(list(fit$theta)),
            loglik = fit$loglik
        )

    } else if (scenario == 5) {
        # Scenario 5: C1-C2-C3 data, Relaxed C3 model (misspecified)
        p <- args$p %||% 0.3
        assumed_alpha <- args$assumed_alpha %||% 1

        sim <- generate_c1_c2_c3(config, p = p)

        fit <- tryCatch(
            mle_exp_series_c1_c2(sim$t, sim$C, sim$delta,
                                  alpha = assumed_alpha, base_p = p),
            error = function(e) list(theta = rep(NA, config$m),
                                      converged = FALSE, loglik = NA_real_)
        )

        results <- data.frame(
            scenario = 5,
            scenario_name = "C1-C2-C3 data, Relaxed-C3 model",
            model = "Relaxed-C3",
            converged = fit$converged,
            theta_true = I(list(sim$true_theta)),
            theta_est = I(list(fit$theta)),
            loglik = fit$loglik
        )

    } else if (scenario == 6) {
        # Scenario 6: Relaxed C3 data, C1-C2-C3 model (misspecified)
        alpha <- args$alpha %||% 1
        base_p <- args$base_p %||% 0.5

        sim <- generate_relaxed_c3(config, alpha = alpha, base_p = base_p)
        fit <- fit_c1_c2_c3(sim)

        results <- data.frame(
            scenario = 6,
            scenario_name = "Relaxed-C3 data, C1-C2-C3 model",
            model = fit$model,
            converged = fit$converged,
            theta_true = I(list(sim$true_theta)),
            theta_est = I(list(fit$theta)),
            loglik = fit$loglik
        )

    } else if (scenario == "6b") {
        # Scenario 6b: Relaxed C3 data, Relaxed C3 model (known alpha)
        alpha <- args$alpha %||% 1
        base_p <- args$base_p %||% 0.5

        sim <- generate_relaxed_c3(config, alpha = alpha, base_p = base_p)

        fit <- tryCatch(
            mle_exp_series_c1_c2(sim$t, sim$C, sim$delta,
                                  alpha = alpha, base_p = base_p),
            error = function(e) list(theta = rep(NA, config$m),
                                      converged = FALSE, loglik = NA_real_)
        )

        results <- data.frame(
            scenario = "6b",
            scenario_name = "Relaxed-C3 data, Relaxed-C3 (known alpha)",
            model = "Relaxed-C3-known-alpha",
            converged = fit$converged,
            theta_true = I(list(sim$true_theta)),
            theta_est = I(list(fit$theta)),
            loglik = fit$loglik
        )
    }

    results
}

# -----------------------------------------------------------------------------
# Main Simulation Functions
# -----------------------------------------------------------------------------

#' Run simulation study for a single scenario
#'
#' @param config Simulation configuration
#' @param scenario Scenario number (1-5)
#' @param ... Additional parameters passed to run_single_replication
#' @param verbose Print progress
#' @return Data frame with all replication results
#' @export
run_simulation_scenario <- function(config, scenario, ..., verbose = TRUE) {
    set.seed(config$seed)

    results_list <- vector("list", config$n_sim)

    for (i in seq_len(config$n_sim)) {
        if (verbose && i %% 10 == 0) {
            cat(sprintf("Scenario %d: Replication %d/%d\n", scenario, i, config$n_sim))
        }
        results_list[[i]] <- run_single_replication(config, scenario, ...)
        results_list[[i]]$replication <- i
    }

    do.call(rbind, results_list)
}

#' Run all simulation scenarios
#'
#' @param config Simulation configuration
#' @param scenarios Which scenarios to run (default: all 1-5)
#' @param ... Additional parameters
#' @param verbose Print progress
#' @return Data frame with all results
#' @export
run_all_simulations <- function(config, scenarios = 1:4, ..., verbose = TRUE) {
    all_results <- list()

    for (s in scenarios) {
        if (verbose) cat(sprintf("\n=== Running Scenario %d ===\n", s))
        all_results[[as.character(s)]] <- run_simulation_scenario(
            config, s, ..., verbose = verbose
        )
    }

    do.call(rbind, all_results)
}

# -----------------------------------------------------------------------------
# Analysis Functions
# -----------------------------------------------------------------------------

#' Compute simulation summary statistics
#'
#' @param results Data frame from run_simulation_scenario or run_all_simulations
#' @return Data frame with bias, variance, MSE, coverage for each scenario
#' @export
summarize_simulation <- function(results) {
    # Handle the list columns
    scenarios <- unique(results$scenario)

    summary_list <- lapply(scenarios, function(s) {
        res_s <- results[results$scenario == s & results$converged, ]

        if (nrow(res_s) == 0) {
            return(data.frame(
                scenario = s,
                scenario_name = NA,
                n_converged = 0,
                convergence_rate = 0
            ))
        }

        # Extract theta estimates and true values
        theta_est <- do.call(rbind, res_s$theta_est)
        theta_true <- res_s$theta_true[[1]]  # Same for all reps
        m <- length(theta_true)

        # Compute statistics for each component
        bias <- colMeans(theta_est) - theta_true
        variance <- apply(theta_est, 2, var)
        mse <- bias^2 + variance
        rmse <- sqrt(mse)

        # Relative bias and RMSE
        rel_bias <- bias / theta_true
        rel_rmse <- rmse / theta_true

        # Coverage probability (95% CI based on asymptotic normality)
        # Check if se_est is available (from individual fits)
        coverage <- rep(NA_real_, m)
        if ("se_est" %in% names(res_s)) {
            se_est <- do.call(rbind, res_s$se_est)
            for (j in seq_len(m)) {
                lower <- theta_est[, j] - 1.96 * se_est[, j]
                upper <- theta_est[, j] + 1.96 * se_est[, j]
                coverage[j] <- mean(lower <= theta_true[j] & theta_true[j] <= upper,
                                     na.rm = TRUE)
            }
        }

        data.frame(
            scenario = s,
            scenario_name = res_s$scenario_name[1],
            n_converged = nrow(res_s),
            convergence_rate = nrow(res_s) / sum(results$scenario == s),
            component = 1:m,
            theta_true = theta_true,
            theta_mean = colMeans(theta_est),
            theta_sd = apply(theta_est, 2, sd),
            bias = bias,
            rel_bias = rel_bias,
            variance = variance,
            mse = mse,
            rmse = rmse,
            rel_rmse = rel_rmse,
            coverage_95 = coverage
        )
    })

    do.call(rbind, summary_list)
}

#' Print simulation summary
#'
#' @param summary_df Output from summarize_simulation
#' @export
print_simulation_summary <- function(summary_df) {
    scenarios <- unique(summary_df$scenario)

    for (s in scenarios) {
        cat(sprintf("\n=== Scenario %s: %s ===\n", as.character(s),
                    summary_df$scenario_name[summary_df$scenario == s][1]))

        sub <- summary_df[summary_df$scenario == s, ]
        cat(sprintf("Convergence rate: %.1f%% (%d/%d)\n",
                    sub$convergence_rate[1] * 100,
                    sub$n_converged[1],
                    round(sub$n_converged[1] / sub$convergence_rate[1])))

        cat("\nComponent-wise results:\n")
        cat(sprintf("%-10s %10s %10s %10s %10s %10s\n",
                    "Component", "True", "Mean Est", "Bias", "Rel.Bias%", "Rel.RMSE%"))
        cat(paste(rep("-", 62), collapse = ""), "\n")

        for (i in seq_len(nrow(sub))) {
            cat(sprintf("%-10d %10.3f %10.3f %10.3f %10.1f%% %10.1f%%\n",
                        sub$component[i],
                        sub$theta_true[i],
                        sub$theta_mean[i],
                        sub$bias[i],
                        sub$rel_bias[i] * 100,
                        sub$rel_rmse[i] * 100))
        }
    }
}

# -----------------------------------------------------------------------------
# Convenience wrapper for quick studies
# -----------------------------------------------------------------------------

#' Run a quick simulation study comparing scenarios 1-4
#'
#' @param n_sim Number of replications (default 50 for quick run)
#' @param n_obs Observations per replication
#' @param theta True parameters
#' @param p Inclusion probability for C1-C2-C3 scenarios
#' @param P P matrix for relaxed C2 scenarios (if NULL, uses default asymmetric)
#' @param seed Random seed
#' @return List with results and summary
#' @export
#' @examples
#' \dontrun{
#' # Quick study with defaults
#' study <- quick_simulation_study(n_sim = 50)
#' print_simulation_summary(study$summary)
#'
#' # Custom parameters
#' study <- quick_simulation_study(
#'     n_sim = 100,
#'     n_obs = 300,
#'     theta = c(1, 1.5, 2),
#'     p = 0.4
#' )
#' }
quick_simulation_study <- function(n_sim = 50, n_obs = 200, theta = c(1, 2),
                                    p = 0.3, P = NULL, seed = 42) {
    config <- sim_config(n_sim = n_sim, n_obs = n_obs, theta = theta, seed = seed)

    if (is.null(P)) {
        # Create asymmetric P matrix that violates C2
        m <- length(theta)
        P <- make_P_matrix(m, "full",
                           values = seq(0.2, 0.5, length.out = m * (m - 1)))
    }

    cat("Running simulation study...\n")
    cat(sprintf("  n_sim = %d, n_obs = %d, m = %d\n", n_sim, n_obs, length(theta)))
    cat(sprintf("  theta = (%s)\n", paste(theta, collapse = ", ")))
    cat(sprintf("  p = %.2f (for C1-C2-C3 scenarios)\n", p))

    results <- run_all_simulations(config, scenarios = 1:4, p = p, P = P)
    summary_df <- summarize_simulation(results)

    list(
        config = config,
        results = results,
        summary = summary_df,
        P = P,
        p = p
    )
}

# -----------------------------------------------------------------------------
# Helper: null-coalescing operator
# -----------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
