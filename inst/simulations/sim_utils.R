# =============================================================================
# Simulation Utilities for Relaxed Candidate Set Models
# =============================================================================
#
# Core infrastructure for Monte Carlo simulation studies examining:
# - KL-divergence effects on MLE efficiency
# - Misspecification bias under violated C2 conditions
# - Identifiability analysis via FIM eigenvalues
#
# Author: Alexander Towell
# Package: md_series_system_relaxed_candidate_set_models
# =============================================================================

# -----------------------------------------------------------------------------
# Dependency Stubs (when md.tools not available)
# -----------------------------------------------------------------------------
# These provide minimal implementations when the full md.tools package
# is not installed, allowing simulations to run standalone.

#' Encode a matrix as prefixed columns in a data frame
#' @param mat Numeric matrix
#' @param prefix Character prefix for column names
#' @return Data frame with columns prefix1, prefix2, ..., prefixm
#' @keywords internal
md_encode_matrix_stub <- function(mat, prefix) {
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    df <- as.data.frame(mat)
    colnames(df) <- paste0(prefix, seq_len(ncol(mat)))
    df
}

#' Decode prefixed columns from a data frame to matrix
#' @param df Data frame
#' @param prefix Character prefix to match
#' @return Numeric matrix or NULL if no matching columns
#' @keywords internal
md_decode_matrix_stub <- function(df, prefix) {
    cols <- grep(paste0("^", prefix, "[0-9]+$"), colnames(df), value = TRUE)
    if (length(cols) == 0) return(NULL)
    # Sort by numeric suffix
    nums <- as.integer(sub(prefix, "", cols))
    cols <- cols[order(nums)]
    as.matrix(df[, cols, drop = FALSE])
}

#' Mark columns as latent (stub - returns df unchanged)
#' @param df Data frame
#' @param cols Column names to mark
#' @return Data frame (unchanged in stub)
#' @keywords internal
md_mark_latent_stub <- function(df, cols) {
    df
}

# Use stubs if md.tools not available
if (!requireNamespace("md.tools", quietly = TRUE)) {
    md_encode_matrix <- md_encode_matrix_stub
    md_decode_matrix <- md_decode_matrix_stub
    md_mark_latent <- md_mark_latent_stub
} else {
    md_encode_matrix <- md.tools::md_encode_matrix
    md_decode_matrix <- md.tools::md_decode_matrix
    md_mark_latent <- md.tools::md_mark_latent
}

# -----------------------------------------------------------------------------
# Core Data Generation Functions
# -----------------------------------------------------------------------------

#' Generate component failure times for exponential series system
#'
#' @param n Sample size (number of systems)
#' @param rates Vector of exponential rate parameters (lambda_1, ..., lambda_m)
#' @return Matrix of component failure times (n x m)
#' @export
generate_component_times <- function(n, rates) {
    m <- length(rates)
    stopifnot(n > 0, m > 0, all(rates > 0))

    Tm <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
        Tm[, j] <- stats::rexp(n, rate = rates[j])
    }
    Tm
}

#' Generate masked data for exponential series system
#'
#' This function generates complete masked data for simulation studies,
#' including component failure times, system lifetimes, censoring,
#' candidate set probabilities, and sampled candidate sets.
#'
#' @param n Sample size (number of systems)
#' @param rates Vector of exponential rate parameters (lambda_1, ..., lambda_m)
#' @param tau Right-censoring time (scalar or vector of length n)
#' @param masking_model Character: "bernoulli" (C1,C2,C3), "kl_divergence", or "informative"
#' @param masking_params List of masking model parameters:
#'   - For "bernoulli": p (probability non-failed in candidate set)
#'   - For "kl_divergence": p, d (target KL-divergence)
#'   - For "informative": alpha, beta (rank-based masking parameters)
#' @param seed Optional random seed for reproducibility
#'
#' @return A tibble with columns:
#'   - t1, ..., tm: Component failure times (latent)
#'   - t: System lifetime (observed, possibly censored)
#'   - k: Failed component index (latent)
#'   - delta: Censoring indicator (TRUE = right-censored)
#'   - q1, ..., qm: Candidate set inclusion probabilities
#'   - x1, ..., xm: Candidate set indicators (sampled)
#'
#' @export
generate_masked_data <- function(n, rates, tau,
                                  masking_model = "bernoulli",
                                  masking_params = list(p = 0.3),
                                  seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    m <- length(rates)
    stopifnot(n > 0, m > 0, all(rates > 0))

    # Generate component failure times
    Tm <- generate_component_times(n, rates)

    # System lifetime and failed component
    sys_time <- apply(Tm, 1, min)
    failed_comp <- apply(Tm, 1, which.min)

    # Apply right-censoring
    tau_vec <- rep(tau, length.out = n)
    delta <- sys_time > tau_vec
    obs_time <- pmin(sys_time, tau_vec)

    # Create base data frame
    md <- tibble::tibble(
        t = obs_time,
        k = failed_comp,
        delta = delta
    )
    md <- dplyr::bind_cols(md, md_encode_matrix(Tm, "t"))

    # Generate masking probabilities based on model
    Q <- matrix(0, nrow = n, ncol = m)

    if (masking_model == "bernoulli") {
        # Standard Bernoulli model satisfying C1, C2, C3
        p <- masking_params$p
        stopifnot(!is.null(p), p >= 0, p <= 1)

        Q <- matrix(p, nrow = n, ncol = m)
        Q[cbind(seq_len(n), failed_comp)] <- 1
        # Set to 0 for censored observations
        Q[delta, ] <- 0

    } else if (masking_model == "kl_divergence") {
        # KL-divergence constrained model
        p <- masking_params$p
        d <- masking_params$d
        alpha0 <- masking_params$alpha0 %||% 5
        beta0 <- masking_params$beta0 %||% p

        stopifnot(!is.null(p), !is.null(d), p >= 0, p <= 1, d >= 0)

        for (i in seq_len(n)) {
            if (delta[i]) {
                # Censored: no candidate set
                Q[i, ] <- 0
            } else {
                res <- generate_relaxed_cand_C1_internal(
                    d = d, ts = Tm[i, ], p = p,
                    alpha0 = alpha0, beta0 = beta0
                )
                Q[i, ] <- res$Q
            }
        }

    } else if (masking_model == "informative") {
        # Rank-based informative masking
        alpha <- masking_params$alpha
        beta <- masking_params$beta

        stopifnot(!is.null(alpha), !is.null(beta),
                  alpha >= 0, beta >= 0, beta <= 1)

        for (i in seq_len(n)) {
            if (delta[i]) {
                Q[i, ] <- 0
            } else {
                Q[i, ] <- informative_masking_by_rank_internal(
                    Tm[i, ], alpha, beta
                )
            }
        }
    } else {
        stop("Unknown masking_model: ", masking_model)
    }

    md <- dplyr::bind_cols(md, md_encode_matrix(Q, "q"))

    # Sample candidate sets from probabilities
    X <- matrix(FALSE, nrow = n, ncol = m)
    for (i in seq_len(n)) {
        X[i, ] <- stats::runif(m) <= Q[i, ]
    }
    md <- dplyr::bind_cols(md, md_encode_matrix(X, "x"))

    # Mark latent columns
    md <- md_mark_latent(md, c(paste0("t", seq_len(m)), "k"))

    md
}

#' Internal implementation of informative_masking_by_rank
#' @keywords internal
informative_masking_by_rank_internal <- function(ts, alpha, beta) {
    m <- length(ts)
    stopifnot(m > 0, alpha >= 0, beta >= 0, beta <= 1)

    ranks <- 0:(m - 2)
    wts <- beta * exp(-alpha * ranks)
    probs <- c(1, wts)

    probs[order(ts)] <- probs
    probs
}

#' Internal implementation of generate_relaxed_cand_C1
#' @keywords internal
generate_relaxed_cand_C1_internal <- function(d, ts, p,
                                               alpha0 = 1, beta0 = p,
                                               eps = 1e-6, max_iter = 5000,
                                               lr = 0.5, lambda = 1) {
    m <- length(ts)
    k <- which.min(ts)
    P <- rep(p, m)
    P[k] <- 1

    if (d == 0) {
        return(list(P = P, Q = P, KL = 0, alpha = Inf, beta = p,
                    converged = TRUE, k = k))
    }

    P_sum <- sum(P)

    # KL divergence for Bernoulli vectors
    kl_bern <- function(p_vec, q_vec) {
        kl_elem <- function(x, y) {
            if (x == 0 && y == 0) return(0)
            if (x == 1 && y == 1) return(0)
            if (x == 0) return((1 - x) * log2((1 - x) / (1 - y)))
            if (x == 1) return(x * log2(x / y))
            x * log2(x / y) + (1 - x) * log2((1 - x) / (1 - y))
        }
        sum(mapply(kl_elem, p_vec, q_vec))
    }

    # Cost function to minimize
    cost <- function(param) {
        alpha <- param[1]
        beta <- param[2]
        Q <- informative_masking_by_rank_internal(ts, alpha, beta)
        delta_sum <- sum(Q) - P_sum
        KL <- kl_bern(P, Q)
        (KL - d)^2 + lambda * delta_sum^2
    }

    # Simple gradient descent with numerical gradient
    x <- c(alpha0, beta0)
    h <- 1e-6

    for (iter in seq_len(max_iter)) {
        # Numerical gradient
        g <- numeric(2)
        fx <- cost(x)
        for (j in 1:2) {
            x_plus <- x
            x_plus[j] <- x_plus[j] + h
            g[j] <- (cost(x_plus) - fx) / h
        }

        grad_norm <- sqrt(sum(g^2))
        if (grad_norm < eps) break

        # Update with projection onto support
        x_new <- x - lr * g
        x_new[1] <- max(0, x_new[1])
        x_new[2] <- max(0, min(1, x_new[2]))
        x <- x_new
    }

    Q <- informative_masking_by_rank_internal(ts, x[1], x[2])
    KL <- kl_bern(P, Q)

    list(P = P, Q = Q, KL = KL, alpha = x[1], beta = x[2],
         converged = (grad_norm < eps), k = k)
}

# -----------------------------------------------------------------------------
# MLE Computation Functions
# -----------------------------------------------------------------------------

#' Compute MLE for exponential series system with masked data (C1,C2,C3)
#'
#' Internal implementation that does not depend on algebraic.mle
#'
#' @param md Masked data frame
#' @param theta0 Initial parameter values (optional)
#' @param sysvar Column name for system lifetime
#' @param setvar Column prefix for candidate sets
#' @param deltavar Column name for censoring indicator
#' @param lower Lower bound for parameters
#' @param upper Upper bound for parameters
#'
#' @return List with components:
#'   - theta_hat: MLE of rate parameters
#'   - loglike: Log-likelihood at MLE
#'   - fim: Observed Fisher information matrix
#'   - vcov: Variance-covariance matrix
#'   - converged: Convergence indicator
#'   - nobs: Number of observations
#'
#' @export
compute_mle_exp_series <- function(md,
                                    theta0 = NULL,
                                    sysvar = "t",
                                    setvar = "x",
                                    deltavar = "delta",
                                    lower = 1e-10,
                                    upper = Inf) {
    n <- nrow(md)
    if (n == 0) stop("md is empty")

    # Extract data
    s <- md[[sysvar]]
    C <- md_decode_matrix(md, setvar)
    if (is.null(C)) stop("No candidate set columns found")
    m <- ncol(C)

    # Handle censoring
    has_delta <- !is.null(deltavar) && deltavar %in% colnames(md)
    delta <- if (has_delta) md[[deltavar]] else rep(FALSE, n)

    # Log-likelihood function
    ll_fn <- function(theta) {
        if (any(theta <= 0)) return(-Inf)

        sum_theta <- sum(theta)
        ll <- -sum(s) * sum_theta

        for (i in seq_len(n)) {
            if (!delta[i]) {
                hazard_sum <- sum(theta[C[i, ]])
                if (hazard_sum <= 0) return(-Inf)
                ll <- ll + log(hazard_sum)
            }
        }
        ll
    }

    # Score function
    score_fn <- function(theta) {
        if (any(theta <= 0)) return(rep(NA_real_, m))

        v <- rep(-sum(s), m)
        for (i in seq_len(n)) {
            if (!delta[i]) {
                C_i <- C[i, ]
                hazard_sum <- sum(theta[C_i])
                if (hazard_sum <= 0) return(rep(NA_real_, m))
                for (j in seq_len(m)) {
                    if (C_i[j]) v[j] <- v[j] + 1 / hazard_sum
                }
            }
        }
        v
    }

    # FIM function
    fim_fn <- function(theta) {
        if (any(theta <= 0)) return(matrix(NA_real_, m, m))

        I <- matrix(0, nrow = m, ncol = m)
        for (i in seq_len(n)) {
            if (!delta[i]) {
                C_i <- C[i, ]
                hazard_sum <- sum(theta[C_i])
                if (hazard_sum <= 0) return(matrix(NA_real_, m, m))
                inv_sq <- 1 / hazard_sum^2
                for (j in seq_len(m)) {
                    if (C_i[j]) {
                        for (k in seq_len(m)) {
                            if (C_i[k]) I[j, k] <- I[j, k] + inv_sq
                        }
                    }
                }
            }
        }
        I
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

    # Compute FIM and variance-covariance
    fim <- fim_fn(result$par)
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
        score = score_fn(result$par)
    )
}

# -----------------------------------------------------------------------------
# Statistical Analysis Functions
# -----------------------------------------------------------------------------

#' Compute MLE statistics (bias, variance, MSE, coverage)
#'
#' @param estimates Matrix of parameter estimates (B x m)
#' @param theta_true True parameter values (vector of length m)
#' @param vcov_list List of variance-covariance matrices (optional)
#' @param alpha Significance level for confidence intervals (default 0.05)
#'
#' @return List with:
#'   - bias: Bias for each parameter
#'   - variance: Variance for each parameter
#'   - mse: Mean squared error for each parameter
#'   - rmse: Root mean squared error
#'   - coverage: Coverage probability (if vcov_list provided)
#'   - mean_ci_width: Average CI width (if vcov_list provided)
#'   - n_valid: Number of valid (converged) estimates
#'
#' @export
compute_mle_stats <- function(estimates, theta_true, vcov_list = NULL, alpha = 0.05) {
    if (!is.matrix(estimates)) estimates <- matrix(estimates, ncol = 1)

    # Remove rows with NA or Inf
    valid_rows <- apply(estimates, 1, function(x) all(is.finite(x)))
    estimates <- estimates[valid_rows, , drop = FALSE]
    n_valid <- nrow(estimates)

    if (n_valid == 0) {
        m <- length(theta_true)
        return(list(
            bias = rep(NA, m),
            variance = rep(NA, m),
            mse = rep(NA, m),
            rmse = rep(NA, m),
            coverage = rep(NA, m),
            mean_ci_width = rep(NA, m),
            n_valid = 0
        ))
    }

    m <- ncol(estimates)
    stopifnot(length(theta_true) == m)

    # Compute statistics
    mean_est <- colMeans(estimates)
    bias <- mean_est - theta_true
    variance <- apply(estimates, 2, stats::var)
    mse <- bias^2 + variance
    rmse <- sqrt(mse)

    # Compute coverage if vcov provided
    coverage <- NULL
    mean_ci_width <- NULL

    if (!is.null(vcov_list)) {
        vcov_list <- vcov_list[valid_rows]
        z_crit <- stats::qnorm(1 - alpha / 2)

        in_ci <- matrix(FALSE, nrow = n_valid, ncol = m)
        ci_widths <- matrix(NA, nrow = n_valid, ncol = m)

        for (i in seq_len(n_valid)) {
            vcov <- vcov_list[[i]]
            if (!is.null(vcov) && all(is.finite(vcov))) {
                se <- sqrt(diag(vcov))
                lower <- estimates[i, ] - z_crit * se
                upper <- estimates[i, ] + z_crit * se
                in_ci[i, ] <- (theta_true >= lower) & (theta_true <= upper)
                ci_widths[i, ] <- 2 * z_crit * se
            }
        }

        coverage <- colMeans(in_ci, na.rm = TRUE)
        mean_ci_width <- colMeans(ci_widths, na.rm = TRUE)
    }

    list(
        bias = bias,
        variance = variance,
        mse = mse,
        rmse = rmse,
        coverage = coverage,
        mean_ci_width = mean_ci_width,
        n_valid = n_valid
    )
}

#' Compute relative efficiency
#'
#' @param mse_test MSE of test estimator
#' @param mse_baseline MSE of baseline estimator
#'
#' @return Relative efficiency (mse_baseline / mse_test)
#' @export
relative_efficiency <- function(mse_test, mse_baseline) {
    ifelse(mse_test > 0, mse_baseline / mse_test, NA)
}

# -----------------------------------------------------------------------------
# Simulation Runner Functions
# -----------------------------------------------------------------------------

#' Run Monte Carlo simulation for a single scenario
#'
#' @param n Sample size
#' @param rates True rate parameters
#' @param tau Right-censoring time
#' @param masking_model Masking model type
#' @param masking_params Masking model parameters
#' @param B Number of replications
#' @param seed Base random seed
#' @param progress Logical, show progress bar
#'
#' @return List with:
#'   - estimates: Matrix of MLE estimates (B x m)
#'   - vcov_list: List of variance-covariance matrices
#'   - converged: Vector of convergence indicators
#'   - stats: Summary statistics from compute_mle_stats
#'
#' @export
run_scenario <- function(n, rates, tau,
                          masking_model = "bernoulli",
                          masking_params = list(p = 0.3),
                          B = 500,
                          seed = NULL,
                          progress = TRUE) {
    m <- length(rates)

    # Storage
    estimates <- matrix(NA_real_, nrow = B, ncol = m)
    vcov_list <- vector("list", B)
    converged <- logical(B)

    # Set base seed
    if (!is.null(seed)) set.seed(seed)

    # Run replications
    for (b in seq_len(B)) {
        if (progress && b %% 50 == 0) {
            cat(sprintf("  Replication %d/%d\n", b, B))
        }

        tryCatch({
            # Generate data
            md <- generate_masked_data(
                n = n,
                rates = rates,
                tau = tau,
                masking_model = masking_model,
                masking_params = masking_params,
                seed = NULL  # Use current RNG state
            )

            # Compute MLE
            mle_result <- compute_mle_exp_series(md, theta0 = rates)

            estimates[b, ] <- mle_result$theta_hat
            vcov_list[[b]] <- mle_result$vcov
            converged[b] <- mle_result$converged

        }, error = function(e) {
            # Leave as NA on error
        })
    }

    # Compute statistics
    stats <- compute_mle_stats(
        estimates = estimates,
        theta_true = rates,
        vcov_list = vcov_list
    )

    list(
        estimates = estimates,
        vcov_list = vcov_list,
        converged = converged,
        stats = stats,
        params = list(
            n = n,
            rates = rates,
            tau = tau,
            masking_model = masking_model,
            masking_params = masking_params,
            B = B
        )
    )
}

#' Run simulation study across multiple scenarios (parallel support)
#'
#' @param scenarios List of scenario specifications, each containing:
#'   - n, rates, tau, masking_model, masking_params
#' @param B Number of replications per scenario
#' @param n_cores Number of cores for parallel execution (1 = sequential)
#' @param seed Base random seed
#' @param progress Show progress
#'
#' @return List of results, one per scenario
#' @export
run_simulation <- function(scenarios,
                            B = 500,
                            n_cores = 1,
                            seed = NULL,
                            progress = TRUE) {
    n_scenarios <- length(scenarios)

    if (progress) {
        cat(sprintf("Running %d scenarios with B=%d replications each\n",
                    n_scenarios, B))
    }

    # Set seeds for reproducibility
    if (!is.null(seed)) {
        set.seed(seed)
        seeds <- sample.int(.Machine$integer.max, n_scenarios)
    } else {
        seeds <- rep(NA, n_scenarios)
    }

    # Run scenarios
    if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
        # Parallel execution
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl))

        # Export necessary functions
        parallel::clusterExport(cl, c(
            "run_scenario", "generate_masked_data", "compute_mle_exp_series",
            "compute_mle_stats", "generate_component_times",
            "informative_masking_by_rank_internal",
            "generate_relaxed_cand_C1_internal",
            "md_encode_matrix", "md_decode_matrix", "md_mark_latent"
        ), envir = environment())

        results <- parallel::parLapply(cl, seq_len(n_scenarios), function(i) {
            sc <- scenarios[[i]]
            if (progress) cat(sprintf("Scenario %d/%d\n", i, n_scenarios))

            run_scenario(
                n = sc$n,
                rates = sc$rates,
                tau = sc$tau,
                masking_model = sc$masking_model %||% "bernoulli",
                masking_params = sc$masking_params %||% list(p = 0.3),
                B = B,
                seed = seeds[i],
                progress = FALSE
            )
        })
    } else {
        # Sequential execution
        results <- lapply(seq_len(n_scenarios), function(i) {
            sc <- scenarios[[i]]
            if (progress) {
                cat(sprintf("\nScenario %d/%d: n=%d, model=%s\n",
                            i, n_scenarios, sc$n,
                            sc$masking_model %||% "bernoulli"))
            }

            run_scenario(
                n = sc$n,
                rates = sc$rates,
                tau = sc$tau,
                masking_model = sc$masking_model %||% "bernoulli",
                masking_params = sc$masking_params %||% list(p = 0.3),
                B = B,
                seed = seeds[i],
                progress = progress
            )
        })
    }

    results
}

# -----------------------------------------------------------------------------
# I/O Functions
# -----------------------------------------------------------------------------

#' Save simulation results to file
#'
#' @param results Simulation results object
#' @param file Path to save file (.rds or .csv)
#' @param format "rds" for full R object, "csv" for summary table
#'
#' @export
save_results <- function(results, file, format = "rds") {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

    if (format == "rds") {
        saveRDS(results, file)
    } else if (format == "csv") {
        # Convert to summary table
        summary_df <- results_to_dataframe(results)
        utils::write.csv(summary_df, file, row.names = FALSE)
    } else {
        stop("Unknown format: ", format)
    }

    invisible(file)
}

#' Load simulation results from file
#'
#' @param file Path to results file
#' @return Simulation results object
#' @export
load_results <- function(file) {
    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    ext <- tools::file_ext(file)
    if (ext == "rds") {
        readRDS(file)
    } else if (ext == "csv") {
        utils::read.csv(file, stringsAsFactors = FALSE)
    } else {
        stop("Unknown file extension: ", ext)
    }
}

#' Convert simulation results to summary data frame
#'
#' @param results List of scenario results from run_simulation
#' @return Data frame with summary statistics
#' @export
results_to_dataframe <- function(results) {
    rows <- lapply(seq_along(results), function(i) {
        r <- results[[i]]
        params <- r$params
        stats <- r$stats
        m <- length(params$rates)

        # Build row for each component
        rows_i <- lapply(seq_len(m), function(j) {
            data.frame(
                scenario = i,
                n = params$n,
                masking_model = params$masking_model,
                component = j,
                theta_true = params$rates[j],
                bias = stats$bias[j],
                variance = stats$variance[j],
                mse = stats$mse[j],
                rmse = stats$rmse[j],
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

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Compute right-censoring time for target censoring proportion
#'
#' @param rates Exponential rate parameters
#' @param target_cens_prop Target censoring proportion (0 to 1)
#' @param method "quantile" uses system lifetime quantile
#'
#' @return Right-censoring time tau
#' @export
compute_tau_for_censoring <- function(rates, target_cens_prop, method = "quantile") {
    # For exponential series, system hazard is sum(rates)
    # System lifetime ~ Exp(sum(rates))
    # P(T > tau) = exp(-sum(rates) * tau) = target_cens_prop
    # tau = -log(target_cens_prop) / sum(rates)

    if (target_cens_prop <= 0) return(Inf)
    if (target_cens_prop >= 1) return(0)

    -log(target_cens_prop) / sum(rates)
}

#' Print simulation summary
#'
#' @param results Results from run_simulation or run_scenario
#' @export
print_simulation_summary <- function(results) {
    if ("stats" %in% names(results)) {
        # Single scenario result
        results <- list(results)
    }

    cat("=== Simulation Study Summary ===\n\n")

    for (i in seq_along(results)) {
        r <- results[[i]]
        params <- r$params
        stats <- r$stats

        cat(sprintf("Scenario %d:\n", i))
        cat(sprintf("  Sample size: %d\n", params$n))
        cat(sprintf("  Masking model: %s\n", params$masking_model))
        cat(sprintf("  Valid replications: %d/%d\n", stats$n_valid, params$B))
        cat("\n")

        cat("  Component | True   | Bias     | RMSE     | Coverage\n")
        cat("  " , rep("-", 55), "\n", sep = "")

        m <- length(params$rates)
        for (j in seq_len(m)) {
            cov_str <- if (!is.null(stats$coverage))
                sprintf("%.3f", stats$coverage[j]) else "N/A"
            cat(sprintf("  %9d | %.4f | %8.5f | %8.5f | %s\n",
                        j, params$rates[j], stats$bias[j],
                        stats$rmse[j], cov_str))
        }
        cat("\n")
    }
}

# -----------------------------------------------------------------------------
# FIM Analysis Functions (for identifiability studies)
# -----------------------------------------------------------------------------

#' Compute expected FIM eigenvalue summary
#'
#' @param fim Fisher information matrix
#' @return List with eigenvalues, condition number, smallest eigenvalue
#' @export
fim_eigenvalue_analysis <- function(fim) {
    if (any(is.na(fim)) || any(!is.finite(fim))) {
        return(list(
            eigenvalues = rep(NA, nrow(fim)),
            condition_number = NA,
            smallest = NA,
            near_singular = TRUE
        ))
    }

    eig <- eigen(fim, symmetric = TRUE, only.values = TRUE)$values

    list(
        eigenvalues = eig,
        condition_number = max(eig) / max(min(eig), .Machine$double.eps),
        smallest = min(eig),
        near_singular = min(eig) < 1e-10
    )
}

#' Compute Monte Carlo expected FIM
#'
#' Average FIM over many random datasets to approximate expected FIM
#'
#' @param n Sample size
#' @param rates True parameters
#' @param tau Censoring time
#' @param masking_params Masking parameters
#' @param n_mc Number of Monte Carlo samples
#' @param seed Random seed
#'
#' @return Average FIM matrix
#' @export
compute_expected_fim_mc <- function(n, rates, tau,
                                     masking_params = list(p = 0.3),
                                     n_mc = 100,
                                     seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    m <- length(rates)
    fim_sum <- matrix(0, m, m)
    n_valid <- 0

    for (i in seq_len(n_mc)) {
        md <- generate_masked_data(
            n = n, rates = rates, tau = tau,
            masking_model = "bernoulli",
            masking_params = masking_params
        )

        mle_result <- compute_mle_exp_series(md, theta0 = rates)

        if (all(is.finite(mle_result$fim))) {
            fim_sum <- fim_sum + mle_result$fim
            n_valid <- n_valid + 1
        }
    }

    if (n_valid > 0) {
        fim_sum / n_valid
    } else {
        matrix(NA, m, m)
    }
}
