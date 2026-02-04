# =============================================================================
# Weibull Series System with Masked Data (C1, C3) - Relaxed C2
# Informative masking where P[j,k] = P(j in C | K = k) can vary with k
# =============================================================================
#
# Parameter convention: theta = c(k1, lambda1, k2, lambda2, ..., km, lambdam)
# where k_j = shape and lambda_j = scale for component j.
#
# The P matrix encodes masking probabilities: P[j,k] = P(j in C | K = k)
# - Diagonal must be 1 (C1 condition)
# - When P is uniform (all off-diagonal equal), this reduces to C1-C2-C3
# - Non-uniform P models informative masking

# -----------------------------------------------------------------------------
# Log-Likelihood
# -----------------------------------------------------------------------------

#' Log-likelihood for Weibull series system with masked data (C1, C3) - Relaxed C2
#'
#' Returns a log-likelihood function for masked data from a Weibull series
#' system under conditions C1 and C3, with relaxed C2 (informative masking).
#'
#' The log-likelihood is:
#' `l(theta) = sum_i [-sum_j (t_i/lambda_j)^k_j] + sum_{i:delta_i=1} log(sum_{j in C_i} h_j(t_i) * pi_j(c_i))`
#'
#' where h_j(t) = (k_j/lambda_j)(t/lambda_j)^(k_j-1) is the Weibull hazard, and
#' pi_j(c_i) = P(C = c_i | K = j) computed from the general Bernoulli model.
#'
#' @param t Numeric vector of system lifetimes
#' @param C Logical matrix of candidate sets (n x m)
#' @param delta Censoring indicators (1=observed, 0=censored)
#' @param P Inclusion probability matrix (m x m): `P[j,k] = P(j in C | K = k)`.
#'        Diagonal must be 1 (C1 condition).
#' @return Function that takes theta = c(k1,lambda1,k2,lambda2,...) and returns log-likelihood
#' @export
loglik_wei_series_c1_c3 <- function(t, C, delta = NULL, P) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    if (nrow(P) != m || ncol(P) != m) stop("P must be m x m matrix")
    if (!all(abs(diag(P) - 1) < 1e-10)) stop("Diagonal of P must be 1 (C1 condition)")

    # Precompute pi values for each uncensored observation
    pi_list <- vector("list", n)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            pi_list[[i]] <- compute_pi_all(C[i, ], P)
        }
    }

    function(theta) {
        params <- unpack_wei_params(theta)
        if (params$m != m) stop("theta implies ", params$m, " components but C has ", m)
        shapes <- params$shapes
        scales <- params$scales
        if (any(shapes <= 0) || any(scales <= 0)) return(-Inf)

        # Survival contribution: -sum_i sum_j (t_i / lambda_j)^k_j
        ll <- 0
        for (j in seq_len(m)) {
            ll <- ll - sum((t / scales[j])^shapes[j])
        }

        # Weighted hazard contribution for uncensored observations
        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                candidates <- which(C[i, ])
                pi_vals <- pi_list[[i]]

                # h_j(t_i) = (k_j / lambda_j) * (t_i / lambda_j)^(k_j - 1)
                hazard_sum <- 0
                for (idx in seq_along(candidates)) {
                    j <- candidates[idx]
                    h_j <- (shapes[j] / scales[j]) * (t[i] / scales[j])^(shapes[j] - 1)
                    hazard_sum <- hazard_sum + h_j * pi_vals[idx]
                }
                if (hazard_sum <= 0) return(-Inf)
                ll <- ll + log(hazard_sum)
            }
        }
        ll
    }
}

#' Score (gradient) for Weibull series system (C1, C3) - Relaxed C2
#'
#' Returns the score function (gradient of log-likelihood) for masked data
#' with informative masking (relaxed C2).
#'
#' @inheritParams loglik_wei_series_c1_c3
#' @return Function that takes theta and returns gradient vector
#' @export
score_wei_series_c1_c3 <- function(t, C, delta = NULL, P) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    if (nrow(P) != m || ncol(P) != m) stop("P must be m x m matrix")

    # Precompute pi values
    pi_list <- vector("list", n)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            pi_list[[i]] <- compute_pi_all(C[i, ], P)
        }
    }

    function(theta) {
        params <- unpack_wei_params(theta)
        shapes <- params$shapes
        scales <- params$scales
        if (any(shapes <= 0) || any(scales <= 0)) return(rep(NA_real_, 2*m))

        g <- numeric(2 * m)

        # Survival derivatives
        for (j in seq_len(m)) {
            z <- t / scales[j]
            z_k <- z^shapes[j]
            log_z <- log(z)
            log_z[t == 0] <- 0

            # d/d(k_j) of -sum (t/lambda_j)^k_j = -sum z^k * log(z)
            g[2*j - 1] <- -sum(z_k * log_z)

            # d/d(lambda_j) of -sum (t/lambda_j)^k_j = sum k_j/lambda_j * z^k_j
            g[2*j] <- sum(shapes[j] / scales[j] * z_k)
        }

        # Weighted hazard derivatives for uncensored observations
        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                candidates <- which(C[i, ])
                pi_vals <- pi_list[[i]]

                # Compute weighted hazard for each candidate
                h_vals <- numeric(m)
                weighted_h <- numeric(m)
                for (idx in seq_along(candidates)) {
                    j <- candidates[idx]
                    z <- t[i] / scales[j]
                    h_vals[j] <- (shapes[j] / scales[j]) * z^(shapes[j] - 1)
                    weighted_h[j] <- h_vals[j] * pi_vals[idx]
                }
                hazard_sum <- sum(weighted_h[candidates])
                if (hazard_sum <= 0) return(rep(NA_real_, 2*m))

                # Derivatives of log(weighted_hazard_sum) w.r.t. theta
                for (idx in seq_along(candidates)) {
                    j <- candidates[idx]
                    z <- t[i] / scales[j]
                    log_z <- log(z)
                    if (t[i] == 0) log_z <- 0

                    # d(h_j)/d(k_j) = h_j * (1/k_j + log(z))
                    dh_dk <- h_vals[j] * (1 / shapes[j] + log_z)

                    # d(h_j)/d(lambda_j) = -h_j * k_j / lambda_j
                    dh_dlam <- -h_vals[j] * shapes[j] / scales[j]

                    # Contribution to score: pi_j * dh_j/d(theta) / weighted_sum
                    g[2*j - 1] <- g[2*j - 1] + pi_vals[idx] * dh_dk / hazard_sum
                    g[2*j] <- g[2*j] + pi_vals[idx] * dh_dlam / hazard_sum
                }
            }
        }
        g
    }
}

#' Fisher Information Matrix for Weibull series system (C1, C3) - Relaxed C2
#'
#' Computed numerically as the negative Hessian of the log-likelihood.
#'
#' @inheritParams loglik_wei_series_c1_c3
#' @return Function that takes theta and returns (2m x 2m) FIM
#' @export
fim_wei_series_c1_c3 <- function(t, C, delta = NULL, P) {
    ll_fn <- loglik_wei_series_c1_c3(t, C, delta, P)

    function(theta) {
        p <- length(theta)
        h <- 1e-5
        H <- matrix(0, p, p)
        for (j in seq_len(p)) {
            for (k in seq_len(p)) {
                theta_pp <- theta_pm <- theta_mp <- theta_mm <- theta
                theta_pp[j] <- theta_pp[j] + h
                theta_pp[k] <- theta_pp[k] + h
                theta_pm[j] <- theta_pm[j] + h
                theta_pm[k] <- theta_pm[k] - h
                theta_mp[j] <- theta_mp[j] - h
                theta_mp[k] <- theta_mp[k] + h
                theta_mm[j] <- theta_mm[j] - h
                theta_mm[k] <- theta_mm[k] - h
                H[j, k] <- (ll_fn(theta_pp) - ll_fn(theta_pm) -
                           ll_fn(theta_mp) + ll_fn(theta_mm)) / (4 * h^2)
            }
        }
        -H
    }
}

# -----------------------------------------------------------------------------
# Maximum Likelihood Estimation
# -----------------------------------------------------------------------------

#' MLE for Weibull series system (C1, C3) - Relaxed C2
#'
#' Computes MLE of Weibull parameters under informative masking with known P.
#'
#' @inheritParams loglik_wei_series_c1_c3
#' @param theta0 Initial parameters c(k1,lambda1,...). If NULL, uses heuristic.
#' @return List with theta, shapes, scales, se, vcov, loglik, converged, fim, m
#' @export
mle_wei_series_c1_c3 <- function(t, C, delta = NULL, P, theta0 = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)

    ll_fn <- loglik_wei_series_c1_c3(t, C, delta, P)
    sc_fn <- score_wei_series_c1_c3(t, C, delta, P)

    # Initialize: use exponential MLE to get rate, then convert
    if (is.null(theta0)) {
        n_unc <- sum(delta == 1)
        total_rate <- n_unc / sum(t)
        rate_each <- total_rate / m
        theta0 <- rep(c(1, 1 / rate_each), m)
    }

    # Lower bounds: all parameters > 0
    lower <- rep(1e-6, 2 * m)

    result <- optim(
        par = theta0,
        fn = ll_fn,
        gr = sc_fn,
        method = "L-BFGS-B",
        lower = lower,
        control = list(fnscale = -1, maxit = 1000)
    )

    # Compute FIM
    fim <- fim_wei_series_c1_c3(t, C, delta, P)(result$par)
    vcov <- tryCatch(solve(fim), error = function(e) matrix(NA, 2*m, 2*m))
    se <- sqrt(diag(vcov))

    params <- unpack_wei_params(result$par)

    list(
        theta = result$par,
        shapes = params$shapes,
        scales = params$scales,
        se = se,
        vcov = vcov,
        loglik = result$value,
        converged = result$convergence == 0,
        fim = fim,
        n = n,
        m = m
    )
}

# -----------------------------------------------------------------------------
# Data Generation
# -----------------------------------------------------------------------------

#' Generate masked data from Weibull series with general Bernoulli model (C1, C3)
#'
#' Generates masked data where the candidate set inclusion probability depends
#' on which component failed (relaxed C2).
#'
#' @param n Sample size
#' @param shapes Shape parameters (k_1, ..., k_m)
#' @param scales Scale parameters (lambda_1, ..., lambda_m)
#' @param P Inclusion probability matrix (m x m): `P[j,k] = P(j in C | K = k)`
#' @param tau Right-censoring time
#' @param keep_latent If TRUE, include latent component times in output
#' @return List with t, delta, C, k, m, shapes, scales, and optionally Tm, P
#' @export
rwei_series_md_c1_c3 <- function(n, shapes, scales, P, tau,
                                  keep_latent = FALSE) {
    m <- length(shapes)
    if (length(scales) != m) stop("shapes and scales must have same length")
    if (any(shapes <= 0) || any(scales <= 0)) stop("shapes and scales must be positive")
    if (nrow(P) != m || ncol(P) != m) stop("P must be m x m")
    if (!all(abs(diag(P) - 1) < 1e-10)) stop("Diagonal of P must be 1 (C1)")

    # Generate component lifetimes
    Tm <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
        Tm[, j] <- rweibull(n, shape = shapes[j], scale = scales[j])
    }

    # System lifetime and failed component
    sys_time <- apply(Tm, 1, min)
    k <- apply(Tm, 1, which.min)

    # Right censoring
    tau <- rep(tau, length.out = n)
    delta <- as.numeric(sys_time <= tau)
    obs_time <- pmin(sys_time, tau)

    # Generate candidate sets using general Bernoulli model
    C <- matrix(FALSE, nrow = n, ncol = m)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            ki <- k[i]
            for (j in seq_len(m)) {
                if (j == ki) {
                    C[i, j] <- TRUE  # C1: failed component always included
                } else {
                    C[i, j] <- runif(1) <= P[j, ki]
                }
            }
        }
        # Censored observations have empty candidate set
    }

    result <- list(
        t = obs_time,
        delta = delta,
        C = C,
        k = k,
        m = m,
        shapes = shapes,
        scales = scales
    )

    if (keep_latent) {
        result$Tm <- Tm
        result$P <- P
    }

    result
}

# -----------------------------------------------------------------------------
# Data Frame Interface
# -----------------------------------------------------------------------------

#' Wrapper: Log-likelihood (Weibull C1-C3) from data frame
#'
#' @param df Data frame with columns t, delta, and x1, x2, ..., xm
#' @param P Inclusion probability matrix (m x m)
#' @param tvar Column name for system lifetime (default "t")
#' @param deltavar Column name for censoring indicator (default "delta")
#' @param setvar Column prefix for candidate set (default "x")
#' @return Log-likelihood function
#' @export
loglik_wei_series_c1_c3_df <- function(df, P, tvar = "t",
                                        deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    loglik_wei_series_c1_c3(t, C, delta, P)
}

#' Wrapper: Score (Weibull C1-C3) from data frame
#' @inheritParams loglik_wei_series_c1_c3_df
#' @return Score function
#' @export
score_wei_series_c1_c3_df <- function(df, P, tvar = "t",
                                       deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    score_wei_series_c1_c3(t, C, delta, P)
}

#' Wrapper: FIM (Weibull C1-C3) from data frame
#' @inheritParams loglik_wei_series_c1_c3_df
#' @return FIM function
#' @export
fim_wei_series_c1_c3_df <- function(df, P, tvar = "t",
                                     deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    fim_wei_series_c1_c3(t, C, delta, P)
}

#' Wrapper: MLE (Weibull C1-C3) from data frame
#' @inheritParams loglik_wei_series_c1_c3_df
#' @param theta0 Initial parameters
#' @return MLE result list
#' @export
mle_wei_series_c1_c3_df <- function(df, P, tvar = "t", deltavar = "delta",
                                     setvar = "x", theta0 = NULL) {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    mle_wei_series_c1_c3(t, C, delta, P, theta0)
}
