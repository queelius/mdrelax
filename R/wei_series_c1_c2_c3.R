# =============================================================================
# Weibull Series System with Masked Data (C1, C2, C3)
# Dependency-free implementation for simulations
# =============================================================================
#
# Parameter convention: theta = c(k1, beta1, k2, beta2, ..., km, betam)
# where k_j = shape and beta_j = scale for component j.
#
# For m components, theta has length 2*m.
# Component j uses theta[2*j - 1] = shape, theta[2*j] = scale.

# -----------------------------------------------------------------------------
# Parameter Helpers
# -----------------------------------------------------------------------------

#' Extract shapes and scales from Weibull parameter vector
#'
#' @param theta Parameter vector c(k1, beta1, k2, beta2, ..., km, betam)
#' @return List with shapes and scales vectors
#' @keywords internal
unpack_wei_params <- function(theta) {
    m <- length(theta) / 2
    if (m != floor(m)) stop("theta must have even length (shape-scale pairs)")
    shapes <- theta[seq(1, 2*m, by = 2)]
    scales <- theta[seq(2, 2*m, by = 2)]
    list(shapes = shapes, scales = scales, m = as.integer(m))
}

#' Pack shapes and scales into Weibull parameter vector
#'
#' @param shapes Shape parameters
#' @param scales Scale parameters
#' @return Parameter vector c(k1, beta1, k2, beta2, ..., km, betam)
#' @keywords internal
pack_wei_params <- function(shapes, scales) {
    m <- length(shapes)
    theta <- numeric(2 * m)
    theta[seq(1, 2*m, by = 2)] <- shapes
    theta[seq(2, 2*m, by = 2)] <- scales
    theta
}

# -----------------------------------------------------------------------------
# Log-Likelihood
# -----------------------------------------------------------------------------

#' Log-likelihood for Weibull series system with masked data (C1, C2, C3)
#'
#' Returns a log-likelihood function for masked data from a Weibull series
#' system under conditions C1, C2, C3.
#'
#' The log-likelihood is:
#' `l(theta) = sum_i [-sum_j (t_i/beta_j)^k_j] + sum_{i:delta_i=1} log(sum_{j in C_i} (k_j/beta_j)(t_i/beta_j)^(k_j-1))`
#'
#' @param t Numeric vector of system lifetimes
#' @param C Logical matrix of candidate sets (n x m)
#' @param delta Censoring indicators (1=observed, 0=censored)
#' @return Function that takes theta = c(k1,b1,k2,b2,...) and returns log-likelihood
#' @export
loglik_wei_series <- function(t, C, delta = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    function(theta) {
        params <- unpack_wei_params(theta)
        if (params$m != m) stop("theta implies ", params$m, " components but C has ", m)
        shapes <- params$shapes
        scales <- params$scales
        if (any(shapes <= 0) || any(scales <= 0)) return(-Inf)

        # Survival contribution: -sum_i sum_j (t_i / beta_j)^k_j
        ll <- 0
        for (j in seq_len(m)) {
            ll <- ll - sum((t / scales[j])^shapes[j])
        }

        # Hazard contribution for uncensored observations
        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                # h_j(t_i) = (k_j / beta_j) * (t_i / beta_j)^(k_j - 1)
                hazard_sum <- 0
                for (j in which(C[i, ])) {
                    hazard_sum <- hazard_sum +
                        (shapes[j] / scales[j]) * (t[i] / scales[j])^(shapes[j] - 1)
                }
                if (hazard_sum <= 0) return(-Inf)
                ll <- ll + log(hazard_sum)
            }
        }
        ll
    }
}

#' Score (gradient) for Weibull series system (C1, C2, C3)
#'
#' @inheritParams loglik_wei_series
#' @return Function that takes theta and returns gradient vector
#' @export
score_wei_series <- function(t, C, delta = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

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
            log_z[t == 0] <- 0  # Handle t=0

            # d/d(k_j) of -sum (t/beta_j)^k_j = -sum z^k * log(z)
            g[2*j - 1] <- -sum(z_k * log_z)

            # d/d(beta_j) of -sum (t/beta_j)^k_j = sum k_j/beta_j * z^k_j
            g[2*j] <- sum(shapes[j] / scales[j] * z_k)
        }

        # Hazard derivatives for uncensored observations
        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                candidates <- which(C[i, ])

                # Compute hazard for each candidate
                h_vals <- numeric(m)
                for (j in candidates) {
                    z <- t[i] / scales[j]
                    h_vals[j] <- (shapes[j] / scales[j]) * z^(shapes[j] - 1)
                }
                hazard_sum <- sum(h_vals[candidates])
                if (hazard_sum <= 0) return(rep(NA_real_, 2*m))

                for (j in candidates) {
                    z <- t[i] / scales[j]
                    log_z <- log(z)
                    if (t[i] == 0) log_z <- 0

                    # d(h_j)/d(k_j) = h_j * (1/k_j + log(z))
                    dh_dk <- h_vals[j] * (1 / shapes[j] + log_z)

                    # d(h_j)/d(beta_j) = -h_j * k_j / beta_j
                    dh_db <- -h_vals[j] * shapes[j] / scales[j]

                    g[2*j - 1] <- g[2*j - 1] + dh_dk / hazard_sum
                    g[2*j] <- g[2*j] + dh_db / hazard_sum
                }
            }
        }
        g
    }
}

#' Fisher Information Matrix for Weibull series system (C1, C2, C3)
#'
#' Computed numerically as the negative Hessian of the log-likelihood.
#'
#' @inheritParams loglik_wei_series
#' @return Function that takes theta and returns (2m x 2m) FIM
#' @export
fim_wei_series <- function(t, C, delta = NULL) {
    ll_fn <- loglik_wei_series(t, C, delta)

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

#' MLE for Weibull series system (C1, C2, C3)
#'
#' @inheritParams loglik_wei_series
#' @param theta0 Initial parameters c(k1,b1,...). If NULL, uses heuristic.
#' @return List with theta, shapes, scales, se, vcov, loglik, converged, fim
#' @export
mle_wei_series <- function(t, C, delta = NULL, theta0 = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)

    ll_fn <- loglik_wei_series(t, C, delta)
    sc_fn <- score_wei_series(t, C, delta)

    # Initialize: use exponential MLE to get rate, then convert
    if (is.null(theta0)) {
        n_unc <- sum(delta == 1)
        total_rate <- n_unc / sum(t)
        rate_each <- total_rate / m
        # Exponential: shape=1, scale=1/rate
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
    fim <- fim_wei_series(t, C, delta)(result$par)
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

#' Generate masked data from Weibull series system
#'
#' @param n Sample size
#' @param shapes Shape parameters (k_1, ..., k_m)
#' @param scales Scale parameters (beta_1, ..., beta_m)
#' @param p Probability each non-failed component is in candidate set
#' @param tau Right-censoring time
#' @param keep_latent If TRUE, include latent component times
#' @return List with t, delta, C, k, m, shapes, scales
#' @export
rwei_series_md <- function(n, shapes, scales, p, tau,
                            keep_latent = FALSE) {
    m <- length(shapes)
    if (length(scales) != m) stop("shapes and scales must have same length")
    if (any(shapes <= 0) || any(scales <= 0)) stop("shapes and scales must be positive")
    if (p < 0 || p > 1) stop("p must be in [0, 1]")

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

    # Generate candidate sets (Bernoulli model, C1-C2-C3)
    C <- matrix(runif(n * m) <= p, nrow = n, ncol = m)
    C[cbind(seq_len(n), k)] <- TRUE
    C[delta == 0, ] <- FALSE

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
    }

    result
}

# -----------------------------------------------------------------------------
# Data Frame Interface
# -----------------------------------------------------------------------------

#' Wrapper: Log-likelihood (Weibull) from data frame
#' @param df Data frame with t, delta, x1, x2, ... columns
#' @param tvar,deltavar,setvar Column name parameters
#' @return Log-likelihood function
#' @export
loglik_wei_series_df <- function(df, tvar = "t", deltavar = "delta",
                                  setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    loglik_wei_series(t, C, delta)
}

#' Wrapper: MLE (Weibull) from data frame
#' @inheritParams loglik_wei_series_df
#' @param theta0 Initial parameters
#' @return MLE result list
#' @export
mle_wei_series_df <- function(df, tvar = "t", deltavar = "delta",
                               setvar = "x", theta0 = NULL) {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    mle_wei_series(t, C, delta, theta0)
}
