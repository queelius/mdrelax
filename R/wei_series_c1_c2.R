# =============================================================================
# Weibull Series System with Masked Data (C1, C2) - Relaxed C3
# Parameter-dependent masking where inclusion prob depends on theta via power weights
# =============================================================================
#
# Parameter convention: theta = c(k1, lambda1, k2, lambda2, ..., km, lambdam)
# where k_j = shape and lambda_j = scale for component j.
#
# For Weibull components, the "failure proneness" measure is the characteristic
# hazard rate k_j/lambda_j (the hazard function evaluated at t = lambda_j).
# The power-weighted model uses:
#   p_j(theta) = base_p * (k_j/lambda_j)^alpha / max((k_l/lambda_l)^alpha)

# -----------------------------------------------------------------------------
# Power Weights for Weibull
# -----------------------------------------------------------------------------

#' Compute power weights for Weibull components
#'
#' For Weibull components, the "failure proneness" measure is the characteristic
#' hazard rate k_j/lambda_j. These weights determine inclusion probabilities
#' under the relaxed C3 model.
#'
#' @param shapes Shape parameters (k_1, ..., k_m)
#' @param scales Scale parameters (lambda_1, ..., lambda_m)
#' @param alpha Power parameter (0 = uniform, higher = more informative)
#' @return Vector of weights that sum to 1
#' @export
wei_power_weights <- function(shapes, scales, alpha) {
    m <- length(shapes)
    if (alpha == 0) return(rep(1/m, m))

    # Characteristic hazard rate: k/lambda = hazard at t = lambda
    char_haz <- shapes / scales
    w <- char_haz^alpha
    w / sum(w)
}

#' Compute normalized inclusion probabilities for Weibull power model
#'
#' @param shapes Shape parameters
#' @param scales Scale parameters
#' @param alpha Power parameter
#' @param base_p Baseline inclusion probability
#' @return Vector of inclusion probabilities (max = base_p)
#' @keywords internal
wei_inclusion_probs <- function(shapes, scales, alpha, base_p) {
    weights <- wei_power_weights(shapes, scales, alpha)
    weights / max(weights) * base_p
}

# -----------------------------------------------------------------------------
# Log-Likelihood
# -----------------------------------------------------------------------------

#' Log-likelihood for Weibull series system (C1, C2) - Relaxed C3
#'
#' Returns a log-likelihood function for masked data where the masking
#' probabilities depend on the model parameters through power-weighted
#' characteristic hazard rates. C2 is assumed to hold, so pi_c factors out.
#'
#' Under this model, the inclusion probability for non-failed component j is:
#' `p_j(theta) = base_p * (k_j/lambda_j)^alpha / max((k/lambda)^alpha)`
#'
#' The log-likelihood includes the masking probability contribution:
#' `l(theta) = sum_i [-sum_j (t_i/lambda_j)^k_j] + sum_{i:delta=1} [log(sum_{j in C_i} h_j(t_i)) + log(pi_c(theta))]`
#'
#' @param t Numeric vector of system lifetimes
#' @param C Logical matrix of candidate sets (n x m)
#' @param delta Censoring indicators (1=observed, 0=censored)
#' @param alpha Power parameter (0 = uniform/uninformative). If NULL, will be
#'        estimated jointly with theta.
#' @param base_p Baseline inclusion probability
#' @return Function that takes theta (and optionally alpha) and returns log-likelihood
#' @export
loglik_wei_series_c1_c2 <- function(t, C, delta = NULL, alpha = NULL,
                                     base_p = 0.5) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    alpha_fixed <- !is.null(alpha)

    if (alpha_fixed) {
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

            # Compute inclusion probabilities
            w_norm <- wei_inclusion_probs(shapes, scales, alpha, base_p)

            # Hazard and masking contributions for uncensored observations
            for (i in seq_len(n)) {
                if (delta[i] == 1) {
                    candidates <- which(C[i, ])

                    # Hazard sum (no pi weighting under C2)
                    hazard_sum <- 0
                    for (j in candidates) {
                        h_j <- (shapes[j] / scales[j]) * (t[i] / scales[j])^(shapes[j] - 1)
                        hazard_sum <- hazard_sum + h_j
                    }
                    if (hazard_sum <= 0) return(-Inf)
                    ll <- ll + log(hazard_sum)

                    # Masking probability contribution
                    # Under C2: log(pi_c) = sum_{j in c\{k}} log(w_norm[j]) + sum_{j not in c} log(1 - w_norm[j])
                    # Since this is same for all k in c, pick k = candidates[1]
                    k <- candidates[1]
                    log_pi <- 0
                    for (j in seq_len(m)) {
                        if (j == k) next
                        if (C[i, j]) {
                            if (w_norm[j] <= 0) return(-Inf)
                            log_pi <- log_pi + log(w_norm[j])
                        } else {
                            if (w_norm[j] >= 1) return(-Inf)
                            log_pi <- log_pi + log(1 - w_norm[j])
                        }
                    }
                    ll <- ll + log_pi
                }
            }
            ll
        }
    } else {
        # alpha is unknown - jointly optimize
        function(theta, alpha) {
            params <- unpack_wei_params(theta)
            if (params$m != m) stop("theta implies ", params$m, " components but C has ", m)
            shapes <- params$shapes
            scales <- params$scales
            if (any(shapes <= 0) || any(scales <= 0)) return(-Inf)
            if (alpha < 0) return(-Inf)

            ll <- 0
            for (j in seq_len(m)) {
                ll <- ll - sum((t / scales[j])^shapes[j])
            }

            w_norm <- wei_inclusion_probs(shapes, scales, alpha, base_p)

            for (i in seq_len(n)) {
                if (delta[i] == 1) {
                    candidates <- which(C[i, ])

                    hazard_sum <- 0
                    for (j in candidates) {
                        h_j <- (shapes[j] / scales[j]) * (t[i] / scales[j])^(shapes[j] - 1)
                        hazard_sum <- hazard_sum + h_j
                    }
                    if (hazard_sum <= 0) return(-Inf)
                    ll <- ll + log(hazard_sum)

                    k <- candidates[1]
                    log_pi <- 0
                    for (j in seq_len(m)) {
                        if (j == k) next
                        if (C[i, j]) {
                            if (w_norm[j] <= 0) return(-Inf)
                            log_pi <- log_pi + log(w_norm[j])
                        } else {
                            if (w_norm[j] >= 1) return(-Inf)
                            log_pi <- log_pi + log(1 - w_norm[j])
                        }
                    }
                    ll <- ll + log_pi
                }
            }
            ll
        }
    }
}

#' Score (gradient) for Weibull series system (C1, C2) - Relaxed C3
#'
#' Returns the score function. When alpha is fixed, this includes the
#' derivative of the masking probability term with respect to theta.
#'
#' @inheritParams loglik_wei_series_c1_c2
#' @return Function that takes theta and returns gradient vector
#' @export
score_wei_series_c1_c2 <- function(t, C, delta = NULL, alpha,
                                    base_p = 0.5) {
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
            log_z[t == 0] <- 0

            # d/d(k_j) of -sum (t/lambda_j)^k_j = -sum z^k * log(z)
            g[2*j - 1] <- -sum(z_k * log_z)

            # d/d(lambda_j) of -sum (t/lambda_j)^k_j = sum k_j/lambda_j * z^k_j
            g[2*j] <- sum(shapes[j] / scales[j] * z_k)
        }

        # Hazard derivatives for uncensored observations
        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                candidates <- which(C[i, ])

                # Compute hazards
                h_vals <- numeric(m)
                for (j in candidates) {
                    z <- t[i] / scales[j]
                    h_vals[j] <- (shapes[j] / scales[j]) * z^(shapes[j] - 1)
                }
                hazard_sum <- sum(h_vals[candidates])
                if (hazard_sum <= 0) return(rep(NA_real_, 2*m))

                # Hazard contribution to score
                for (j in candidates) {
                    z <- t[i] / scales[j]
                    log_z <- log(z)
                    if (t[i] == 0) log_z <- 0

                    dh_dk <- h_vals[j] * (1 / shapes[j] + log_z)
                    dh_dlam <- -h_vals[j] * shapes[j] / scales[j]

                    g[2*j - 1] <- g[2*j - 1] + dh_dk / hazard_sum
                    g[2*j] <- g[2*j] + dh_dlam / hazard_sum
                }
            }
        }

        # Add masking probability gradient numerically
        # (analytical form is complex due to max() in normalization)
        if (alpha != 0) {
            ll_mask <- function(th) {
                params <- unpack_wei_params(th)
                wn <- wei_inclusion_probs(params$shapes, params$scales, alpha, base_p)
                val <- 0
                for (i in seq_len(n)) {
                    if (delta[i] == 1) {
                        cands <- which(C[i, ])
                        k <- cands[1]
                        for (j in seq_len(m)) {
                            if (j == k) next
                            if (C[i, j]) {
                                if (wn[j] <= 0) return(-Inf)
                                val <- val + log(wn[j])
                            } else {
                                if (wn[j] >= 1) return(-Inf)
                                val <- val + log(1 - wn[j])
                            }
                        }
                    }
                }
                val
            }

            h <- 1e-7
            for (l in seq_len(2*m)) {
                theta_plus <- theta
                theta_plus[l] <- theta_plus[l] + h
                g[l] <- g[l] + (ll_mask(theta_plus) - ll_mask(theta)) / h
            }
        }

        g
    }
}

#' Fisher Information Matrix for Weibull series (C1, C2) - Relaxed C3
#'
#' Computed numerically as the negative Hessian of the log-likelihood.
#'
#' @inheritParams loglik_wei_series_c1_c2
#' @return Function that takes theta and returns (2m x 2m) FIM
#' @export
fim_wei_series_c1_c2 <- function(t, C, delta = NULL, alpha,
                                  base_p = 0.5) {
    ll_fn <- loglik_wei_series_c1_c2(t, C, delta, alpha, base_p)

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

#' MLE for Weibull series system (C1, C2) - Relaxed C3
#'
#' Computes MLE of Weibull parameters under power-weighted masking.
#'
#' @inheritParams loglik_wei_series_c1_c2
#' @param theta0 Initial parameter values. If NULL, uses heuristic.
#' @return List with theta, shapes, scales, alpha, se, vcov, loglik, converged, fim, m
#' @export
mle_wei_series_c1_c2 <- function(t, C, delta = NULL, alpha = NULL,
                                  base_p = 0.5, theta0 = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)

    if (is.null(theta0)) {
        n_unc <- sum(delta == 1)
        total_rate <- n_unc / sum(t)
        rate_each <- total_rate / m
        theta0 <- rep(c(1, 1 / rate_each), m)
    }

    if (!is.null(alpha)) {
        # Fixed alpha: optimize theta only
        ll_fn <- loglik_wei_series_c1_c2(t, C, delta, alpha, base_p)
        sc_fn <- score_wei_series_c1_c2(t, C, delta, alpha, base_p)

        result <- optim(
            par = theta0,
            fn = ll_fn,
            gr = sc_fn,
            method = "L-BFGS-B",
            lower = rep(1e-6, 2*m),
            control = list(fnscale = -1, maxit = 1000)
        )

        fim <- fim_wei_series_c1_c2(t, C, delta, alpha, base_p)(result$par)
        vcov <- tryCatch(solve(fim), error = function(e) matrix(NA, 2*m, 2*m))
        se <- sqrt(diag(vcov))

        params <- unpack_wei_params(result$par)

        return(list(
            theta = result$par,
            shapes = params$shapes,
            scales = params$scales,
            alpha = alpha,
            se = se,
            vcov = vcov,
            loglik = result$value,
            converged = result$convergence == 0,
            fim = fim,
            n = n,
            m = m
        ))
    } else {
        # Joint optimization of theta and alpha
        ll_fn <- loglik_wei_series_c1_c2(t, C, delta, alpha = NULL, base_p)

        obj <- function(par) {
            th <- par[1:(2*m)]
            a <- par[2*m + 1]
            ll_fn(th, a)
        }

        par0 <- c(theta0, 1)  # Start with alpha = 1
        lower <- c(rep(1e-6, 2*m), 0)
        upper <- c(rep(Inf, 2*m), 10)

        result <- optim(
            par = par0,
            fn = obj,
            method = "L-BFGS-B",
            lower = lower,
            upper = upper,
            control = list(fnscale = -1, maxit = 1000)
        )

        theta_hat <- result$par[1:(2*m)]
        alpha_hat <- result$par[2*m + 1]

        # FIM for theta only (profile over alpha)
        fim <- fim_wei_series_c1_c2(t, C, delta, alpha_hat, base_p)(theta_hat)
        vcov <- tryCatch(solve(fim), error = function(e) matrix(NA, 2*m, 2*m))
        se <- sqrt(diag(vcov))

        params <- unpack_wei_params(theta_hat)

        return(list(
            theta = theta_hat,
            shapes = params$shapes,
            scales = params$scales,
            alpha = alpha_hat,
            se = se,
            vcov = vcov,
            loglik = result$value,
            converged = result$convergence == 0,
            fim = fim,
            n = n,
            m = m
        ))
    }
}

# -----------------------------------------------------------------------------
# Data Generation
# -----------------------------------------------------------------------------

#' Generate masked data with power-weighted masking (Weibull, relaxed C3)
#'
#' Generates data where the candidate set inclusion probability depends on
#' the component parameters through a power-weighted characteristic hazard model.
#'
#' @param n Sample size
#' @param shapes Shape parameters for Weibull components
#' @param scales Scale parameters for Weibull components
#' @param alpha Power parameter (0 = uniform masking)
#' @param base_p Baseline inclusion probability
#' @param tau Right-censoring time
#' @param keep_latent If TRUE, include latent variables in output
#' @return List with t, delta, C, k, m, shapes, scales
#' @export
rwei_series_md_c1_c2 <- function(n, shapes, scales, alpha, base_p = 0.5, tau,
                                  keep_latent = FALSE) {
    m <- length(shapes)
    if (length(scales) != m) stop("shapes and scales must have same length")
    if (any(shapes <= 0) || any(scales <= 0)) stop("shapes and scales must be positive")
    if (alpha < 0) stop("alpha must be non-negative")

    # Compute inclusion probabilities
    w_norm <- wei_inclusion_probs(shapes, scales, alpha, base_p)

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

    # Generate candidate sets
    C <- matrix(FALSE, nrow = n, ncol = m)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            ki <- k[i]
            for (j in seq_len(m)) {
                if (j == ki) {
                    C[i, j] <- TRUE  # C1: failed component always in
                } else {
                    C[i, j] <- runif(1) <= w_norm[j]
                }
            }
        }
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
        result$alpha <- alpha
        result$base_p <- base_p
        result$w_norm <- w_norm
    }

    result
}

# -----------------------------------------------------------------------------
# Data Frame Interface
# -----------------------------------------------------------------------------

#' Wrapper: Log-likelihood (Weibull C1-C2) from data frame
#'
#' @param df Data frame with columns t, delta, and x1, x2, ..., xm
#' @param alpha Power parameter (NULL to estimate)
#' @param base_p Baseline inclusion probability
#' @param tvar Column name for system lifetime (default "t")
#' @param deltavar Column name for censoring indicator (default "delta")
#' @param setvar Column prefix for candidate set (default "x")
#' @return Log-likelihood function
#' @export
loglik_wei_series_c1_c2_df <- function(df, alpha = NULL, base_p = 0.5,
                                        tvar = "t", deltavar = "delta",
                                        setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    loglik_wei_series_c1_c2(t, C, delta, alpha, base_p)
}

#' Wrapper: Score (Weibull C1-C2) from data frame
#' @inheritParams loglik_wei_series_c1_c2_df
#' @return Score function
#' @export
score_wei_series_c1_c2_df <- function(df, alpha, base_p = 0.5,
                                       tvar = "t", deltavar = "delta",
                                       setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    score_wei_series_c1_c2(t, C, delta, alpha, base_p)
}

#' Wrapper: FIM (Weibull C1-C2) from data frame
#' @inheritParams loglik_wei_series_c1_c2_df
#' @return FIM function
#' @export
fim_wei_series_c1_c2_df <- function(df, alpha, base_p = 0.5,
                                     tvar = "t", deltavar = "delta",
                                     setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    fim_wei_series_c1_c2(t, C, delta, alpha, base_p)
}

#' Wrapper: MLE (Weibull C1-C2) from data frame
#' @inheritParams loglik_wei_series_c1_c2_df
#' @param theta0 Initial parameters
#' @return MLE result list
#' @export
mle_wei_series_c1_c2_df <- function(df, alpha = NULL, base_p = 0.5,
                                     theta0 = NULL, tvar = "t",
                                     deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    mle_wei_series_c1_c2(t, C, delta, alpha, base_p, theta0)
}
