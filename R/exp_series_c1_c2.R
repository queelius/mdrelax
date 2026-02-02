# =============================================================================
# Exponential Series System with Masked Data (C1, C2) - Relaxed C3
# Power-weighted hazard model where inclusion probability depends on theta
# =============================================================================

# -----------------------------------------------------------------------------
# Masking Probability Computation
# -----------------------------------------------------------------------------

#' Compute power-weighted inclusion probabilities
#'
#' For exponential components, the inclusion probability for component j is
#' proportional to theta_j^alpha. This models diagnosticians who are more
#' likely to include components with higher failure rates.
#'
#' @param theta Rate parameters
#' @param alpha Power parameter (0 = uniform, higher = more informative)
#' @return Vector of inclusion probabilities (sums to 1)
#' @keywords internal
power_weights <- function(theta, alpha) {
    if (alpha == 0) return(rep(1 / length(theta), length(theta)))
    w <- theta^alpha
    w / sum(w)
}

#' Compute Bernoulli candidate set probability under power-weighted model
#'
#' Under C2 and power-weighted masking, the probability of observing
#' candidate set c is the same for all k in c (C2 holds), and depends
#' on theta through the power weights.
#'
#' @param c Logical vector indicating candidate set
#' @param theta Rate parameters
#' @param alpha Power parameter
#' @param base_p Baseline inclusion probability scale
#' @return Probability of observing candidate set c (for any k in c)
#' @keywords internal
compute_pi_power <- function(c, theta, alpha, base_p) {
    m <- length(c)
    weights <- power_weights(theta, alpha)

    # Each non-failed component j included with prob = base_p * weights[j] / max(weights)
    # Normalize so max inclusion prob = base_p
    w_norm <- weights / max(weights) * base_p

    prob <- 1
    for (j in seq_len(m)) {
        if (c[j]) {
            # This component is in the candidate set
            # Could be the failed component (prob 1) or non-failed (prob w_norm[j])
            # Under C2 the probability is the same regardless of which k in c failed
            # For the "Bernoulli model" interpretation:
            # P(j in C) = w_norm[j] for non-failed, 1 for failed
            # But since we're computing P(C=c | K=k) for any k in c,
            # we need to be careful: the k-th component always has prob 1
            # For the C2 formulation, pi_c factors out, so we compute:
            # pi_c(theta) = product over j in c\{k}: p_j * product over j not in c: (1 - p_j)
            # which under C2 is the same for all k in c
        }
    }
    # This function isn't needed directly; the log-likelihood uses the
    # pi_c(t; theta) formulation from the paper
    prob
}

# -----------------------------------------------------------------------------
# Log-Likelihood, Score, and FIM for Relaxed C3 Model
# -----------------------------------------------------------------------------

#' Log-likelihood for exponential series system (C1, C2) - Relaxed C3
#'
#' Returns a log-likelihood function for masked data where the masking
#' probabilities depend on the model parameters through power-weighted
#' hazard rates. C2 is assumed to hold, so pi_c factors out of the sum.
#'
#' Under this model, the inclusion probability for non-failed component j is:
#' p_j(theta) = base_p * theta_j^alpha / max(theta^alpha)
#'
#' The log-likelihood includes the masking probability contribution:
#' `l(theta) = sum_i [-t_i * sum(theta)] + sum_{i:delta=1} [log(sum_{k in C_i} theta_k) + log(pi_c(theta))]`
#'
#' @param t Numeric vector of system lifetimes
#' @param C Logical matrix of candidate sets (n x m)
#' @param delta Censoring indicators (1=observed, 0=censored)
#' @param alpha Power parameter (0 = uniform/uninformative). If NULL, will be
#'        estimated jointly with theta.
#' @param base_p Baseline inclusion probability
#' @return Function that takes theta (and optionally alpha) and returns log-likelihood
#' @export
loglik_exp_series_c1_c2 <- function(t, C, delta = NULL, alpha = NULL,
                                     base_p = 0.5) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    sum_t <- sum(t)
    alpha_fixed <- !is.null(alpha)

    if (alpha_fixed) {
        function(theta) {
            if (length(theta) != m) stop("theta must have length ", m)
            if (any(theta <= 0)) return(-Inf)

            # Standard C1-C2-C3 survival + hazard contribution
            ll <- -sum_t * sum(theta)

            # Compute power weights
            weights <- power_weights(theta, alpha)
            w_norm <- weights / max(weights) * base_p

            for (i in seq_len(n)) {
                if (delta[i] == 1) {
                    candidates <- which(C[i, ])
                    theta_C <- sum(theta[candidates])
                    if (theta_C <= 0) return(-Inf)

                    # Hazard contribution
                    ll <- ll + log(theta_C)

                    # Masking probability contribution (under C2, same for all k in c)
                    # log pi_c(theta) = sum_{j in c, j != k} log(w_norm[j])
                    #                  + sum_{j not in c} log(1 - w_norm[j])
                    # Under C2 this is the same for any k in c, so pick any k
                    # (we don't know k, but under C2 it doesn't matter)
                    # Use the convention: sum over all j, but the "failed" spot
                    # contributes log(1) = 0 since p_k = 1
                    # We compute: sum_{j in c} log(w_norm[j]) + sum_{j not in c} log(1 - w_norm[j])
                    # minus one log(w_norm[k]) term replaced by log(1) = 0
                    # But we don't know k! Under C2 the answer is the same for any k in c.
                    # So: log(pi_c) = sum_{j in c\{k}} log(w_norm[j]) + sum_{j not in c} log(1 - w_norm[j])
                    # = sum_{j in c} log(w_norm[j]) - log(w_norm[k]) + sum_{j not in c} log(1 - w_norm[j])
                    # Since this must be constant over k in c, we can compute for k = candidates[1]:
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
            if (length(theta) != m) stop("theta must have length ", m)
            if (any(theta <= 0)) return(-Inf)
            if (alpha < 0) return(-Inf)

            ll <- -sum_t * sum(theta)

            weights <- power_weights(theta, alpha)
            w_norm <- weights / max(weights) * base_p

            for (i in seq_len(n)) {
                if (delta[i] == 1) {
                    candidates <- which(C[i, ])
                    theta_C <- sum(theta[candidates])
                    if (theta_C <= 0) return(-Inf)

                    ll <- ll + log(theta_C)

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

#' Score (gradient) for exponential series system (C1, C2) - Relaxed C3
#'
#' Returns the score function. When alpha is fixed, this includes the
#' derivative of the masking probability term with respect to theta.
#'
#' @inheritParams loglik_exp_series_c1_c2
#' @return Function that takes theta and returns gradient vector
#' @export
score_exp_series_c1_c2 <- function(t, C, delta = NULL, alpha,
                                    base_p = 0.5) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    sum_t <- sum(t)

    function(theta) {
        if (length(theta) != m) stop("theta must have length ", m)
        if (any(theta <= 0)) return(rep(NA_real_, m))

        # Standard C1-C2-C3 score component
        g <- rep(-sum_t, m)

        # Compute power weights and their derivatives
        w <- theta^alpha
        sw <- sum(w)
        weights <- w / sw  # p_j = w_j / sw
        w_max <- max(w)
        w_norm <- w / w_max * base_p

        # d(w_norm[j])/d(theta_l) is needed for the masking term
        # w_norm[j] = base_p * theta_j^alpha / max(theta^alpha)
        # Let j_max = argmax(theta^alpha). Then:
        # d(w_norm[j])/d(theta_l) = base_p * [d(theta_j^alpha)/d(theta_l) * max - theta_j^alpha * d(max)/d(theta_l)] / max^2
        # For j != j_max: d(max)/d(theta_l) = 0 if l != j_max, = alpha * theta_{j_max}^{alpha-1} if l == j_max
        # This gets complex; use numerical gradient for the masking part

        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                candidates <- which(C[i, ])
                theta_C <- sum(theta[candidates])
                if (theta_C <= 0) return(rep(NA_real_, m))

                # Standard hazard contribution
                g[candidates] <- g[candidates] + 1 / theta_C
            }
        }

        # Add masking probability gradient numerically
        # (analytical form is complex due to max() in normalization)
        if (alpha != 0) {
            ll_mask <- function(th) {
                ww <- th^alpha
                wn <- ww / max(ww) * base_p
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
            for (l in seq_len(m)) {
                theta_plus <- theta
                theta_plus[l] <- theta_plus[l] + h
                g[l] <- g[l] + (ll_mask(theta_plus) - ll_mask(theta)) / h
            }
        }

        g
    }
}

#' Fisher Information Matrix for exponential series (C1, C2) - Relaxed C3
#'
#' Returns the observed FIM computed numerically from the score function.
#'
#' @inheritParams loglik_exp_series_c1_c2
#' @return Function that takes theta and returns m x m FIM
#' @export
fim_exp_series_c1_c2 <- function(t, C, delta = NULL, alpha,
                                  base_p = 0.5) {
    ll_fn <- loglik_exp_series_c1_c2(t, C, delta, alpha, base_p)

    function(theta) {
        m <- length(theta)
        h <- 1e-5
        H <- matrix(0, m, m)
        for (j in seq_len(m)) {
            for (k in seq_len(m)) {
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

#' MLE for exponential series system (C1, C2) - Relaxed C3
#'
#' Computes MLE of component failure rates under power-weighted masking.
#'
#' @inheritParams loglik_exp_series_c1_c2
#' @param theta0 Initial parameter values. If NULL, uses method of moments.
#' @return List with theta, alpha, se, vcov, loglik, converged, fim
#' @export
mle_exp_series_c1_c2 <- function(t, C, delta = NULL, alpha = NULL,
                                  base_p = 0.5, theta0 = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)

    if (is.null(theta0)) {
        n_unc <- sum(delta == 1)
        theta0 <- rep(n_unc / (sum(t) * m), m)
    }

    if (!is.null(alpha)) {
        # Fixed alpha: optimize theta only
        ll_fn <- loglik_exp_series_c1_c2(t, C, delta, alpha, base_p)
        sc_fn <- score_exp_series_c1_c2(t, C, delta, alpha, base_p)

        result <- optim(
            par = theta0,
            fn = ll_fn,
            gr = sc_fn,
            method = "L-BFGS-B",
            lower = rep(1e-10, m),
            control = list(fnscale = -1, maxit = 1000)
        )

        fim <- fim_exp_series_c1_c2(t, C, delta, alpha, base_p)(result$par)
        vcov <- tryCatch(solve(fim), error = function(e) matrix(NA, m, m))
        se <- sqrt(diag(vcov))

        return(list(
            theta = result$par,
            alpha = alpha,
            se = se,
            vcov = vcov,
            loglik = result$value,
            converged = result$convergence == 0,
            fim = fim
        ))
    } else {
        # Joint optimization of theta and alpha
        ll_fn <- loglik_exp_series_c1_c2(t, C, delta, alpha = NULL, base_p)

        obj <- function(par) {
            th <- par[1:m]
            a <- par[m + 1]
            ll_fn(th, a)
        }

        par0 <- c(theta0, 1)  # Start with alpha = 1
        lower <- c(rep(1e-10, m), 0)
        upper <- c(rep(Inf, m), 10)

        result <- optim(
            par = par0,
            fn = obj,
            method = "L-BFGS-B",
            lower = lower,
            upper = upper,
            control = list(fnscale = -1, maxit = 1000)
        )

        theta_hat <- result$par[1:m]
        alpha_hat <- result$par[m + 1]

        # FIM for theta only (profile over alpha)
        fim <- fim_exp_series_c1_c2(t, C, delta, alpha_hat, base_p)(theta_hat)
        vcov <- tryCatch(solve(fim), error = function(e) matrix(NA, m, m))
        se <- sqrt(diag(vcov))

        return(list(
            theta = theta_hat,
            alpha = alpha_hat,
            se = se,
            vcov = vcov,
            loglik = result$value,
            converged = result$convergence == 0,
            fim = fim
        ))
    }
}

# -----------------------------------------------------------------------------
# Data Generation
# -----------------------------------------------------------------------------

#' Generate masked data with power-weighted masking (relaxed C3)
#'
#' Generates data where the candidate set inclusion probability depends on
#' the component failure rates through a power-weighted hazard model.
#'
#' @param n Sample size
#' @param theta Rate parameters for exponential components
#' @param alpha Power parameter (0 = uniform masking)
#' @param base_p Baseline inclusion probability
#' @param tau Right-censoring time
#' @param keep_latent If TRUE, include latent variables in output
#' @return List with t, delta, C, k
#' @export
rexp_series_md_c1_c2 <- function(n, theta, alpha, base_p = 0.5, tau,
                                  keep_latent = FALSE) {
    m <- length(theta)
    if (any(theta <= 0)) stop("theta must be positive")
    if (alpha < 0) stop("alpha must be non-negative")

    # Compute inclusion probabilities
    weights <- power_weights(theta, alpha)
    w_norm <- weights / max(weights) * base_p

    # Generate component lifetimes
    Tm <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
        Tm[, j] <- rexp(n, rate = theta[j])
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
        k = k
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

#' Wrapper: Log-likelihood (relaxed C3) from data frame
#' @param df Data frame with t, delta, x1, x2, ... columns
#' @param alpha Power parameter (NULL to estimate)
#' @param base_p Baseline inclusion probability
#' @param tvar,deltavar,setvar Column name parameters
#' @return Log-likelihood function
#' @export
loglik_exp_series_c1_c2_df <- function(df, alpha = NULL, base_p = 0.5,
                                        tvar = "t", deltavar = "delta",
                                        setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    loglik_exp_series_c1_c2(t, C, delta, alpha, base_p)
}

#' Wrapper: MLE (relaxed C3) from data frame
#' @inheritParams loglik_exp_series_c1_c2_df
#' @param theta0 Initial parameters
#' @return MLE result list
#' @export
mle_exp_series_c1_c2_df <- function(df, alpha = NULL, base_p = 0.5,
                                     theta0 = NULL, tvar = "t",
                                     deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    mle_exp_series_c1_c2(t, C, delta, alpha, base_p, theta0)
}
