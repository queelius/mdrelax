# =============================================================================
# Exponential Series System with Masked Data - Relaxed C1
# P(K in C) < 1: Failed component may not be in candidate set
# =============================================================================

# -----------------------------------------------------------------------------
# Log-Likelihood for Relaxed C1 Model
# -----------------------------------------------------------------------------

#' Log-likelihood for exponential series system with relaxed C1
#'
#' When C1 is relaxed, the failed component may not appear in the candidate
#' set. The P matrix diagonal can be less than 1: `P[k,k] = P(k in C | K=k)`.
#'
#' The likelihood must sum over ALL possible failed components k=1,...,m
#' (not just those in c), since any component could have failed but not
#' been included in the candidate set.
#'
#' @param t Numeric vector of system lifetimes
#' @param C Logical matrix of candidate sets (n x m)
#' @param delta Censoring indicators (1=observed, 0=censored)
#' @param P Inclusion probability matrix (m x m): `P[j,k] = P(j in C | K = k)`.
#'        Diagonal can be < 1 (relaxed C1).
#' @return Function that takes theta and returns log-likelihood
#' @export
loglik_exp_series_relaxed_c1 <- function(t, C, delta = NULL, P) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    sum_t <- sum(t)

    # Precompute pi_k(c_i) for ALL k = 1,...,m (not just k in c_i)
    # pi_k(c) = P(C=c | K=k) = prod_{j in c} P[j,k] * prod_{j not in c} (1 - P[j,k])
    # For k not in c: pi_k(c) includes P[k,k] in the "not in c" product
    pi_list <- vector("list", n)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            pi_vals <- numeric(m)
            for (k in seq_len(m)) {
                prob <- 1
                for (j in seq_len(m)) {
                    if (C[i, j]) {
                        prob <- prob * P[j, k]
                    } else {
                        prob <- prob * (1 - P[j, k])
                    }
                }
                pi_vals[k] <- prob
            }
            pi_list[[i]] <- pi_vals
        }
    }

    function(theta) {
        if (length(theta) != m) stop("theta must have length ", m)
        if (any(theta <= 0)) return(-Inf)

        ll <- -sum_t * sum(theta)

        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                pi_vals <- pi_list[[i]]

                # Sum over ALL k = 1,...,m (not just candidates)
                hazard_sum <- sum(theta * pi_vals)
                if (hazard_sum <= 0) return(-Inf)
                ll <- ll + log(hazard_sum)
            }
        }
        ll
    }
}

#' Score for exponential series system with relaxed C1
#'
#' @inheritParams loglik_exp_series_relaxed_c1
#' @return Function that takes theta and returns gradient vector
#' @export
score_exp_series_relaxed_c1 <- function(t, C, delta = NULL, P) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    sum_t <- sum(t)

    # Precompute pi_k(c_i) for ALL k
    pi_list <- vector("list", n)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            pi_vals <- numeric(m)
            for (k in seq_len(m)) {
                prob <- 1
                for (j in seq_len(m)) {
                    if (C[i, j]) {
                        prob <- prob * P[j, k]
                    } else {
                        prob <- prob * (1 - P[j, k])
                    }
                }
                pi_vals[k] <- prob
            }
            pi_list[[i]] <- pi_vals
        }
    }

    function(theta) {
        if (length(theta) != m) stop("theta must have length ", m)
        if (any(theta <= 0)) return(rep(NA_real_, m))

        g <- rep(-sum_t, m)

        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                pi_vals <- pi_list[[i]]
                denom <- sum(theta * pi_vals)
                if (denom <= 0) return(rep(NA_real_, m))

                # d/d(theta_j) log(sum_k theta_k * pi_k) = pi_j / denom
                g <- g + pi_vals / denom
            }
        }
        g
    }
}

#' FIM for exponential series system with relaxed C1
#'
#' @inheritParams loglik_exp_series_relaxed_c1
#' @return Function that takes theta and returns m x m FIM
#' @export
fim_exp_series_relaxed_c1 <- function(t, C, delta = NULL, P) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    # Precompute pi_k(c_i) for ALL k
    pi_list <- vector("list", n)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            pi_vals <- numeric(m)
            for (k in seq_len(m)) {
                prob <- 1
                for (j in seq_len(m)) {
                    if (C[i, j]) {
                        prob <- prob * P[j, k]
                    } else {
                        prob <- prob * (1 - P[j, k])
                    }
                }
                pi_vals[k] <- prob
            }
            pi_list[[i]] <- pi_vals
        }
    }

    function(theta) {
        if (length(theta) != m) stop("theta must have length ", m)
        if (any(theta <= 0)) return(matrix(NA_real_, m, m))

        I <- matrix(0, m, m)

        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                pi_vals <- pi_list[[i]]
                denom <- sum(theta * pi_vals)
                if (denom <= 0) return(matrix(NA_real_, m, m))

                # -d^2/d(theta_j)d(theta_k) log(sum_l theta_l * pi_l)
                # = pi_j * pi_k / denom^2
                I <- I + outer(pi_vals, pi_vals) / denom^2
            }
        }
        I
    }
}

# -----------------------------------------------------------------------------
# MLE for Relaxed C1
# -----------------------------------------------------------------------------

#' MLE for exponential series system with relaxed C1
#'
#' @inheritParams loglik_exp_series_relaxed_c1
#' @param theta0 Initial parameter values
#' @return List with theta, se, vcov, loglik, converged, fim
#' @export
mle_exp_series_relaxed_c1 <- function(t, C, delta = NULL, P,
                                       theta0 = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)

    ll_fn <- loglik_exp_series_relaxed_c1(t, C, delta, P)
    sc_fn <- score_exp_series_relaxed_c1(t, C, delta, P)
    fim_fn <- fim_exp_series_relaxed_c1(t, C, delta, P)

    if (is.null(theta0)) {
        n_unc <- sum(delta == 1)
        theta0 <- rep(n_unc / (sum(t) * m), m)
    }

    result <- optim(
        par = theta0,
        fn = ll_fn,
        gr = sc_fn,
        method = "L-BFGS-B",
        lower = rep(1e-10, m),
        control = list(fnscale = -1, maxit = 1000)
    )

    fim <- fim_fn(result$par)
    vcov <- tryCatch(solve(fim), error = function(e) matrix(NA, m, m))
    se <- sqrt(diag(vcov))

    list(
        theta = result$par,
        P = P,
        se = se,
        vcov = vcov,
        loglik = result$value,
        converged = result$convergence == 0,
        fim = fim
    )
}

# -----------------------------------------------------------------------------
# Data Generation
# -----------------------------------------------------------------------------

#' Generate masked data with relaxed C1
#'
#' Generates data where the failed component may not be in the candidate set.
#' `P[k,k] < 1` means the failed component k is only included with that probability.
#'
#' @param n Sample size
#' @param theta Rate parameters for exponential components
#' @param P Inclusion probability matrix (m x m). Diagonal can be < 1.
#' @param tau Right-censoring time
#' @param keep_latent If TRUE, include latent variables in output
#' @return List with t, delta, C, k
#' @export
rexp_series_md_relaxed_c1 <- function(n, theta, P, tau,
                                       keep_latent = FALSE) {
    m <- length(theta)
    if (any(theta <= 0)) stop("theta must be positive")
    if (nrow(P) != m || ncol(P) != m) stop("P must be m x m")
    if (any(P < 0) || any(P > 1)) stop("P entries must be in [0, 1]")

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

    # Generate candidate sets (relaxed C1: failed component NOT guaranteed)
    C <- matrix(FALSE, nrow = n, ncol = m)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            ki <- k[i]
            for (j in seq_len(m)) {
                C[i, j] <- runif(1) <= P[j, ki]
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
        result$P <- P
    }

    result
}
