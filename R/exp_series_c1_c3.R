# =============================================================================
# Exponential Series System with Masked Data (C1, C3) - Relaxed C2
# General Bernoulli model where p_j(k) can vary with k
# =============================================================================

#' Compute pi_k(c) = P(C = c | K = k) for the general Bernoulli model
#'
#' @param c Logical vector indicating which components are in candidate set
#' @param k Index of failed component
#' @param P Inclusion probability matrix: `P[j,k] = P(j in C | K = k)`
#' @return Probability of observing candidate set c given component k failed
#' @export
compute_pi <- function(c, k, P) {
    m <- length(c)
    if (!c[k]) return(0)  # C1 violation

    prob <- 1
    for (j in seq_len(m)) {
        if (j == k) next  # p_k(k) = 1 by C1
        if (c[j]) {
            prob <- prob * P[j, k]
        } else {
            prob <- prob * (1 - P[j, k])
        }
    }
    prob
}

#' Compute pi_k(c) for all k in c
#'
#' @param c Logical vector indicating candidate set
#' @param P Inclusion probability matrix
#' @return Named vector of pi_k values for k in c
#' @export
compute_pi_all <- function(c, P) {
    candidates <- which(c)
    pi_vals <- sapply(candidates, function(k) compute_pi(c, k, P))
    names(pi_vals) <- candidates
    pi_vals
}

# -----------------------------------------------------------------------------
# Log-Likelihood, Score, and Fisher Information for Relaxed C2 Model
# -----------------------------------------------------------------------------

#' Log-likelihood for exponential series system (C1, C3) - Relaxed C2
#'
#' Returns a log-likelihood function for masked data where the masking
#' probabilities p_j(k) can depend on which component k failed.
#'
#' @param t Numeric vector of system lifetimes
#' @param C Logical matrix of candidate sets (n x m)
#' @param delta Censoring indicators (1=observed, 0=censored)
#' @param P Inclusion probability matrix (m x m): `P[j,k] = P(j in C | K = k)`.
#'        Diagonal must be 1 (C1). If NULL, treated as unknown parameters.
#' @return Function that takes (theta, P) and returns log-likelihood.
#'         If P was provided, function takes only theta.
#' @export
loglik_exp_series_c1_c3 <- function(t, C, delta = NULL, P = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    sum_t <- sum(t)
    P_fixed <- !is.null(P)

    if (P_fixed) {
        # P is known - precompute pi values
        pi_list <- vector("list", n)
        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                pi_list[[i]] <- compute_pi_all(C[i, ], P)
            }
        }

        function(theta) {
            if (length(theta) != m) stop("theta must have length ", m)
            if (any(theta <= 0)) return(-Inf)

            ll <- -sum_t * sum(theta)

            for (i in seq_len(n)) {
                if (delta[i] == 1) {
                    candidates <- which(C[i, ])
                    pi_vals <- pi_list[[i]]

                    # Weighted hazard sum
                    hazard_sum <- sum(theta[candidates] * pi_vals)
                    if (hazard_sum <= 0) return(-Inf)
                    ll <- ll + log(hazard_sum)
                }
            }
            ll
        }
    } else {
        # P is unknown - will be passed as parameter
        function(theta, P) {
            if (length(theta) != m) stop("theta must have length ", m)
            if (any(theta <= 0)) return(-Inf)
            if (!all(diag(P) == 1)) stop("Diagonal of P must be 1 (C1)")
            if (any(P < 0) || any(P > 1)) return(-Inf)

            ll <- -sum_t * sum(theta)

            for (i in seq_len(n)) {
                if (delta[i] == 1) {
                    candidates <- which(C[i, ])
                    pi_vals <- sapply(candidates, function(k) compute_pi(C[i,], k, P))

                    hazard_sum <- sum(theta[candidates] * pi_vals)
                    if (hazard_sum <= 0) return(-Inf)
                    ll <- ll + log(hazard_sum)
                }
            }
            ll
        }
    }
}

#' Score (gradient) for exponential series system (C1, C3) - Relaxed C2
#'
#' Returns the score function (gradient of log-likelihood) for masked data
#' with informative masking (relaxed C2). The P matrix must be known.
#'
#' @inheritParams loglik_exp_series_c1_c3
#' @return Function that takes theta and returns gradient vector
#' @export
score_exp_series_c1_c3 <- function(t, C, delta = NULL, P) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    sum_t <- sum(t)

    # Precompute pi values for each uncensored observation
    pi_list <- vector("list", n)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            pi_list[[i]] <- compute_pi_all(C[i, ], P)
        }
    }

    function(theta) {
        if (length(theta) != m) stop("theta must have length ", m)
        if (any(theta <= 0)) return(rep(NA_real_, m))

        g <- rep(-sum_t, m)

        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                candidates <- which(C[i, ])
                pi_vals <- pi_list[[i]]

                # Weighted hazard sum: sum_{k in c_i} theta_k * pi_k(c_i)
                denom <- sum(theta[candidates] * pi_vals)
                if (denom <= 0) return(rep(NA_real_, m))

                # For each j in c_i: add pi_j(c_i) / denom
                for (idx in seq_along(candidates)) {
                    j <- candidates[idx]
                    g[j] <- g[j] + pi_vals[idx] / denom
                }
            }
        }
        g
    }
}

#' Fisher Information Matrix for exponential series system (C1, C3) - Relaxed C2
#'
#' Returns the observed Fisher information matrix (negative Hessian) for
#' masked data with informative masking. The P matrix must be known.
#'
#' @inheritParams loglik_exp_series_c1_c3
#' @return Function that takes theta and returns m x m FIM
#' @export
fim_exp_series_c1_c3 <- function(t, C, delta = NULL, P) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    # Precompute pi values
    pi_list <- vector("list", n)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            pi_list[[i]] <- compute_pi_all(C[i, ], P)
        }
    }

    function(theta) {
        if (length(theta) != m) stop("theta must have length ", m)
        if (any(theta <= 0)) return(matrix(NA_real_, m, m))

        I <- matrix(0, m, m)

        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                candidates <- which(C[i, ])
                pi_vals <- pi_list[[i]]

                # Weighted hazard sum
                denom <- sum(theta[candidates] * pi_vals)
                if (denom <= 0) return(matrix(NA_real_, m, m))

                inv_sq <- 1 / denom^2

                # FIM contribution: pi_j * pi_k / denom^2 for j,k in c_i
                for (a in seq_along(candidates)) {
                    j <- candidates[a]
                    for (b in seq_along(candidates)) {
                        k <- candidates[b]
                        I[j, k] <- I[j, k] + pi_vals[a] * pi_vals[b] * inv_sq
                    }
                }
            }
        }
        I
    }
}

#' MLE for exponential series (C1, C3) with unknown P matrix
#'
#' Jointly estimates rate parameters theta and inclusion probability matrix P.
#'
#' @param t Numeric vector of system lifetimes
#' @param C Logical matrix of candidate sets (n x m)
#' @param delta Censoring indicators
#' @param theta0 Initial rate parameters (if NULL, uses method of moments)
#' @param P0 Initial P matrix (if NULL, uses uniform p = 0.5)
#' @param fixed_P If not NULL, a fixed P matrix (only estimate theta)
#' @return List with theta, P, loglik, converged, se, vcov
#' @export
mle_exp_series_c1_c3 <- function(t, C, delta = NULL, theta0 = NULL, P0 = NULL,
                                  fixed_P = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)

    # If P is fixed, use analytical score and FIM
    if (!is.null(fixed_P)) {
        ll_fn <- loglik_exp_series_c1_c3(t, C, delta, P = fixed_P)
        sc_fn <- score_exp_series_c1_c3(t, C, delta, P = fixed_P)
        fim_fn <- fim_exp_series_c1_c3(t, C, delta, P = fixed_P)

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
            control = list(fnscale = -1)
        )

        # Compute analytical FIM at MLE
        fim <- fim_fn(result$par)
        vcov <- tryCatch(solve(fim), error = function(e) matrix(NA, m, m))
        se <- sqrt(diag(vcov))

        return(list(
            theta = result$par,
            P = fixed_P,
            se = se,
            vcov = vcov,
            loglik = result$value,
            converged = result$convergence == 0,
            fim = fim
        ))
    }

    # Joint optimization of theta and P
    ll_fn <- loglik_exp_series_c1_c3(t, C, delta, P = NULL)

    # Initialize
    if (is.null(theta0)) {
        n_unc <- sum(delta == 1)
        theta0 <- rep(n_unc / (sum(t) * m), m)
    }
    if (is.null(P0)) {
        P0 <- matrix(0.5, m, m)
        diag(P0) <- 1
    }

    # Pack parameters: theta (m), then off-diagonal of P (m*(m-1))
    pack_params <- function(theta, P) {
        P_offdiag <- P[row(P) != col(P)]  # Off-diagonal elements
        c(theta, P_offdiag)
    }

    unpack_params <- function(par) {
        theta <- par[1:m]
        P_offdiag <- par[(m+1):length(par)]
        P <- matrix(0, m, m)
        P[row(P) != col(P)] <- P_offdiag
        diag(P) <- 1
        list(theta = theta, P = P)
    }

    par0 <- pack_params(theta0, P0)
    n_par <- length(par0)

    # Objective function
    obj <- function(par) {
        params <- unpack_params(par)
        ll_fn(params$theta, params$P)
    }

    # Bounds
    lower <- c(rep(1e-10, m), rep(0, m * (m - 1)))
    upper <- c(rep(Inf, m), rep(1, m * (m - 1)))

    result <- optim(
        par = par0,
        fn = obj,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(fnscale = -1, maxit = 1000)
    )

    params <- unpack_params(result$par)

    # Compute numerical Hessian for theta only (for standard errors)
    h <- 1e-5
    H_theta <- matrix(0, m, m)
    for (j in seq_len(m)) {
        for (k in seq_len(m)) {
            theta_pp <- theta_pm <- theta_mp <- theta_mm <- params$theta
            theta_pp[j] <- theta_pp[j] + h
            theta_pp[k] <- theta_pp[k] + h
            theta_pm[j] <- theta_pm[j] + h
            theta_pm[k] <- theta_pm[k] - h
            theta_mp[j] <- theta_mp[j] - h
            theta_mp[k] <- theta_mp[k] + h
            theta_mm[j] <- theta_mm[j] - h
            theta_mm[k] <- theta_mm[k] - h
            H_theta[j, k] <- (ll_fn(theta_pp, params$P) - ll_fn(theta_pm, params$P) -
                              ll_fn(theta_mp, params$P) + ll_fn(theta_mm, params$P)) / (4 * h^2)
        }
    }
    fim_theta <- -H_theta
    vcov_theta <- tryCatch(solve(fim_theta), error = function(e) matrix(NA, m, m))
    se_theta <- sqrt(diag(vcov_theta))

    list(
        theta = params$theta,
        P = params$P,
        se = se_theta,
        vcov = vcov_theta,
        loglik = result$value,
        converged = result$convergence == 0,
        fim = fim_theta,
        n_params = n_par
    )
}

# -----------------------------------------------------------------------------
# Data Generation with Relaxed C2
# -----------------------------------------------------------------------------

#' Generate masked data with general Bernoulli model (relaxed C2)
#'
#' @param n Sample size
#' @param theta Rate parameters for exponential components
#' @param P Inclusion probability matrix (m x m): `P[j,k] = P(j in C | K = k)`
#' @param tau Right-censoring time
#' @param keep_latent If TRUE, include latent variables in output
#' @return List with t, delta, C, k, and optionally Tm, P
#' @export
rexp_series_md_c1_c3 <- function(n, theta, P, tau, keep_latent = FALSE) {
    m <- length(theta)
    if (any(theta <= 0)) stop("theta must be positive")
    if (nrow(P) != m || ncol(P) != m) stop("P must be m x m")
    if (!all(diag(P) == 1)) stop("Diagonal of P must be 1 (C1)")

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

    # Generate candidate sets using general Bernoulli model
    C <- matrix(FALSE, nrow = n, ncol = m)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            ki <- k[i]
            for (j in seq_len(m)) {
                if (j == ki) {
                    C[i, j] <- TRUE  # C1: failed component always in
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
        k = k
    )

    if (keep_latent) {
        result$Tm <- Tm
        result$P <- P
    }

    result
}

#' Create a P matrix for the general Bernoulli model
#'
#' Helper function to create inclusion probability matrices.
#'
#' @param m Number of components
#' @param type Type of matrix: "uniform" (all off-diag = p),
#'        "symmetric" (`P[j,k] = P[k,j]`), or "full" (all different)
#' @param p Baseline probability (for uniform type)
#' @param values For "full" type: vector of m*(m-1) off-diagonal values
#'        in column-major order
#' @return m x m inclusion probability matrix with diagonal = 1
#' @export
#' @examples
#' # Uniform (satisfies C2)
#' P <- make_P_matrix(3, "uniform", p = 0.3)
#'
#' # Full specification
#' P <- make_P_matrix(3, "full", values = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
make_P_matrix <- function(m, type = "uniform", p = 0.5, values = NULL) {
    P <- matrix(0, m, m)
    diag(P) <- 1

    if (type == "uniform") {
        P[row(P) != col(P)] <- p
    } else if (type == "symmetric") {
        if (is.null(values)) {
            values <- rep(p, m * (m - 1) / 2)
        }
        idx <- 1
        for (j in 1:(m-1)) {
            for (k in (j+1):m) {
                P[j, k] <- values[idx]
                P[k, j] <- values[idx]
                idx <- idx + 1
            }
        }
    } else if (type == "full") {
        if (is.null(values)) stop("values required for full type")
        if (length(values) != m * (m - 1)) {
            stop("values must have length m*(m-1) = ", m * (m - 1))
        }
        P[row(P) != col(P)] <- values
    }

    P
}

#' Check if P matrix satisfies C2
#'
#' @param P Inclusion probability matrix
#' @param tol Tolerance for equality check
#' @return TRUE if C2 is satisfied (rows have constant off-diagonal)
#' @export
satisfies_C2 <- function(P, tol = 1e-10) {
    m <- nrow(P)
    for (j in seq_len(m)) {
        offdiag <- P[j, -j]
        if (max(offdiag) - min(offdiag) > tol) {
            return(FALSE)
        }
    }
    TRUE
}

# -----------------------------------------------------------------------------
# Data Frame Interface for Relaxed C2 Model
# -----------------------------------------------------------------------------

#' Wrapper: Log-likelihood (relaxed C2) from data frame
#'
#' @param df Data frame with columns t, delta, and x1, x2, ..., xm
#' @param P Inclusion probability matrix (m x m). If NULL, treated as unknown.
#' @param tvar Column name for system lifetime (default "t")
#' @param deltavar Column name for censoring indicator (default "delta")
#' @param setvar Column prefix for candidate set (default "x")
#' @return Log-likelihood function
#' @export
loglik_exp_series_c1_c3_df <- function(df, P = NULL, tvar = "t",
                                        deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    loglik_exp_series_c1_c3(t, C, delta, P)
}

#' Wrapper: Score (relaxed C2) from data frame
#' @inheritParams loglik_exp_series_c1_c3_df
#' @return Score function
#' @export
score_exp_series_c1_c3_df <- function(df, P, tvar = "t",
                                       deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    score_exp_series_c1_c3(t, C, delta, P)
}

#' Wrapper: FIM (relaxed C2) from data frame
#' @inheritParams loglik_exp_series_c1_c3_df
#' @return FIM function
#' @export
fim_exp_series_c1_c3_df <- function(df, P, tvar = "t",
                                     deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    fim_exp_series_c1_c3(t, C, delta, P)
}

#' Wrapper: MLE (relaxed C2) from data frame
#' @inheritParams loglik_exp_series_c1_c3_df
#' @param theta0 Initial rate parameters
#' @param P0 Initial P matrix (for joint estimation)
#' @param fixed_P If not NULL, a fixed P matrix (only estimate theta)
#' @return MLE result list
#' @export
mle_exp_series_c1_c3_df <- function(df, tvar = "t", deltavar = "delta",
                                     setvar = "x", theta0 = NULL,
                                     P0 = NULL, fixed_P = NULL) {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    mle_exp_series_c1_c3(t, C, delta, theta0, P0, fixed_P)
}
