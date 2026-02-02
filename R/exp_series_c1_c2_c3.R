# =============================================================================
# Exponential Series System with Masked Data (C1, C2, C3)
# Minimal, dependency-free implementation for simulations
# =============================================================================

# -----------------------------------------------------------------------------
# Helper Functions (no external dependencies)
# -----------------------------------------------------------------------------

#' Decode matrix from prefixed columns
#'
#' Extracts columns with a common prefix (e.g., x1, x2, x3) into a matrix.
#'
#' @param df Data frame
#' @param prefix Column prefix (e.g., "x" matches x1, x2, ...)
#' @return Matrix with columns in numeric order, or NULL if no matches
#' @export
decode_matrix <- function(df, prefix) {
    pattern <- paste0("^", prefix, "([0-9]+)$")
    matches <- grep(pattern, names(df), value = TRUE)
    if (length(matches) == 0) return(NULL)

    # Extract numeric suffixes and sort
    nums <- as.integer(sub(pattern, "\\1", matches))
    ord <- order(nums)
    matches <- matches[ord]

    as.matrix(df[, matches, drop = FALSE])
}

#' Encode matrix to prefixed columns
#'
#' Converts a matrix to named columns with a prefix (e.g., x1, x2, x3).
#'
#' @param mat Matrix
#' @param prefix Column prefix
#' @return Data frame with prefixed column names
#' @export
encode_matrix <- function(mat, prefix) {
    df <- as.data.frame(mat)
    names(df) <- paste0(prefix, seq_len(ncol(mat)))
    df
}

# -----------------------------------------------------------------------------
# Log-Likelihood, Score, and Fisher Information
# -----------------------------------------------------------------------------

#' Log-likelihood for exponential series system (C1, C2, C3)
#'
#' Returns a log-likelihood function for masked data from an exponential series
#' system under conditions C1, C2, C3.
#'
#' @param t Numeric vector of system lifetimes
#' @param C Logical matrix of candidate sets (n x m), TRUE if component in set
#' @param delta Numeric/logical vector of censoring indicators (1=observed, 0=censored).
#'        If NULL, assumes all observations are uncensored.
#' @return Function that takes theta (rate parameters) and returns log-likelihood
#' @export
#' @examples
#' t <- c(0.5, 1.0, 0.8)
#' C <- matrix(c(TRUE, TRUE, FALSE,
#'               TRUE, FALSE, TRUE,
#'               FALSE, TRUE, TRUE), nrow = 3, byrow = TRUE)
#' ll <- loglik_exp_series(t, C)
#' ll(c(1, 1.5, 2))
loglik_exp_series <- function(t, C, delta = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    # Pre-compute sum of times
    sum_t <- sum(t)

    function(theta) {
        if (length(theta) != m) stop("theta must have length ", m)
        if (any(theta <= 0)) return(-Inf)

        # Survival contribution: -sum(t) * sum(theta)
        ll <- -sum_t * sum(theta)

        # Hazard contribution for uncensored observations
        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                theta_C <- sum(theta[C[i, ]])
                if (theta_C <= 0) return(-Inf)
                ll <- ll + log(theta_C)
            }
        }
        ll
    }
}

#' Score (gradient) for exponential series system (C1, C2, C3)
#'
#' Returns the score function (gradient of log-likelihood) for masked data.
#'
#' @inheritParams loglik_exp_series
#' @return Function that takes theta and returns gradient vector
#' @export
score_exp_series <- function(t, C, delta = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    sum_t <- sum(t)

    function(theta) {
        if (length(theta) != m) stop("theta must have length ", m)
        if (any(theta <= 0)) return(rep(NA_real_, m))

        # Survival contribution: -sum(t) for each component
        g <- rep(-sum_t, m)

        # Hazard contribution
        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                C_i <- C[i, ]
                theta_C <- sum(theta[C_i])
                if (theta_C <= 0) return(rep(NA_real_, m))
                # Add 1/theta_C for each j in C_i
                g[C_i] <- g[C_i] + 1 / theta_C
            }
        }
        g
    }
}

#' Fisher Information Matrix for exponential series system (C1, C2, C3)
#'
#' Returns the observed Fisher information matrix (negative Hessian).
#'
#' @inheritParams loglik_exp_series
#' @return Function that takes theta and returns m x m FIM
#' @export
fim_exp_series <- function(t, C, delta = NULL) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)
    delta <- as.numeric(delta)

    function(theta) {
        if (length(theta) != m) stop("theta must have length ", m)
        if (any(theta <= 0)) return(matrix(NA_real_, m, m))

        I <- matrix(0, m, m)

        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                C_i <- C[i, ]
                theta_C <- sum(theta[C_i])
                if (theta_C <= 0) return(matrix(NA_real_, m, m))

                # FIM contribution: outer product scaled by 1/theta_C^2
                idx <- which(C_i)
                for (j in idx) {
                    for (k in idx) {
                        I[j, k] <- I[j, k] + 1 / theta_C^2
                    }
                }
            }
        }
        I
    }
}

# -----------------------------------------------------------------------------
# Maximum Likelihood Estimation
# -----------------------------------------------------------------------------

#' MLE for exponential series system (C1, C2, C3)
#'
#' Computes maximum likelihood estimates of component failure rates.
#'
#' @inheritParams loglik_exp_series
#' @param theta0 Initial parameter values. If NULL, uses method of moments.
#' @param lower Lower bound for parameters (default 1e-10)
#' @param upper Upper bound for parameters (default Inf)
#' @return List with components:
#'   - `theta`: MLE of rate parameters
#'   - `se`: Standard errors (from FIM)
#'   - `vcov`: Variance-covariance matrix
#'   - `loglik`: Log-likelihood at MLE
#'   - `converged`: Logical, did optimization converge?
#'   - `fim`: Fisher information matrix at MLE
#' @export
#' @examples
#' \dontrun{
#' # Generate data
#' set.seed(42)
#' sim <- rexp_series_md(n = 200, theta = c(1, 1.5, 2), p = 0.3, tau = 2)
#' fit <- mle_exp_series(sim$t, sim$C, sim$delta)
#' print(fit$theta)
#' }
mle_exp_series <- function(t, C, delta = NULL, theta0 = NULL,
                            lower = 1e-10, upper = Inf) {
    n <- length(t)
    m <- ncol(C)
    if (is.null(delta)) delta <- rep(1, n)

    # Get likelihood and score functions
    ll <- loglik_exp_series(t, C, delta)
    sc <- score_exp_series(t, C, delta)
    fim_fn <- fim_exp_series(t, C, delta)

    # Initialize theta0 if not provided (method of moments)
    if (is.null(theta0)) {
        n_uncensored <- sum(delta == 1)
        total_hazard <- n_uncensored / sum(t)
        theta0 <- rep(total_hazard / m, m)
    }
    theta0 <- pmax(theta0, lower)

    # Optimize using L-BFGS-B
    result <- optim(
        par = theta0,
        fn = ll,
        gr = sc,
        method = "L-BFGS-B",
        lower = rep(lower, m),
        upper = rep(upper, m),
        control = list(fnscale = -1, maxit = 1000)
    )

    # Compute FIM and variance-covariance
    fim <- fim_fn(result$par)
    vcov <- tryCatch(
        solve(fim),
        error = function(e) matrix(NA_real_, m, m)
    )
    se <- sqrt(diag(vcov))

    list(
        theta = result$par,
        se = se,
        vcov = vcov,
        loglik = result$value,
        converged = (result$convergence == 0),
        fim = fim,
        n = n,
        m = m
    )
}

# -----------------------------------------------------------------------------
# Data Generation
# -----------------------------------------------------------------------------

#' Generate masked data from exponential series system
#'
#' Simulates component lifetimes, applies series system structure, right-censoring,
#' and generates candidate sets using Bernoulli model (C1, C2, C3).
#'
#' @param n Sample size
#' @param theta Rate parameters for exponential component lifetimes
#' @param p Probability each non-failed component is in candidate set
#' @param tau Right-censoring time (scalar or vector of length n)
#' @param keep_latent If TRUE, include latent component times in output
#' @return List with components:
#'   - `t`: System lifetimes (possibly censored)
#'   - `delta`: Censoring indicator (1 = observed failure, 0 = censored)
#'   - `C`: Candidate set matrix (n x m logical)
#'   - `k`: True failed component (latent)
#'   - `Tm`: Component lifetimes matrix (if keep_latent = TRUE)
#' @export
#' @examples
#' set.seed(42)
#' sim <- rexp_series_md(n = 100, theta = c(1, 1.5, 2), p = 0.3, tau = 2)
#' str(sim)
rexp_series_md <- function(n, theta, p, tau, keep_latent = FALSE) {
    m <- length(theta)
    if (any(theta <= 0)) stop("theta must be positive")
    if (p < 0 || p > 1) stop("p must be in [0, 1]")

    # Generate component lifetimes
    Tm <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
        Tm[, j] <- rexp(n, rate = theta[j])
    }

    # System lifetime = min of component lifetimes
    sys_time <- apply(Tm, 1, min)
    k <- apply(Tm, 1, which.min)  # Failed component

    # Right censoring
    tau <- rep(tau, length.out = n)
    delta <- as.numeric(sys_time <= tau)
    obs_time <- pmin(sys_time, tau)

    # Generate candidate sets (Bernoulli model, C1-C2-C3)
    # Failed component always in set (prob = 1), others with prob = p
    C <- matrix(runif(n * m) <= p, nrow = n, ncol = m)
    C[cbind(seq_len(n), k)] <- TRUE  # Ensure failed component in set

    # For censored observations, no candidate set (all FALSE)
    C[delta == 0, ] <- FALSE

    result <- list(
        t = obs_time,
        delta = delta,
        C = C,
        k = k
    )

    if (keep_latent) {
        result$Tm <- Tm
    }

    result
}

# -----------------------------------------------------------------------------
# Data Frame Interface (convenience wrappers)
# -----------------------------------------------------------------------------

#' Wrapper: Log-likelihood from data frame
#'
#' @param df Data frame with columns t, delta, and x1, x2, ..., xm
#' @param tvar Column name for system lifetime (default "t")
#' @param deltavar Column name for censoring indicator (default "delta")
#' @param setvar Column prefix for candidate set (default "x")
#' @return Log-likelihood function
#' @export
loglik_exp_series_df <- function(df, tvar = "t", deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    if (is.null(C)) stop("No candidate set columns found with prefix '", setvar, "'")
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    loglik_exp_series(t, C, delta)
}

#' Wrapper: Score from data frame
#' @inheritParams loglik_exp_series_df
#' @return Score function
#' @export
score_exp_series_df <- function(df, tvar = "t", deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    score_exp_series(t, C, delta)
}

#' Wrapper: FIM from data frame
#' @inheritParams loglik_exp_series_df
#' @return FIM function
#' @export
fim_exp_series_df <- function(df, tvar = "t", deltavar = "delta", setvar = "x") {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    fim_exp_series(t, C, delta)
}

#' Wrapper: MLE from data frame
#' @inheritParams loglik_exp_series_df
#' @inheritParams mle_exp_series
#' @return MLE result list
#' @export
mle_exp_series_df <- function(df, tvar = "t", deltavar = "delta", setvar = "x",
                               theta0 = NULL, lower = 1e-10, upper = Inf) {
    t <- df[[tvar]]
    C <- decode_matrix(df, setvar)
    delta <- if (deltavar %in% names(df)) df[[deltavar]] else NULL
    mle_exp_series(t, C, delta, theta0, lower, upper)
}

#' Convert simulation output to data frame
#'
#' @param sim Output from rexp_series_md()
#' @return Data frame with t, delta, x1, x2, ..., xm columns
#' @export
as_dataframe <- function(sim) {
    df <- data.frame(t = sim$t, delta = sim$delta)
    C_df <- encode_matrix(sim$C, "x")
    cbind(df, C_df)
}
