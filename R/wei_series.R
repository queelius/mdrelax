# =============================================================================
# Weibull Series System Distribution Functions
# Dependency-free implementation following R conventions (d/p/q/r)
# =============================================================================

# -----------------------------------------------------------------------------
# Distribution Functions
# -----------------------------------------------------------------------------

#' PDF of Weibull series system
#'
#' The system fails when the first component fails. The PDF is:
#' f(t) = h_sys(t) * R_sys(t) where h_sys = sum of component hazards
#' and R_sys = product of component reliabilities.
#'
#' @param t Time points (non-negative)
#' @param shapes Shape parameters (k_1, ..., k_m)
#' @param scales Scale parameters (beta_1, ..., beta_m)
#' @param log If TRUE, return log-density
#' @return Density values
#' @export
dwei_series <- function(t, shapes, scales, log = FALSE) {
    m <- length(shapes)
    if (length(scales) != m) stop("shapes and scales must have same length")

    # Compute log-density for numerical stability
    # log f(t) = log(h_sys(t)) + log(R_sys(t))
    # log(R_sys(t)) = -sum_j (t/beta_j)^k_j
    # h_sys(t) = sum_j (k_j/beta_j) * (t/beta_j)^(k_j - 1)

    log_surv <- rep(0, length(t))
    hazard <- rep(0, length(t))

    for (j in seq_len(m)) {
        z <- t / scales[j]
        log_surv <- log_surv - z^shapes[j]
        hazard <- hazard + (shapes[j] / scales[j]) * z^(shapes[j] - 1)
    }

    ld <- ifelse(t > 0, log(hazard) + log_surv, -Inf)
    ld[t == 0 & all(shapes >= 1)] <- -Inf  # Handle t=0

    if (log) ld else exp(ld)
}

#' CDF of Weibull series system
#'
#' @param q Quantile values
#' @param shapes Shape parameters
#' @param scales Scale parameters
#' @param lower.tail If TRUE (default), P(T <= q)
#' @param log.p If TRUE, return log probability
#' @return Probabilities
#' @export
pwei_series <- function(q, shapes, scales, lower.tail = TRUE, log.p = FALSE) {
    m <- length(shapes)
    if (length(scales) != m) stop("shapes and scales must have same length")

    # R_sys(t) = exp(-sum_j (t/beta_j)^k_j)
    log_surv <- rep(0, length(q))
    for (j in seq_len(m)) {
        log_surv <- log_surv - (pmax(q, 0) / scales[j])^shapes[j]
    }

    if (lower.tail) {
        # P(T <= q) = 1 - R(q)
        p <- -expm1(log_surv)  # More numerically stable than 1 - exp(log_surv)
        p[q <= 0] <- 0
    } else {
        p <- exp(log_surv)
        p[q <= 0] <- 1
    }

    if (log.p) log(p) else p
}

#' Quantile function of Weibull series system
#'
#' Uses Newton's method to invert the CDF.
#'
#' @param p Probabilities
#' @param shapes Shape parameters
#' @param scales Scale parameters
#' @param lower.tail If TRUE (default), find q such that P(T <= q) = p
#' @param tol Convergence tolerance
#' @param max_iter Maximum iterations
#' @return Quantile values
#' @export
qwei_series <- function(p, shapes, scales, lower.tail = TRUE,
                         tol = 1e-10, max_iter = 100) {
    if (!lower.tail) p <- 1 - p

    sapply(p, function(pp) {
        if (pp <= 0) return(0)
        if (pp >= 1) return(Inf)

        # Initial guess: use mean of individual Weibull quantiles
        q <- min(sapply(seq_along(shapes), function(j) {
            scales[j] * (-log(1 - pp))^(1 / shapes[j])
        }))

        for (iter in seq_len(max_iter)) {
            F_q <- pwei_series(q, shapes, scales)
            f_q <- dwei_series(q, shapes, scales)
            if (f_q <= 0) break

            delta <- (F_q - pp) / f_q
            q <- max(q - delta, q / 2)  # Ensure q stays positive

            if (abs(delta) < tol * q) break
        }
        q
    })
}

#' Random generation from Weibull series system
#'
#' @param n Number of observations
#' @param shapes Shape parameters
#' @param scales Scale parameters
#' @return Vector of system lifetimes
#' @export
rwei_series <- function(n, shapes, scales) {
    m <- length(shapes)
    if (length(scales) != m) stop("shapes and scales must have same length")

    # Generate component lifetimes and take minimum
    Tm <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
        Tm[, j] <- rweibull(n, shape = shapes[j], scale = scales[j])
    }
    apply(Tm, 1, min)
}

# -----------------------------------------------------------------------------
# Hazard and Survival Functions
# -----------------------------------------------------------------------------

#' Hazard function of Weibull series system
#'
#' @param t Time points
#' @param shapes Shape parameters
#' @param scales Scale parameters
#' @return Hazard values
#' @export
hazard_wei_series <- function(t, shapes, scales) {
    m <- length(shapes)
    h <- rep(0, length(t))
    for (j in seq_len(m)) {
        z <- t / scales[j]
        h <- h + (shapes[j] / scales[j]) * z^(shapes[j] - 1)
    }
    h
}

#' Survival function of Weibull series system
#'
#' @param t Time points
#' @param shapes Shape parameters
#' @param scales Scale parameters
#' @param log.p If TRUE, return log survival probability
#' @return Survival probabilities
#' @export
surv_wei_series <- function(t, shapes, scales, log.p = FALSE) {
    m <- length(shapes)
    log_surv <- rep(0, length(t))
    for (j in seq_len(m)) {
        log_surv <- log_surv - (pmax(t, 0) / scales[j])^shapes[j]
    }
    if (log.p) log_surv else exp(log_surv)
}

# -----------------------------------------------------------------------------
# System Statistics
# -----------------------------------------------------------------------------

#' Mean time to failure (MTTF) of Weibull series system
#'
#' Computed by numerical integration of the survival function.
#'
#' @param shapes Shape parameters
#' @param scales Scale parameters
#' @param upper Upper integration limit (default: large value)
#' @return MTTF
#' @export
wei_series_mttf <- function(shapes, scales, upper = NULL) {
    if (is.null(upper)) {
        # Choose reasonable upper bound
        upper <- max(scales * gamma(1 + 1/shapes)) * 10
    }
    integrate(function(t) surv_wei_series(t, shapes, scales),
              lower = 0, upper = upper)$value
}
