# Robustness Interval Computation for Masked Series Systems
#
# Given a fitted C1-C2-C3 model and a reliability estimand of interest,
# compute the maximum violation severity at which the estimate stays within
# a stated tolerance of the baseline. This is the practitioner-facing
# applied tool that connects the sensitivity framework to engineering
# conclusions: "my system-MTTF estimate is robust to C1 violations up to
# alpha = 0.27" (rather than just "the first-order coefficient is X").
#
# Two implementations are provided:
#
# 1. ri_simulation: Monte-Carlo robustness interval. For each candidate
#    alpha, simulates B datasets under the violation, refits the
#    misspecified MLE, computes the estimand, and reports the alpha at
#    which the mean deviation from baseline exceeds the tolerance.
#    Accurate but expensive.
#
# 2. ri_first_order: Closed-form robustness interval from the local
#    sensitivity index ISNI. For exponential components, the ISNI is
#    available analytically (and is zero for the system hazard, giving
#    infinite robustness). For Weibull, the ISNI is computed by
#    perturbation. Fast but only first-order accurate.

# -----------------------------------------------------------------------------
# Monte-Carlo robustness interval
# -----------------------------------------------------------------------------

#' Monte-Carlo robustness interval for a Weibull series-system MLE
#'
#' Computes the robustness interval of a reliability estimand by simulation.
#' For each candidate severity \code{alpha}, the function simulates \code{B}
#' datasets under the specified coarsening-condition violation, refits the
#' standard C1-C2-C3 Weibull series MLE, evaluates the estimand, and records
#' the mean deviation from the baseline value. The robustness interval is the
#' largest contiguous prefix \code{[0, alpha*]} of the severity grid on which
#' the mean deviation stays within \code{tolerance}. This is the
#' practitioner-facing applied tool that connects the sensitivity framework to
#' an engineering conclusion (for example, "the system-MTTF estimate is robust
#' to C1 violations up to alpha = 0.27"). It is accurate but expensive; for a
#' fast first-order alternative see \code{\link{ri_first_order}}.
#'
#' @param shapes_hat Shape MLE estimates (length m).
#' @param scales_hat Scale MLE estimates (length m).
#' @param estimand_fn Function taking \code{(shapes, scales)} and returning a
#'   scalar estimand value. Examples:
#'   \code{function(s, l) hazard_wei_series(t_ref, s, l)} (system hazard);
#'   \code{function(s, l) wei_series_mttf(s, l)} (mean time to failure).
#' @param violation One of \code{"C1"}, \code{"C2"}, or \code{"C3"}.
#' @param tolerance Acceptable absolute deviation from the baseline
#'   estimand value.
#' @param alpha_grid Grid of severity values to evaluate. The robustness
#'   interval is the largest contiguous prefix \code{[0, alpha*]} on which the
#'   mean simulated deviation stays within tolerance.
#' @param n Sample size per simulated dataset.
#' @param B Number of Monte-Carlo replicates per alpha (default 100).
#' @param tau Right-censoring time.
#' @param base_p Baseline off-diagonal masking probability (default 0.5).
#' @param C2_direction For C2 violations, the direction matrix D such that
#'   \code{P(s) = P_0 + s * D}. Default: a fixed cyclic structure.
#' @param verbose If TRUE, print per-alpha progress.
#' @return A list with elements \code{alpha_grid}, \code{mean_deviation},
#'   \code{n_converged}, \code{RI} (the largest alpha within tolerance),
#'   \code{baseline_value}, \code{tolerance}, and \code{violation}.
#' @seealso \code{\link{ri_first_order}} for the fast first-order variant.
#' @export
ri_simulation <- function(shapes_hat, scales_hat, estimand_fn,
                          violation = c("C1", "C2", "C3"),
                          tolerance,
                          alpha_grid,
                          n,
                          B = 100,
                          tau,
                          base_p = 0.5,
                          C2_direction = NULL,
                          verbose = FALSE) {
    violation <- match.arg(violation)
    m <- length(shapes_hat)

    baseline <- estimand_fn(shapes_hat, scales_hat)
    theta_init <- as.vector(rbind(shapes_hat, scales_hat))

    if (violation == "C2" && is.null(C2_direction)) {
        C2_direction <- default_c2_direction(m)
    }

    mean_dev <- numeric(length(alpha_grid))
    n_conv   <- integer(length(alpha_grid))

    for (ai in seq_along(alpha_grid)) {
        alpha <- alpha_grid[ai]
        if (verbose) cat(sprintf("  alpha = %.3f: ", alpha))

        ests <- numeric(B)
        ok   <- 0L
        for (b in seq_len(B)) {
            sim <- simulate_violation(violation, alpha, n,
                                       shapes_hat, scales_hat,
                                       base_p, tau, C2_direction)
            fit <- tryCatch(
                mle_wei_series(sim$t, sim$C, sim$delta, theta0 = theta_init),
                error = function(e) NULL
            )
            if (!is.null(fit) && fit$converged) {
                ok <- ok + 1L
                k_hat <- fit$theta[seq(1, 2 * m, by = 2)]
                l_hat <- fit$theta[seq(2, 2 * m, by = 2)]
                ests[ok] <- estimand_fn(k_hat, l_hat)
            }
        }
        ests <- ests[seq_len(ok)]
        mean_dev[ai] <- if (ok > 0) mean(ests) - baseline else NA_real_
        n_conv[ai]   <- ok
        if (verbose) cat(sprintf("%d/%d converged, mean deviation %+.4f\n",
                                  ok, B, mean_dev[ai]))
    }

    within <- !is.na(mean_dev) & abs(mean_dev) <= tolerance
    if (!any(within)) {
        RI <- NA_real_
    } else if (all(within)) {
        RI <- max(alpha_grid)
    } else {
        first_break <- which(!within)[1]
        RI <- alpha_grid[max(1, first_break - 1L)]
    }

    list(alpha_grid     = alpha_grid,
         mean_deviation = mean_dev,
         n_converged    = n_conv,
         RI             = RI,
         baseline_value = baseline,
         tolerance      = tolerance,
         violation      = violation)
}

# -----------------------------------------------------------------------------
# Closed-form first-order robustness interval (via numerical ISNI)
# -----------------------------------------------------------------------------

#' First-order robustness interval via the Index of Local Sensitivity to
#' Nonignorability (ISNI)
#'
#' Computes a fast, first-order robustness interval for a reliability estimand.
#' The local sensitivity index (ISNI) is the slope of the estimand's deviation
#' from baseline with respect to the violation severity at severity zero, and
#' the robustness interval is \code{tolerance / |ISNI|}. For exponential
#' components the ISNI for the total system hazard is exactly zero under any of
#' the three violations (the all-orders preservation result), giving an
#' infinite robustness interval. For Weibull components the ISNI is estimated
#' numerically by perturbing the masking severity, refitting at a small grid of
#' probe severities, and fitting a linear coefficient. This is the fast
#' counterpart to the Monte-Carlo \code{\link{ri_simulation}}; it is only
#' first-order accurate.
#'
#' @param shapes_hat Shape MLE estimates.
#' @param scales_hat Scale MLE estimates.
#' @param estimand_fn Function \code{(shapes, scales)} returning a scalar
#'   estimand value (for example \code{function(s, l) wei_series_mttf(s, l)}).
#' @param violation One of \code{"C1"}, \code{"C2"}, or \code{"C3"}.
#' @param tolerance Absolute deviation tolerance.
#' @param n Sample size for ISNI estimation.
#' @param B Number of Monte-Carlo replicates per probe (default 50).
#' @param tau Censoring time.
#' @param base_p Off-diagonal base masking probability (default 0.5).
#' @param probe_alpha Two small alpha values used to estimate the slope
#'   of the deviation curve at alpha = 0. Default \code{c(0.05, 0.10)}.
#' @param C2_direction For C2 violations, the direction matrix D such that
#'   \code{P(s) = P_0 + s * D}. Default: a fixed cyclic structure.
#' @return A list with elements \code{ISNI} (the estimated first-order
#'   coefficient), \code{RI} (\code{= tolerance / |ISNI|} if \code{ISNI != 0},
#'   \code{Inf} otherwise), \code{baseline_value}, \code{tolerance},
#'   \code{violation}, \code{probe_alpha}, and \code{probe_deviations}.
#' @seealso \code{\link{ri_simulation}} for the Monte-Carlo variant.
#' @export
ri_first_order <- function(shapes_hat, scales_hat, estimand_fn,
                            violation = c("C1", "C2", "C3"),
                            tolerance,
                            n,
                            B = 50,
                            tau,
                            base_p = 0.5,
                            probe_alpha = c(0.05, 0.10),
                            C2_direction = NULL) {
    violation <- match.arg(violation)
    m <- length(shapes_hat)
    baseline <- estimand_fn(shapes_hat, scales_hat)
    theta_init <- as.vector(rbind(shapes_hat, scales_hat))
    if (violation == "C2" && is.null(C2_direction)) {
        C2_direction <- default_c2_direction(m)
    }

    # Estimate the deviation curve at the probe points
    devs <- numeric(length(probe_alpha))
    for (ai in seq_along(probe_alpha)) {
        alpha <- probe_alpha[ai]
        ests <- numeric(B)
        ok <- 0L
        for (b in seq_len(B)) {
            sim <- simulate_violation(violation, alpha, n,
                                       shapes_hat, scales_hat,
                                       base_p, tau, C2_direction)
            fit <- tryCatch(
                mle_wei_series(sim$t, sim$C, sim$delta, theta0 = theta_init),
                error = function(e) NULL
            )
            if (!is.null(fit) && fit$converged) {
                ok <- ok + 1L
                k_hat <- fit$theta[seq(1, 2 * m, by = 2)]
                l_hat <- fit$theta[seq(2, 2 * m, by = 2)]
                ests[ok] <- estimand_fn(k_hat, l_hat)
            }
        }
        devs[ai] <- if (ok > 0) mean(ests[seq_len(ok)]) - baseline else NA_real_
    }

    # Linear fit to estimate ISNI at alpha = 0
    if (any(is.na(devs))) {
        return(list(ISNI = NA_real_, RI = NA_real_,
                    baseline_value = baseline, tolerance = tolerance,
                    violation = violation))
    }
    isni <- (devs[length(devs)] - devs[1]) /
            (probe_alpha[length(probe_alpha)] - probe_alpha[1])
    RI <- if (isni == 0) Inf else tolerance / abs(isni)

    list(ISNI = isni, RI = RI,
         baseline_value = baseline, tolerance = tolerance,
         violation = violation,
         probe_alpha = probe_alpha, probe_deviations = devs)
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

#' Default direction matrix for C2 perturbation (matches main paper).
default_c2_direction <- function(m) {
    if (m == 5) {
        matrix(c(
             0,    0.3, -0.2,  0.1, -0.2,
            -0.2,  0,    0.3, -0.2,  0.1,
             0.1, -0.2,  0,    0.3, -0.2,
            -0.2,  0.1, -0.2,  0,    0.3,
             0.3, -0.2,  0.1, -0.2,  0
        ), nrow = 5, byrow = TRUE)
    } else {
        D <- matrix(0, nrow = m, ncol = m)
        for (i in seq_len(m)) {
            for (j in seq_len(m)) {
                if (i != j) {
                    D[i, j] <- 0.2 * cos(2 * pi * (i + j) / m)
                }
            }
        }
        D
    }
}

#' Generate one simulated dataset under a violation parameter.
simulate_violation <- function(violation, alpha, n,
                                shapes, scales, base_p, tau, C2_direction) {
    m <- length(shapes)
    if (violation == "C2") {
        P <- make_P_matrix(m, "uniform", p = base_p) + alpha * C2_direction
        offdiag <- row(P) != col(P)
        P[offdiag] <- pmax(0.05, pmin(0.95, P[offdiag]))
        diag(P) <- 1
        return(rwei_series_md_c1_c3(n = n, shapes = shapes, scales = scales,
                                      P = P, tau = tau))
    }
    if (violation == "C3") {
        return(rwei_series_md_c1_c2(n = n, shapes = shapes, scales = scales,
                                      alpha = alpha, base_p = base_p,
                                      tau = tau))
    }
    if (violation == "C1") {
        Tm <- matrix(nrow = n, ncol = m)
        for (j in seq_len(m)) {
            Tm[, j] <- rweibull(n, shape = shapes[j], scale = scales[j])
        }
        sys_time <- apply(Tm, 1, min)
        k_true   <- apply(Tm, 1, which.min)
        tau_v    <- rep(tau, length.out = n)
        delta    <- as.numeric(sys_time <= tau_v)
        obs_time <- pmin(sys_time, tau_v)

        diag_p <- 1 - alpha
        P <- matrix(base_p, nrow = m, ncol = m)
        diag(P) <- diag_p
        C <- matrix(FALSE, nrow = n, ncol = m)
        for (i in seq_len(n)) {
            if (delta[i] == 1) {
                ki <- k_true[i]
                for (j in seq_len(m)) {
                    C[i, j] <- runif(1) <= P[j, ki]
                }
            }
        }
        has_set <- rowSums(C) >= 1
        keep    <- (delta == 0) | has_set
        return(list(t = obs_time[keep], delta = delta[keep],
                     C = C[keep, , drop = FALSE], k = k_true[keep]))
    }
    stop("unknown violation: ", violation)
}
