#!/usr/bin/env Rscript
# =============================================================================
# Sensitivity Sweep: How does C1 violation severity affect C1-C2-C3 MLE?
# =============================================================================
#
# Generates data under C1 violation (true cause may not be in candidate set)
# with varying severity alpha, fits the misspecified C1-C2-C3 model, and
# measures bias, RMSE, and coverage as a function of alpha.
#
# Theoretical question tested: For Weibull components, does the total system
# hazard remain robust under C1 violation? The exponential case is proved
# exactly preserved in sensitivity_framework.tex (Cor. c1-hazard-general)
# via the profile-likelihood factorization that does not apply to Weibull.
#
# C1 violation parameterization: p_k(k) = 1 - alpha_C1 (probability of
# including the failed component in the candidate set); p_j(k) = p_0 for
# j != k (off-diagonal masking unaffected). Observations with empty
# candidate sets are dropped to match the standard analyst workflow.
#
# System: 5-component Weibull series (matches the C2 and C3 sweeps)
#   shapes = (2.0, 1.5, 1.2, 1.8, 1.0)
#   scales = (3.0, 4.0, 5.0, 3.5, 4.5)
#   tau = 5, n = 500, B = 200
#
# Usage:
#   Rscript paper/run_sensitivity_sweep_c1.R [--quick]
#
# Outputs:
#   paper/data/sensitivity_sweep_c1.rds
#   inst/simulations/figures/fig_c1_bias_vs_alpha.pdf
#   inst/simulations/figures/fig_c1_rmse_vs_alpha.pdf
#   inst/simulations/figures/fig_c1_coverage_vs_alpha.pdf
#   inst/simulations/figures/fig_c1_median_bias_vs_alpha.pdf
#   inst/simulations/figures/fig_c1_mad_vs_alpha.pdf
#
# =============================================================================

devtools::load_all(quiet = TRUE)

set.seed(42)

# -----------------------------------------------------------------------------
# Parse Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args
n_arg <- grep("^--n=", args, value = TRUE)
n_cli <- if (length(n_arg) > 0) as.integer(sub("^--n=", "", n_arg)) else NA_integer_
b_arg <- grep("^--B=", args, value = TRUE)
b_cli <- if (length(b_arg) > 0) as.integer(sub("^--B=", "", b_arg)) else NA_integer_

# Output directories
script_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error = function(e) normalizePath("paper")
)
data_dir <- file.path(script_dir, "data")
pkg_root <- dirname(script_dir)
fig_dir  <- file.path(pkg_root, "inst", "simulations", "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("=============================================================\n")
cat("Sensitivity Sweep: C1 Violation Severity vs C1-C2-C3 MLE\n")
cat("  5-component Weibull series system\n")
cat("  Tests whether total-hazard preservation extends from\n")
cat("  exponential (exact) to Weibull (first-order) under C1\n")
cat("=============================================================\n")
cat(sprintf("Quick mode:  %s\n", quick_mode))
cat(sprintf("Data dir:    %s\n", data_dir))
cat(sprintf("Figure dir:  %s\n", fig_dir))
cat("=============================================================\n\n")

# -----------------------------------------------------------------------------
# Weibull C1 violation data generator (inline; not in mdrelax package)
# -----------------------------------------------------------------------------

#' Generate Weibull series-system data under Bernoulli C1 violation.
#'
#' Mirrors rexp_series_md_relaxed_c1 from the package but uses Weibull
#' component lifetimes. The candidate set is generated component-by-component
#' with p_k(k) = 1 - alpha_C1 (failed component) and p_j(k) = base_p
#' (non-failed components). Empty candidate sets are returned as-is; the
#' caller is expected to filter them out for the standard MLE.
#'
#' @param n Sample size
#' @param shapes Weibull shapes (length m)
#' @param scales Weibull scales (length m)
#' @param alpha_C1 Probability that the failed component is NOT in the
#'        candidate set (in [0, 0.5])
#' @param base_p Off-diagonal inclusion probability (constant under C2)
#' @param tau Right-censoring time
#' @return list(t, delta, C, k)
rwei_series_md_relaxed_c1_inline <- function(n, shapes, scales, alpha_C1,
                                              base_p, tau) {
    m <- length(shapes)

    Tm <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
        Tm[, j] <- rweibull(n, shape = shapes[j], scale = scales[j])
    }

    sys_time <- apply(Tm, 1, min)
    k        <- apply(Tm, 1, which.min)

    tau_v    <- rep(tau, length.out = n)
    delta    <- as.numeric(sys_time <= tau_v)
    obs_time <- pmin(sys_time, tau_v)

    # Construct P with p_k(k) = 1 - alpha_C1, p_j(k) = base_p for j != k
    diag_p <- 1 - alpha_C1
    P <- matrix(base_p, nrow = m, ncol = m)
    diag(P) <- diag_p

    C <- matrix(FALSE, nrow = n, ncol = m)
    for (i in seq_len(n)) {
        if (delta[i] == 1) {
            ki <- k[i]
            for (j in seq_len(m)) {
                C[i, j] <- runif(1) <= P[j, ki]
            }
        }
    }

    list(t = obs_time, delta = delta, C = C, k = k)
}

# =============================================================================
# Weibull C1 Sensitivity Sweep
# =============================================================================

run_sweep <- function(quick_mode) {
    cat("\n=== Weibull C1 Sensitivity Sweep (m=5) ===\n")

    m         <- 5
    shapes    <- c(2.0, 1.5, 1.2, 1.8, 1.0)
    scales    <- c(3.0, 4.0, 5.0, 3.5, 4.5)
    tau       <- 5
    n         <- if (!is.na(n_cli)) n_cli else 500L
    B         <- if (!is.na(b_cli)) b_cli else (if (quick_mode) 10L else 200L)
    base_p    <- 0.5
    # alpha_C1 in [0, 0.5]: probability the failed component is excluded
    # from the candidate set. Above 0.5 the candidate set becomes anti-
    # informative and the modeling regime is qualitatively different.
    alpha_grid <- if (quick_mode) c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
                  else            seq(0, 0.5, by = 0.05)

    theta_true <- as.vector(rbind(shapes, scales))
    theta0     <- theta_true * 0.9

    t_ref      <- 2.0
    comp_names <- c(paste0(rep(c("k", "lambda"), m),
                           rep(seq_len(m), each = 2)),
                    "sys_hazard")
    true_haz   <- hazard_wei_series(t_ref, shapes, scales)

    cat(sprintf("  True params: shapes = (%s), scales = (%s)\n",
                paste(shapes, collapse = ", "), paste(scales, collapse = ", ")))
    cat(sprintf("  System hazard at t=%.1f: %.4f\n", t_ref, true_haz))
    cat(sprintf("  n = %d, B = %d, tau = %.1f, base_p = %.1f\n",
                n, B, tau, base_p))
    cat(sprintf("  alpha_C1 grid: %s\n\n",
                paste(sprintf("%.2f", alpha_grid), collapse = ", ")))

    all_results <- vector("list", length(alpha_grid))

    for (ai in seq_along(alpha_grid)) {
        alpha_C1 <- alpha_grid[ai]
        cat(sprintf("  alpha_C1 = %.2f: ", alpha_C1))

        theta_hat_mat <- matrix(NA_real_, nrow = B, ncol = 2 * m)
        se_mat        <- matrix(NA_real_, nrow = B, ncol = 2 * m)
        fim_mat_list  <- vector("list", B)
        converged     <- logical(B)
        n_retained_v  <- numeric(B)

        for (b in seq_len(B)) {
            if (b %% 10 == 0) cat(".")

            sim <- rwei_series_md_relaxed_c1_inline(
                n = n, shapes = shapes, scales = scales,
                alpha_C1 = alpha_C1, base_p = base_p, tau = tau
            )

            # Retention: drop observations with empty candidate sets among
            # failures (the standard C1-C2-C3 MLE is undefined on those).
            has_set <- rowSums(sim$C) >= 1
            keep    <- (sim$delta == 0) | has_set
            sim_kept <- list(
                t     = sim$t[keep],
                C     = sim$C[keep, , drop = FALSE],
                delta = sim$delta[keep]
            )
            n_retained_v[b] <- length(sim_kept$t)

            fit <- tryCatch(
                mle_wei_series(sim_kept$t, sim_kept$C, sim_kept$delta,
                               theta0 = theta0),
                error = function(e) NULL
            )

            if (!is.null(fit) && fit$converged) {
                converged[b]       <- TRUE
                theta_hat_mat[b, ] <- fit$theta
                se_mat[b, ]        <- fit$se
                fim_mat_list[[b]]  <- fit$fim
            }
        }
        cat(sprintf(" %d/%d converged  (mean retained = %.0f / %d)\n",
                    sum(converged), B, mean(n_retained_v), n))

        ok     <- which(converged)
        n_conv <- length(ok)

        rows <- vector("list", length(comp_names))
        for (ci in seq_along(comp_names)) {
            if (ci <= 2 * m) {
                ests <- theta_hat_mat[ok, ci]
                ses  <- se_mat[ok, ci]
                tv   <- theta_true[ci]
            } else {
                ests <- numeric(n_conv)
                ses  <- numeric(n_conv)
                for (r in seq_len(n_conv)) {
                    idx <- ok[r]
                    th  <- theta_hat_mat[idx, ]
                    k_hat <- th[seq(1, 2 * m, by = 2)]
                    l_hat <- th[seq(2, 2 * m, by = 2)]
                    ests[r] <- hazard_wei_series(t_ref, k_hat, l_hat)

                    h_fn <- function(theta) {
                        ks <- theta[seq(1, 2 * m, by = 2)]
                        ls <- theta[seq(2, 2 * m, by = 2)]
                        hazard_wei_series(t_ref, ks, ls)
                    }
                    grad_h <- tryCatch(numDeriv::grad(h_fn, th),
                                       error = function(e) rep(NA_real_, 2 * m))
                    cov_mat <- tryCatch(solve(fim_mat_list[[idx]]),
                                        error = function(e) matrix(NA_real_, 2 * m, 2 * m))
                    ses[r] <- tryCatch(
                        sqrt(as.numeric(t(grad_h) %*% cov_mat %*% grad_h)),
                        error = function(e) NA_real_
                    )
                }
                tv <- true_haz
            }

            bias     <- mean(ests) - tv
            rmse     <- sqrt(mean((ests - tv)^2))
            med_bias <- median(ests) - tv
            mad_est  <- mad(ests)

            coverage <- mean(ests - 1.96 * ses <= tv &
                             tv <= ests + 1.96 * ses, na.rm = TRUE)

            rows[[ci]] <- data.frame(
                alpha       = alpha_C1,
                component   = comp_names[ci],
                bias        = bias,
                rmse        = rmse,
                coverage    = coverage,
                median_bias = med_bias,
                mad         = mad_est,
                mean_est    = mean(ests),
                true_value  = tv,
                n_converged = n_conv,
                n_reps      = B,
                stringsAsFactors = FALSE
            )
        }
        all_results[[ai]] <- do.call(rbind, rows)
    }

    do.call(rbind, all_results)
}

# =============================================================================
# Figure Generation (reuses palette and plot style from C3 sweep)
# =============================================================================

get_palette <- function(n) {
    cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
              "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB")
    cols[seq_len(min(n, length(cols)))]
}

plot_metric <- function(df, metric, ylabel, fig_path, nominal_line = NULL) {
    sub <- df[!is.na(df[[metric]]), ]
    comps <- unique(sub$component)
    cols  <- get_palette(length(comps))

    pdf(fig_path, width = 8, height = 6)
    par(mar = c(4.5, 4.5, 2, 1), cex.lab = 1.2, cex.axis = 1)

    yvals <- sub[[metric]]
    ylim  <- range(yvals, na.rm = TRUE)
    if (!is.null(nominal_line)) ylim <- range(c(ylim, nominal_line))
    pad <- diff(ylim) * 0.05
    ylim <- ylim + c(-pad, pad)

    plot(NULL,
         xlim = range(sub$alpha),
         ylim = ylim,
         xlab = expression(paste("C1 Violation Severity (", alpha[C1], ")")),
         ylab = ylabel,
         main = "Weibull Series (m = 5)")

    for (ci in seq_along(comps)) {
        d <- sub[sub$component == comps[ci], ]
        lty <- if (grepl("^k", comps[ci])) 1 else if (grepl("^lambda", comps[ci])) 2 else 3
        pch <- if (grepl("sys", comps[ci])) 17 else 16
        lines(d$alpha, d[[metric]], col = cols[ci], lwd = 2, type = "o",
              pch = pch, cex = 0.8, lty = lty)
    }

    if (!is.null(nominal_line)) {
        abline(h = nominal_line, lty = 2, col = "gray40", lwd = 1.5)
    }

    legend("topleft", legend = comps, col = cols,
           lwd = 2, pch = ifelse(grepl("sys", comps), 17, 16),
           lty = ifelse(grepl("^k", comps), 1, ifelse(grepl("^lambda", comps), 2, 3)),
           cex = 0.7, bty = "n", ncol = 2)

    dev.off()
    cat(sprintf("  Saved: %s\n", fig_path))
}

# =============================================================================
# Slope test battery (first-order preservation diagnostic)
# =============================================================================

slope_test_battery <- function(df, label) {
    sys <- df[df$component == "sys_hazard", ]

    cat(sprintf("\n=== %s: Total Hazard Preservation Check ===\n", label))
    cat(sprintf("  True system hazard: %.4f\n", sys$true_value[1]))
    cat("  System hazard bias as function of alpha:\n")
    for (i in seq_len(nrow(sys))) {
        rel_bias <- sys$bias[i] / sys$true_value[i]
        cat(sprintf("    alpha = %.2f: bias = %+.5f  (%+.2f%% relative)\n",
                    sys$alpha[i], sys$bias[i], 100 * rel_bias))
    }

    # (1) Full-range linear
    if (nrow(sys) >= 3) {
        fit1 <- lm(bias ~ alpha, data = sys)
        cat("\n  (1) Full-range linear:\n")
        cat(sprintf("      intercept = %+.5f (p = %.3f)\n",
                    coef(fit1)[1], summary(fit1)$coefficients[1, 4]))
        cat(sprintf("      slope     = %+.5f (p = %.3f)\n",
                    coef(fit1)[2], summary(fit1)$coefficients[2, 4]))
    }

    # (2) Small-alpha linear
    cutoff <- max(sys$alpha) / 4
    small  <- sys[sys$alpha <= cutoff, ]
    if (nrow(small) >= 3) {
        fit2 <- lm(bias ~ alpha, data = small)
        cat(sprintf("\n  (2) Small-alpha linear (alpha <= %.2f, n = %d):\n",
                    cutoff, nrow(small)))
        cat(sprintf("      slope = %+.5f (p = %.3f)\n",
                    coef(fit2)[2], summary(fit2)$coefficients[2, 4]))
        cat("      Tests local first-order behavior at alpha = 0.\n")
    }

    # (3) Quadratic fit
    if (nrow(sys) >= 4) {
        fit3 <- lm(bias ~ alpha + I(alpha^2), data = sys)
        cs <- summary(fit3)$coefficients
        cat("\n  (3) Quadratic fit:\n")
        cat(sprintf("      linear term  = %+.5f (p = %.3f)\n", cs[2,1], cs[2,4]))
        cat(sprintf("      quadratic    = %+.5f (p = %.3f)\n", cs[3,1], cs[3,4]))
        cat("      First-order preservation: linear term not significantly nonzero.\n")
    }
}

# =============================================================================
# Main
# =============================================================================

results <- run_sweep(quick_mode)

n_used <- if (!is.na(n_cli)) n_cli else 500L
suffix <- if (n_used != 500L) sprintf("_n%d", n_used) else ""
rds_path <- file.path(data_dir, sprintf("sensitivity_sweep_c1%s.rds", suffix))
saveRDS(results, rds_path)
cat(sprintf("\nSaved: %s\n", rds_path))

slope_test_battery(results, "C1 Violation")

cat("\n=== Generating Figures ===\n")
plot_metric(results, "bias", "Bias",
            file.path(fig_dir, sprintf("fig_c1_bias_vs_alpha%s.pdf", suffix)))
plot_metric(results, "rmse", "RMSE",
            file.path(fig_dir, sprintf("fig_c1_rmse_vs_alpha%s.pdf", suffix)))
plot_metric(results, "coverage", "Coverage",
            file.path(fig_dir, sprintf("fig_c1_coverage_vs_alpha%s.pdf", suffix)),
            nominal_line = 0.95)
plot_metric(results, "median_bias", "Median Bias",
            file.path(fig_dir, sprintf("fig_c1_median_bias_vs_alpha%s.pdf", suffix)))
plot_metric(results, "mad", "MAD",
            file.path(fig_dir, sprintf("fig_c1_mad_vs_alpha%s.pdf", suffix)))

cat("\n=== C1 sensitivity sweep complete ===\n")
