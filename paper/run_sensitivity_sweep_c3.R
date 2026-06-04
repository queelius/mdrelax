#!/usr/bin/env Rscript
# =============================================================================
# Sensitivity Sweep: How does C3 violation severity affect C1-C2-C3 MLE?
# =============================================================================
#
# Generates data under parameter-dependent masking (relaxed C3, power-weight
# model) with varying severity alpha, fits the misspecified C1-C2-C3 model,
# and measures bias, RMSE, and coverage as a function of alpha.
#
# Theoretical question tested: For Weibull components, is the total system
# hazard preserved at first order under C3 violation?
#
# Background: Theorem (Exact Total-Hazard Preservation, Exponential Components)
# in sensitivity_framework.tex shows that for exponential components, the
# C1-C2-C3 MLE preserves the total hazard exactly under any masking violation,
# including C3. The argument rests on the additive factorization of the
# exponential log-likelihood as ell(S, phi) = ell_S(S) + ell_phi(phi), which
# does not hold for Weibull components. This sweep empirically tests whether
# Weibull total-hazard preservation holds at least at first order under C3,
# or whether the parameter-dependent masking biases the time data enough to
# break first-order robustness.
#
# System: 5-component Weibull series (matches the C2 sweep for comparability)
#   shapes = (2.0, 1.5, 1.2, 1.8, 1.0)
#   scales = (3.0, 4.0, 5.0, 3.5, 4.5)
#   tau = 5, n = 500, B = 200
#
# Usage:
#   Rscript paper/run_sensitivity_sweep_c3.R [--quick]
#
# Outputs:
#   paper/data/sensitivity_sweep_c3.rds
#   inst/simulations/figures/fig_c3_bias_vs_alpha.pdf
#   inst/simulations/figures/fig_c3_rmse_vs_alpha.pdf
#   inst/simulations/figures/fig_c3_coverage_vs_alpha.pdf
#   inst/simulations/figures/fig_c3_median_bias_vs_alpha.pdf
#   inst/simulations/figures/fig_c3_mad_vs_alpha.pdf
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
cat("Sensitivity Sweep: C3 Violation Severity vs C1-C2-C3 MLE\n")
cat("  5-component Weibull series system\n")
cat("  Tests whether total-hazard preservation holds for Weibull\n")
cat("  under parameter-dependent masking\n")
cat("=============================================================\n")
cat(sprintf("Quick mode:  %s\n", quick_mode))
cat(sprintf("Data dir:    %s\n", data_dir))
cat(sprintf("Figure dir:  %s\n", fig_dir))
cat("=============================================================\n\n")

# =============================================================================
# Weibull C3 Sensitivity Sweep
# =============================================================================

run_sweep <- function(quick_mode) {
    cat("\n=== Weibull C3 Sensitivity Sweep (m=5) ===\n")

    m         <- 5
    shapes    <- c(2.0, 1.5, 1.2, 1.8, 1.0)
    scales    <- c(3.0, 4.0, 5.0, 3.5, 4.5)
    tau       <- 5
    n         <- if (!is.na(n_cli)) n_cli else 500L
    B         <- if (!is.na(b_cli)) b_cli else (if (quick_mode) 10L else 200L)
    base_p    <- 0.5
    # alpha = 0 is C3-respecting (uniform masking); higher values weight
    # high-rate components more heavily. Practical sweep range: [0, 2].
    alpha_grid <- if (quick_mode) c(0, 0.5, 1.0, 1.5, 2.0)
                  else            seq(0, 2.0, by = 0.25)

    theta_true <- as.vector(rbind(shapes, scales))  # k1,l1,k2,l2,...,k5,l5
    theta0     <- theta_true * 0.9

    # Reference time for system hazard (approx median; matches C2 sweep)
    t_ref <- 2.0

    comp_names <- c(paste0(rep(c("k", "lambda"), m),
                           rep(seq_len(m), each = 2)),
                    "sys_hazard")
    true_haz   <- hazard_wei_series(t_ref, shapes, scales)

    cat(sprintf("  True params: shapes = (%s), scales = (%s)\n",
                paste(shapes, collapse = ", "), paste(scales, collapse = ", ")))
    cat(sprintf("  System hazard at t=%.1f: %.4f\n", t_ref, true_haz))
    cat(sprintf("  n = %d, B = %d, tau = %.1f, base_p = %.1f\n",
                n, B, tau, base_p))
    cat(sprintf("  alpha grid: %s\n\n",
                paste(sprintf("%.2f", alpha_grid), collapse = ", ")))

    all_results <- vector("list", length(alpha_grid))

    for (ai in seq_along(alpha_grid)) {
        alpha <- alpha_grid[ai]
        cat(sprintf("  alpha = %.2f: ", alpha))

        theta_hat_mat <- matrix(NA_real_, nrow = B, ncol = 2 * m)
        se_mat        <- matrix(NA_real_, nrow = B, ncol = 2 * m)
        fim_mat_list  <- vector("list", B)
        converged     <- logical(B)

        for (b in seq_len(B)) {
            if (b %% 10 == 0) cat(".")

            sim <- rwei_series_md_c1_c2(
                n = n, shapes = shapes, scales = scales,
                alpha = alpha, base_p = base_p, tau = tau
            )

            fit <- tryCatch(
                mle_wei_series(sim$t, sim$C, sim$delta, theta0 = theta0),
                error = function(e) NULL
            )

            if (!is.null(fit) && fit$converged) {
                converged[b]       <- TRUE
                theta_hat_mat[b, ] <- fit$theta
                se_mat[b, ]        <- fit$se
                fim_mat_list[[b]]  <- fit$fim
            }
        }
        cat(sprintf(" %d/%d converged\n", sum(converged), B))

        ok     <- which(converged)
        n_conv <- length(ok)

        rows <- vector("list", length(comp_names))
        for (ci in seq_along(comp_names)) {
            if (ci <= 2 * m) {
                ests <- theta_hat_mat[ok, ci]
                ses  <- se_mat[ok, ci]
                tv   <- theta_true[ci]
            } else {
                # System hazard at t_ref via delta method
                ests <- numeric(n_conv)
                ses  <- numeric(n_conv)
                for (r in seq_len(n_conv)) {
                    idx <- ok[r]
                    th  <- theta_hat_mat[idx, ]
                    k_hat <- th[seq(1, 2 * m, by = 2)]
                    l_hat <- th[seq(2, 2 * m, by = 2)]
                    ests[r] <- hazard_wei_series(t_ref, k_hat, l_hat)

                    # Delta method SE for system hazard
                    h_fn <- function(theta) {
                        ks <- theta[seq(1, 2 * m, by = 2)]
                        ls <- theta[seq(2, 2 * m, by = 2)]
                        hazard_wei_series(t_ref, ks, ls)
                    }
                    grad_h <- tryCatch(
                        numDeriv::grad(h_fn, th),
                        error = function(e) rep(NA_real_, 2 * m)
                    )
                    cov_mat <- tryCatch(
                        solve(fim_mat_list[[idx]]),
                        error = function(e) matrix(NA_real_, 2 * m, 2 * m)
                    )
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
                alpha       = alpha,
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
# Figure Generation
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
         xlab = expression(paste("C3 Violation Severity (", alpha, ")")),
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
# Main
# =============================================================================

results <- run_sweep(quick_mode)

# Save RDS
n_used <- if (!is.na(n_cli)) n_cli else 500L
suffix <- if (n_used != 500L) sprintf("_n%d", n_used) else ""
rds_path <- file.path(data_dir, sprintf("sensitivity_sweep_c3%s.rds", suffix))
saveRDS(results, rds_path)
cat(sprintf("\nSaved: %s\n", rds_path))

# Summary: total hazard preservation check
cat("\n=== Total Hazard Preservation Check ===\n")
sys_results <- results[results$component == "sys_hazard", ]
cat(sprintf("  Reference time: t = 2.0, True hazard: %.4f\n",
            sys_results$true_value[1]))
cat("  System hazard bias as function of alpha:\n")
for (i in seq_len(nrow(sys_results))) {
    rel_bias <- sys_results$bias[i] / sys_results$true_value[i]
    cat(sprintf("    alpha = %.2f: bias = %+.5f  (%+.2f%% relative)\n",
                sys_results$alpha[i], sys_results$bias[i], 100 * rel_bias))
}

# Slope test for first-order preservation
if (nrow(sys_results) >= 3) {
    fit_lm <- lm(bias ~ alpha, data = sys_results)
    cat(sprintf("\n  Linear fit of system hazard bias vs alpha:\n"))
    cat(sprintf("    intercept = %+.5f (p = %.3f)\n",
                coef(fit_lm)[1], summary(fit_lm)$coefficients[1, 4]))
    cat(sprintf("    slope     = %+.5f (p = %.3f)\n",
                coef(fit_lm)[2], summary(fit_lm)$coefficients[2, 4]))
    cat("  If slope is significantly nonzero, first-order preservation fails.\n")
    cat("  If slope is indistinguishable from zero, first-order preservation\n")
    cat("  holds even for Weibull under C3 (matching the exponential case).\n")
}

# Generate figures
cat("\n=== Generating Figures ===\n")

plot_metric(results, "bias", "Bias",
            file.path(fig_dir, sprintf("fig_c3_bias_vs_alpha%s.pdf", suffix)))

plot_metric(results, "rmse", "RMSE",
            file.path(fig_dir, sprintf("fig_c3_rmse_vs_alpha%s.pdf", suffix)))

plot_metric(results, "coverage", "Coverage",
            file.path(fig_dir, sprintf("fig_c3_coverage_vs_alpha%s.pdf", suffix)),
            nominal_line = 0.95)

plot_metric(results, "median_bias", "Median Bias",
            file.path(fig_dir, sprintf("fig_c3_median_bias_vs_alpha%s.pdf", suffix)))

plot_metric(results, "mad", "MAD",
            file.path(fig_dir, sprintf("fig_c3_mad_vs_alpha%s.pdf", suffix)))

cat("\n=== C3 sensitivity sweep complete ===\n")
