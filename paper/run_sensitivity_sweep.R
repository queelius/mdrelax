#!/usr/bin/env Rscript
# =============================================================================
# Sensitivity Sweep: How does C2 violation severity affect C1-C2-C3 MLE?
# =============================================================================
#
# Generates data under informative masking (relaxed C2) with varying severity,
# fits the misspecified C1-C2-C3 model, and measures bias, RMSE, and coverage
# as a function of severity.
#
# System: 5-component Weibull series
#   shapes = (2.0, 1.5, 1.2, 1.8, 1.0)
#   scales = (3.0, 4.0, 5.0, 3.5, 4.5)
#   tau = 5, n = 500, B = 200
#
# Usage:
#   Rscript paper/run_sensitivity_sweep.R [--quick]
#
# Outputs:
#   paper/data/sensitivity_sweep.rds
#   inst/simulations/figures/fig_bias_vs_severity.pdf
#   inst/simulations/figures/fig_rmse_vs_severity.pdf
#   inst/simulations/figures/fig_coverage_vs_severity.pdf
#
# =============================================================================

devtools::load_all(quiet = TRUE)

set.seed(42)

# -----------------------------------------------------------------------------
# Parse Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args

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
cat("Sensitivity Sweep: C2 Violation Severity vs C1-C2-C3 MLE\n")
cat("  5-component Weibull series system\n")
cat("=============================================================\n")
cat(sprintf("Quick mode:  %s\n", quick_mode))
cat(sprintf("Data dir:    %s\n", data_dir))
cat(sprintf("Figure dir:  %s\n", fig_dir))
cat("=============================================================\n\n")

# -----------------------------------------------------------------------------
# Direction matrix for P(s) = P_uniform + s * D
# -----------------------------------------------------------------------------

#' Create a 5x5 P matrix parameterized by severity s in [0, 1]
#'
#' P(s) = P_uniform(base_p) + s * D, with off-diag clamped to [0.05, 0.95]
#' and diagonal fixed at 1.
#'
#' The direction matrix D has a cyclic structure with rows summing to zero,
#' meaning the perturbation redistributes masking probability across
#' non-failed components without changing the total.
#'
#' @param s Severity parameter in [0, 1]
#' @param base_p Base off-diagonal probability for the uniform P matrix
#' @return 5 x 5 P matrix
make_sweep_P <- function(s, base_p) {
    D <- matrix(c(
         0,    0.3, -0.2,  0.1, -0.2,
        -0.2,  0,    0.3, -0.2,  0.1,
         0.1, -0.2,  0,    0.3, -0.2,
        -0.2,  0.1, -0.2,  0,    0.3,
         0.3, -0.2,  0.1, -0.2,  0
    ), nrow = 5, byrow = TRUE)

    P <- make_P_matrix(5, type = "uniform", p = base_p)
    P <- P + s * D

    # Clamp off-diagonal to [0.05, 0.95], keep diagonal = 1
    offdiag <- row(P) != col(P)
    P[offdiag] <- pmax(0.05, pmin(0.95, P[offdiag]))
    diag(P) <- 1

    P
}

# =============================================================================
# Weibull Sensitivity Sweep
# =============================================================================

run_sweep <- function(quick_mode) {
    cat("\n=== Weibull Sensitivity Sweep (m=5) ===\n")

    m       <- 5
    shapes  <- c(2.0, 1.5, 1.2, 1.8, 1.0)
    scales  <- c(3.0, 4.0, 5.0, 3.5, 4.5)
    tau     <- 5
    n       <- 500
    B       <- if (quick_mode) 10 else 200
    base_p  <- 0.5
    s_grid  <- seq(0, 1, by = 0.1)

    theta_true <- as.vector(rbind(shapes, scales))  # k1,l1,k2,l2,...,k5,l5
    theta0     <- theta_true * 0.9

    # Reference time for system hazard (approx median)
    t_ref <- 2.0

    comp_names <- c(paste0(rep(c("k", "lambda"), m),
                           rep(seq_len(m), each = 2)),
                    "sys_hazard")
    true_haz   <- hazard_wei_series(t_ref, shapes, scales)

    cat(sprintf("  True params: shapes = (%s), scales = (%s)\n",
                paste(shapes, collapse = ", "), paste(scales, collapse = ", ")))
    cat(sprintf("  System hazard at t=%.1f: %.4f\n", t_ref, true_haz))
    cat(sprintf("  n = %d, B = %d, tau = %.1f, base_p = %.1f\n\n", n, B, tau, base_p))

    all_results <- vector("list", length(s_grid))

    for (si in seq_along(s_grid)) {
        s <- s_grid[si]
        P <- make_sweep_P(s, base_p)
        cat(sprintf("  s = %.1f: ", s))

        theta_hat_mat <- matrix(NA_real_, nrow = B, ncol = 2 * m)
        se_mat        <- matrix(NA_real_, nrow = B, ncol = 2 * m)
        converged     <- logical(B)

        for (b in seq_len(B)) {
            if (b %% 10 == 0) cat(".")

            sim <- rwei_series_md_c1_c3(
                n = n, shapes = shapes, scales = scales, P = P, tau = tau
            )

            fit <- tryCatch(
                mle_wei_series(sim$t, sim$C, sim$delta, theta0 = theta0),
                error = function(e) NULL
            )

            if (!is.null(fit) && fit$converged) {
                converged[b]       <- TRUE
                theta_hat_mat[b, ] <- fit$theta
                se_mat[b, ]        <- fit$se
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
                # System hazard at t_ref
                ests <- numeric(n_conv)
                ses  <- rep(NA_real_, n_conv)
                for (r in seq_len(n_conv)) {
                    idx <- ok[r]
                    k_hat <- theta_hat_mat[idx, seq(1, 2 * m, by = 2)]
                    l_hat <- theta_hat_mat[idx, seq(2, 2 * m, by = 2)]
                    ests[r] <- hazard_wei_series(t_ref, k_hat, l_hat)
                }
                tv <- true_haz
            }

            bias <- mean(ests) - tv
            rmse <- sqrt(mean((ests - tv)^2))

            if (ci <= 2 * m) {
                coverage <- mean(ests - 1.96 * ses <= tv &
                                 tv <= ests + 1.96 * ses)
            } else {
                coverage <- NA_real_
            }

            rows[[ci]] <- data.frame(
                severity    = s,
                component   = comp_names[ci],
                bias        = bias,
                rmse        = rmse,
                coverage    = coverage,
                mean_est    = mean(ests),
                true_value  = tv,
                n_converged = n_conv,
                n_reps      = B,
                stringsAsFactors = FALSE
            )
        }
        all_results[[si]] <- do.call(rbind, rows)
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

    plot(NULL, xlim = c(0, 1), ylim = ylim,
         xlab = "C2 Violation Severity (s)", ylab = ylabel,
         main = "Weibull Series (m = 5)")

    for (ci in seq_along(comps)) {
        d <- sub[sub$component == comps[ci], ]
        lty <- if (grepl("^k", comps[ci])) 1 else if (grepl("^lambda", comps[ci])) 2 else 3
        pch <- if (grepl("sys", comps[ci])) 17 else 16
        lines(d$severity, d[[metric]], col = cols[ci], lwd = 2, type = "o",
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
rds_path <- file.path(data_dir, "sensitivity_sweep.rds")
saveRDS(results, rds_path)
cat(sprintf("\nSaved: %s\n", rds_path))

# Generate figures
cat("\n=== Generating Figures ===\n")

plot_metric(results, "bias", "Bias",
            file.path(fig_dir, "fig_bias_vs_severity.pdf"))

plot_metric(results, "rmse", "RMSE",
            file.path(fig_dir, "fig_rmse_vs_severity.pdf"))

plot_metric(results, "coverage", "Coverage",
            file.path(fig_dir, "fig_coverage_vs_severity.pdf"),
            nominal_line = 0.95)

cat("\n=== Sensitivity sweep complete ===\n")
