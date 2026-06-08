#!/usr/bin/env Rscript
# =============================================================================
# Regenerate the Application section numbers (Guo et al. turbine-engine data)
# from the package's canonical objects, so the manuscript stays consistent
# with the shipped package.
#
# Produces:
#   1. The s = 0 baseline row of the Guo sensitivity table, taken directly
#      from guo_weibull_series_mle (the package's committed canonical fit),
#      so the printed baseline scales match the package and not a stale rerun.
#   2. A worked first-order robustness interval for the system hazard under a
#      C2 violation, via ri_first_order(). This exercises contribution 3 on
#      data. The seed is pinned so the printed RI is exactly reproducible.
#
# n = 30 is small, so the robustness interval is an illustration of the
# workflow, not a precise estimate; the qualitative verdict (the system
# hazard tolerates substantial C2 severity) is what the table and the
# simulation study independently corroborate.
#
# Usage:
#   Rscript paper/regenerate_application.R
# =============================================================================

suppressMessages({
    if (requireNamespace("devtools", quietly = TRUE)) {
        devtools::load_all(file.path(dirname(
            tryCatch(normalizePath(sys.frame(1)$ofile),
                     error = function(e) "paper/x")), ".."),
            quiet = TRUE)
    } else {
        library(mdrelax)
    }
})

# -----------------------------------------------------------------------------
# Canonical baseline (s = 0 row)
# -----------------------------------------------------------------------------

mle    <- guo_weibull_series_mle$mle
shapes <- mle[c(1, 3, 5)]
scales <- mle[c(2, 4, 6)]

t0 <- 249.5                      # median failure time used in the manuscript
h0 <- hazard_wei_series(t0, shapes, scales)

cat("=== Guo et al. canonical baseline (s = 0) ===\n")
cat(sprintf("  shapes  k_hat      = (%.2f, %.2f, %.2f)\n",
            shapes[1], shapes[2], shapes[3]))
cat(sprintf("  scales  lambda_hat = (%.0f, %.0f, %.0f)\n",
            scales[1], scales[2], scales[3]))
cat(sprintf("  system hazard h_T(%.1f) = %.5f\n\n", t0, h0))

# -----------------------------------------------------------------------------
# First-order robustness interval for the system hazard under C2
# -----------------------------------------------------------------------------

tol_frac <- 0.05                 # tolerance: 5% of the baseline system hazard
tolerance <- tol_frac * h0

set.seed(20260608)
ri <- suppressWarnings(ri_first_order(
    shapes_hat  = shapes,
    scales_hat  = scales,
    estimand_fn = function(s, l) hazard_wei_series(t0, s, l),
    violation   = "C2",
    tolerance   = tolerance,
    n           = 30,
    B           = 300,
    tau         = 2000,
    probe_alpha = c(0.05, 0.10)
))

cat("=== First-order robustness interval (system hazard, C2) ===\n")
cat(sprintf("  tolerance         = %.0f%% of baseline (= %.6f)\n",
            100 * tol_frac, tolerance))
cat(sprintf("  ISNI (slope at 0) = %+.5f\n", ri$ISNI))
cat(sprintf("  robustness interval RI = %.2f\n", ri$RI))
cat(sprintf(
    "  -> the system-hazard estimate stays within %.0f%% of baseline for C2\n",
    100 * tol_frac))
cat(sprintf(
    "     severities up to s = %.2f (well beyond the sweep range s in [0,1]).\n",
    ri$RI))
