#!/usr/bin/env Rscript
# =============================================================================
# Specification Test Validation: Joint Null {C1, C2, C3}
# =============================================================================
#
# Implements and validates the singleton specification test sketched in
# paper/sections/specification_test_sketch.tex. The test statistic
#
#   T_n = G^T V^- G   with   G_j = sum_i [ 1{C_i = {j}} - pi_j R_j(T_i) ]
#
# should be asymptotically chi-squared on (m-1) degrees of freedom under
# H_0: {C1, C2, C3}. This script verifies the null distribution and
# computes power against parametric C1, C2, C3 alternatives.
#
# System: 5-component Weibull series (matches the C1, C2, C3 sweeps)
#   shapes = (2.0, 1.5, 1.2, 1.8, 1.0)
#   scales = (3.0, 4.0, 5.0, 3.5, 4.5)
#   tau = 5, n = 500, B = 200
#
# Usage:
#   Rscript paper/specification_test_validation.R [--quick]
#
# Outputs:
#   paper/data/spec_test_null.rds                (null-distribution sample)
#   paper/data/spec_test_power_c1.rds            (power vs alpha_C1)
#   paper/data/spec_test_power_c2.rds            (power vs s_C2)
#   paper/data/spec_test_power_c3.rds            (power vs alpha_C3)
#   inst/simulations/figures/fig_spec_test_null_qq.pdf
#   inst/simulations/figures/fig_spec_test_power.pdf
#
# =============================================================================

devtools::load_all(quiet = TRUE)
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args

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
cat("Specification Test Validation: H_0 = {C1, C2, C3}\n")
cat("  Singleton statistic, chi-squared on (m-1) df under null\n")
cat("  5-component Weibull series system\n")
cat("=============================================================\n\n")

# Shared system
m         <- 5
shapes    <- c(2.0, 1.5, 1.2, 1.8, 1.0)
scales    <- c(3.0, 4.0, 5.0, 3.5, 4.5)
tau       <- 5
n         <- 500
base_p    <- 0.5
B         <- if (quick_mode) 100 else 500

theta_true <- as.vector(rbind(shapes, scales))
theta0     <- theta_true * 0.9

# -----------------------------------------------------------------------------
# Specification test statistic (singleton version)
# -----------------------------------------------------------------------------

#' Compute the singleton specification test statistic via bin-based
#' Pearson chi-squared.
#'
#' Tests whether Pr(C = {j} | T = t) factors as pi_j * R_j(t) for some
#' constant pi_j (the testable content of {C1, C2, C3} per the
#' factorization lemma). The test bins observations by time and compares
#' observed singleton-{j} counts per bin against expected counts under
#' the factored form. Aggregate-level moment estimators absorb the global
#' rate but leave bin-level deviations that the test detects.
#'
#' @param t Vector of (possibly censored) observation times.
#' @param C Logical matrix (n x m) of candidate sets.
#' @param delta Censoring indicators (1 = failure, 0 = censored).
#' @param shapes_hat Weibull shape MLE estimates (length m).
#' @param scales_hat Weibull scale MLE estimates (length m).
#' @param L Number of time bins (default 5; quantile-based).
#' @return list(T_stat, df, p_value, pi_hat, n_failures, L)
spec_test_singleton <- function(t, C, delta, shapes_hat, scales_hat, L = 5) {
    m  <- length(shapes_hat)
    ok <- which(delta == 1)
    if (length(ok) < m * L + 5) {
        return(list(T_stat = NA_real_, df = NA_integer_, p_value = NA_real_,
                    pi_hat = rep(NA_real_, m), n_failures = length(ok),
                    L = L))
    }
    t_f <- t[ok]
    C_f <- C[ok, , drop = FALSE]
    n_f <- length(t_f)

    # Component hazards at failure times
    H <- matrix(NA_real_, nrow = n_f, ncol = m)
    for (j in seq_len(m)) {
        H[, j] <- (shapes_hat[j] / scales_hat[j]) *
                  (t_f / scales_hat[j])^(shapes_hat[j] - 1)
    }
    S_t      <- rowSums(H)
    R_single <- H / S_t

    # Identify singletons
    set_size  <- rowSums(C_f)
    is_single <- set_size == 1
    obs_j     <- rep(NA_integer_, n_f)
    for (i in which(is_single)) {
        w <- which(C_f[i, ])
        if (length(w) == 1) obs_j[i] <- w
    }

    # Bin failures by time quantiles
    bin_breaks <- unique(quantile(t_f, probs = seq(0, 1, length.out = L + 1)))
    if (length(bin_breaks) < L) {
        return(list(T_stat = NA_real_, df = NA_integer_, p_value = NA_real_,
                    pi_hat = rep(NA_real_, m), n_failures = n_f,
                    L = length(bin_breaks) - 1))
    }
    L_actual <- length(bin_breaks) - 1
    bin_idx  <- cut(t_f, breaks = bin_breaks, include.lowest = TRUE,
                    labels = FALSE)

    # Global moment estimator for pi_j (constant across bins under H_0)
    pi_hat <- numeric(m)
    for (j in seq_len(m)) {
        denom_global <- sum(R_single[, j])
        if (denom_global > 0) {
            pi_hat[j] <- sum(obs_j == j, na.rm = TRUE) / denom_global
        }
    }

    # Pearson chi-squared: bin-level deviations from the factored prediction
    T_stat <- 0
    valid_terms <- 0
    for (j in seq_len(m)) {
        if (pi_hat[j] == 0) next
        for (l in seq_len(L_actual)) {
            in_bin <- bin_idx == l
            n_obs  <- sum(obs_j[in_bin] == j, na.rm = TRUE)
            e_exp  <- pi_hat[j] * sum(R_single[in_bin, j])
            if (e_exp >= 0.5) {
                T_stat <- T_stat + (n_obs - e_exp)^2 / e_exp
                valid_terms <- valid_terms + 1
            }
        }
    }

    # Approximate df: m * (L - 1), reduced by valid_terms shortfall
    df_base <- m * (L_actual - 1)
    df      <- max(1, min(valid_terms - m, df_base))
    p_value <- pchisq(T_stat, df = df, lower.tail = FALSE)

    list(T_stat = T_stat, df = df, p_value = p_value,
         pi_hat = pi_hat, n_failures = n_f, L = L_actual)
}

# -----------------------------------------------------------------------------
# Data generators for null and three alternatives
# -----------------------------------------------------------------------------

# Null: data under C1+C2+C3 (uniform Bernoulli with p_k(k) = 1)
gen_null <- function() {
    P_uni <- make_P_matrix(m, "uniform", p = base_p)
    sim <- rwei_series_md_c1_c3(n = n, shapes = shapes, scales = scales,
                                  P = P_uni, tau = tau)
    sim
}

# C1 alternative: relaxed C1 with p_k(k) = 1 - alpha_C1; failures with
# empty candidate sets are dropped to match practitioner workflow.
gen_c1 <- function(alpha_C1) {
    Tm <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
        Tm[, j] <- rweibull(n, shape = shapes[j], scale = scales[j])
    }
    sys_time <- apply(Tm, 1, min)
    k_true   <- apply(Tm, 1, which.min)
    tau_v    <- rep(tau, length.out = n)
    delta    <- as.numeric(sys_time <= tau_v)
    obs_time <- pmin(sys_time, tau_v)

    diag_p <- 1 - alpha_C1
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
    list(t = obs_time[keep], delta = delta[keep],
         C = C[keep, , drop = FALSE], k = k_true[keep])
}

# C2 alternative: directional P perturbation (matches main paper's sweep)
gen_c2 <- function(s_C2) {
    D <- matrix(c(
         0,    0.3, -0.2,  0.1, -0.2,
        -0.2,  0,    0.3, -0.2,  0.1,
         0.1, -0.2,  0,    0.3, -0.2,
        -0.2,  0.1, -0.2,  0,    0.3,
         0.3, -0.2,  0.1, -0.2,  0
    ), nrow = 5, byrow = TRUE)
    P <- make_P_matrix(m, "uniform", p = base_p) + s_C2 * D
    offdiag <- row(P) != col(P)
    P[offdiag] <- pmax(0.05, pmin(0.95, P[offdiag]))
    diag(P) <- 1
    rwei_series_md_c1_c3(n = n, shapes = shapes, scales = scales,
                          P = P, tau = tau)
}

# C3 alternative: power-weight masking
gen_c3 <- function(alpha_C3) {
    rwei_series_md_c1_c2(n = n, shapes = shapes, scales = scales,
                          alpha = alpha_C3, base_p = base_p, tau = tau)
}

# -----------------------------------------------------------------------------
# Single-replicate runner: fit MLE under H_0, compute test statistic
# -----------------------------------------------------------------------------

run_one <- function(sim) {
    fit <- tryCatch(
        mle_wei_series(sim$t, sim$C, sim$delta, theta0 = theta0),
        error = function(e) NULL
    )
    if (is.null(fit) || !fit$converged) return(NULL)
    k_hat <- fit$theta[seq(1, 2 * m, by = 2)]
    l_hat <- fit$theta[seq(2, 2 * m, by = 2)]
    spec_test_singleton(sim$t, sim$C, sim$delta, k_hat, l_hat)
}

# -----------------------------------------------------------------------------
# 1. Null distribution
# -----------------------------------------------------------------------------

cat("=== Null distribution validation ===\n")
cat(sprintf("  Generating B = %d replicates under H_0 ...\n", B))

T_null <- numeric(B)
df_used <- NA_integer_
n_valid <- 0
for (b in seq_len(B)) {
    if (b %% 50 == 0) cat(sprintf("    rep %d / %d\n", b, B))
    res <- run_one(gen_null())
    if (!is.null(res) && !is.na(res$T_stat)) {
        n_valid <- n_valid + 1
        T_null[n_valid] <- res$T_stat
        df_used <- res$df
    }
}
T_null <- T_null[seq_len(n_valid)]
cat(sprintf("  Valid replicates: %d / %d\n", n_valid, B))

saveRDS(list(T_null = T_null, df = df_used, n_valid = n_valid),
        file.path(data_dir, "spec_test_null.rds"))

# Compare empirical quantiles to chi-squared
qs   <- c(0.50, 0.75, 0.90, 0.95, 0.99)
emp  <- quantile(T_null, probs = qs)
ref  <- qchisq(qs, df = df_used)
cat("\n  Quantile comparison (empirical vs chi-squared on", df_used, "df):\n")
for (i in seq_along(qs)) {
    cat(sprintf("    q = %.2f: empirical = %6.3f, chi-sq = %6.3f, ratio = %.3f\n",
                qs[i], emp[i], ref[i], emp[i] / ref[i]))
}

# Kolmogorov-Smirnov test of fit
ks_p <- suppressWarnings(ks.test(T_null, "pchisq", df_used))$p.value
cat(sprintf("\n  Kolmogorov-Smirnov vs chi-squared(%d): p = %.3f\n",
            df_used, ks_p))
cat("  (large p indicates the empirical distribution matches chi-squared)\n")

# Type I error at nominal 0.05
crit_005 <- qchisq(0.95, df = df_used)
type1    <- mean(T_null > crit_005)
cat(sprintf("\n  Empirical Type I error at nominal 0.05: %.3f\n", type1))

# Null QQ plot
pdf(file.path(fig_dir, "fig_spec_test_null_qq.pdf"), width = 6, height = 6)
par(mar = c(4.5, 4.5, 2, 1), cex.lab = 1.2, cex.axis = 1)
theo <- qchisq(ppoints(length(T_null)), df = df_used)
plot(theo, sort(T_null), pch = 16, cex = 0.5,
     xlab = bquote(paste("Theoretical ", chi[.(df_used)]^2, " quantile")),
     ylab = bquote(paste("Empirical ", T[n], " quantile")),
     main = "Null-distribution Q-Q plot")
abline(0, 1, col = "red", lwd = 2, lty = 2)
dev.off()
cat(sprintf("\n  Saved: %s\n",
            file.path(fig_dir, "fig_spec_test_null_qq.pdf")))

# -----------------------------------------------------------------------------
# 2. Power against C1, C2, C3 alternatives
# -----------------------------------------------------------------------------

cat("\n=== Power against parametric alternatives ===\n")

if (quick_mode) {
    alpha_C1_grid <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
    s_C2_grid     <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
    alpha_C3_grid <- c(0, 0.5, 1.0, 1.5, 2.0)
} else {
    alpha_C1_grid <- seq(0, 0.5, by = 0.05)
    s_C2_grid     <- seq(0, 1.0, by = 0.1)
    alpha_C3_grid <- seq(0, 2.0, by = 0.25)
}

power_sweep <- function(grid, gen_fn, label) {
    cat(sprintf("\n  -- %s --\n", label))
    out <- data.frame()
    for (alpha in grid) {
        cat(sprintf("    %s = %.2f: ", label, alpha))
        rej <- 0
        valid <- 0
        for (b in seq_len(B)) {
            if (b %% 50 == 0) cat(".")
            res <- run_one(gen_fn(alpha))
            if (!is.null(res) && !is.na(res$T_stat)) {
                valid <- valid + 1
                if (res$T_stat > crit_005) rej <- rej + 1
            }
        }
        power <- if (valid > 0) rej / valid else NA_real_
        cat(sprintf(" %d/%d valid, power = %.3f\n", valid, B, power))
        out <- rbind(out, data.frame(
            severity = alpha, label = label, n_valid = valid,
            n_reject = rej, power = power, stringsAsFactors = FALSE
        ))
    }
    out
}

pwr_c1 <- power_sweep(alpha_C1_grid, gen_c1, "alpha_C1")
pwr_c2 <- power_sweep(s_C2_grid,     gen_c2, "s_C2")
pwr_c3 <- power_sweep(alpha_C3_grid, gen_c3, "alpha_C3")

saveRDS(pwr_c1, file.path(data_dir, "spec_test_power_c1.rds"))
saveRDS(pwr_c2, file.path(data_dir, "spec_test_power_c2.rds"))
saveRDS(pwr_c3, file.path(data_dir, "spec_test_power_c3.rds"))

# Power figure
pdf(file.path(fig_dir, "fig_spec_test_power.pdf"), width = 8, height = 6)
par(mar = c(4.5, 4.5, 2, 1), cex.lab = 1.2, cex.axis = 1)
plot(NULL, xlim = c(0, max(c(alpha_C1_grid, s_C2_grid, alpha_C3_grid))),
     ylim = c(0, 1),
     xlab = "Violation severity",
     ylab = "Power (at 0.05 nominal level)",
     main = "Specification test power")
lines(pwr_c1$severity, pwr_c1$power, col = "#E41A1C", lwd = 2,
      type = "o", pch = 16)
lines(pwr_c2$severity, pwr_c2$power, col = "#377EB8", lwd = 2,
      type = "o", pch = 17)
lines(pwr_c3$severity, pwr_c3$power, col = "#4DAF4A", lwd = 2,
      type = "o", pch = 18)
abline(h = 0.05, lty = 2, col = "gray40")
legend("bottomright",
       legend = c("C1 violation", "C2 violation", "C3 violation",
                  "Nominal 0.05"),
       col = c("#E41A1C", "#377EB8", "#4DAF4A", "gray40"),
       lwd = c(2, 2, 2, 1), lty = c(1, 1, 1, 2),
       pch = c(16, 17, 18, NA), bty = "n")
dev.off()
cat(sprintf("\n  Saved: %s\n",
            file.path(fig_dir, "fig_spec_test_power.pdf")))

cat("\n=== Specification test validation complete ===\n")
