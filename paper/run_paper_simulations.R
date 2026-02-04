#!/usr/bin/env Rscript
# =============================================================================
# Paper Simulations: Studies 4 and 5
# =============================================================================
#
# Study 4: C3 Misspecification Bias (Exponential)
#   - Scenario 5:  C1-C2-C3 data -> Relaxed C3 model (overfit check)
#   - Scenario 6:  Relaxed C3 data -> C1-C2-C3 model (misspecified)
#   - Scenario 6b: Relaxed C3 data -> Relaxed C3 model (correct, known alpha)
#
# Study 5: Weibull Series Systems
#   - W1: C1-C2-C3 baseline
#   - W3: Relaxed C2 data -> C1-C2-C3 (C2 misspecification)
#   - W4: Relaxed C2 data -> Relaxed C2 (correct)
#   - W6: Relaxed C3 data -> C1-C2-C3 (C3 misspecification)
#   - W7: Relaxed C3 data -> Relaxed C3 (correct)
#
# Usage:
#   Rscript paper/run_paper_simulations.R [--quick] [--study=4|5|all]
#
# =============================================================================

library(ggplot2)

# Load the package
devtools::load_all(quiet = TRUE)

# -----------------------------------------------------------------------------
# Parse Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

study_arg <- grep("^--study=", args, value = TRUE)
run_study <- if (length(study_arg) > 0) sub("^--study=", "", study_arg) else "all"

# Output directory — detect script location robustly
script_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error = function(e) {
        # Fallback: assume run from package root
        normalizePath("paper")
    }
)
data_dir <- file.path(script_dir, "data")
pkg_root <- dirname(script_dir)
fig_dir <- file.path(pkg_root, "inst", "simulations", "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("=============================================================\n")
cat("Paper Simulations: Studies 4 and 5\n")
cat("=============================================================\n")
cat(sprintf("Quick mode: %s\n", quick))
cat(sprintf("Study: %s\n", run_study))
cat(sprintf("Data dir: %s\n", data_dir))
cat(sprintf("Figure dir: %s\n", fig_dir))
cat("=============================================================\n\n")

# =============================================================================
# Study 4: C3 Misspecification (Exponential)
# =============================================================================

run_study_4 <- function(quick = FALSE) {
    cat("\n=== Study 4: C3 Misspecification Bias ===\n")

    # Parameters matching existing Studies 1-3
    m <- 3
    theta <- c(1, 1.5, 2)
    tau <- 3
    n <- 200
    B <- if (quick) 50 else 200
    base_p <- 0.5
    alpha_levels <- c(0, 0.5, 1, 2)

    config <- sim_config(n_sim = B, n_obs = n, theta = theta, tau = tau, seed = 42)

    results <- list()

    # Scenario 5: C1-C2-C3 data -> Relaxed C3 model (overfit check)
    cat("  Scenario 5: C1-C2-C3 data -> Relaxed C3 model\n")
    for (alpha in alpha_levels) {
        cat(sprintf("    alpha = %.1f: ", alpha))
        res5 <- run_simulation_scenario(config, scenario = 5,
                                         p = base_p, assumed_alpha = alpha,
                                         verbose = FALSE)
        res5$alpha <- alpha
        key <- paste0("S5_alpha_", alpha)
        results[[key]] <- res5
        cat(sprintf("%d/%d converged\n",
                    sum(res5$converged), nrow(res5)))
    }

    # Scenario 6: Relaxed C3 data -> C1-C2-C3 model (misspecified)
    cat("  Scenario 6: Relaxed C3 data -> C1-C2-C3 model\n")
    for (alpha in alpha_levels) {
        cat(sprintf("    alpha = %.1f: ", alpha))
        res6 <- run_simulation_scenario(config, scenario = 6,
                                         alpha = alpha, base_p = base_p,
                                         verbose = FALSE)
        res6$alpha <- alpha
        key <- paste0("S6_alpha_", alpha)
        results[[key]] <- res6
        cat(sprintf("%d/%d converged\n",
                    sum(res6$converged), nrow(res6)))
    }

    # Scenario 6b: Relaxed C3 data -> Relaxed C3 model (known alpha)
    cat("  Scenario 6b: Relaxed C3 data -> Relaxed C3 model (known alpha)\n")
    for (alpha in alpha_levels) {
        cat(sprintf("    alpha = %.1f: ", alpha))
        res6b <- run_simulation_scenario(config, scenario = "6b",
                                          alpha = alpha, base_p = base_p,
                                          verbose = FALSE)
        res6b$alpha <- alpha
        key <- paste0("S6b_alpha_", alpha)
        results[[key]] <- res6b
        cat(sprintf("%d/%d converged\n",
                    sum(res6b$converged), nrow(res6b)))
    }

    # Combine into summary
    all_res <- do.call(rbind, results)

    # Build comparison table
    comparison <- data.frame()
    for (alpha in alpha_levels) {
        # Scenario 6 (misspecified)
        res6 <- results[[paste0("S6_alpha_", alpha)]]
        res6_conv <- res6[res6$converged, ]

        # Scenario 6b (correct)
        res6b <- results[[paste0("S6b_alpha_", alpha)]]
        res6b_conv <- res6b[res6b$converged, ]

        if (nrow(res6_conv) > 0 && nrow(res6b_conv) > 0) {
            theta_est_6 <- do.call(rbind, res6_conv$theta_est)
            theta_est_6b <- do.call(rbind, res6b_conv$theta_est)

            for (j in seq_len(m)) {
                bias_misspec <- mean(theta_est_6[, j]) - theta[j]
                bias_correct <- mean(theta_est_6b[, j]) - theta[j]
                rmse_misspec <- sqrt(mean((theta_est_6[, j] - theta[j])^2))
                rmse_correct <- sqrt(mean((theta_est_6b[, j] - theta[j])^2))

                comparison <- rbind(comparison, data.frame(
                    alpha = alpha,
                    component = j,
                    bias_misspec = bias_misspec,
                    bias_correct = bias_correct,
                    rmse_misspec = rmse_misspec,
                    rmse_correct = rmse_correct,
                    rmse_ratio = rmse_misspec / rmse_correct
                ))
            }
        }
    }

    # Scenario 5 summary (overfit check)
    overfit <- data.frame()
    for (alpha in alpha_levels) {
        res5 <- results[[paste0("S5_alpha_", alpha)]]
        res5_conv <- res5[res5$converged, ]
        if (nrow(res5_conv) > 0) {
            theta_est_5 <- do.call(rbind, res5_conv$theta_est)
            for (j in seq_len(m)) {
                overfit <- rbind(overfit, data.frame(
                    alpha = alpha,
                    component = j,
                    bias = mean(theta_est_5[, j]) - theta[j],
                    rmse = sqrt(mean((theta_est_5[, j] - theta[j])^2)),
                    converged = nrow(res5_conv),
                    total = nrow(res5)
                ))
            }
        }
    }

    study4 <- list(
        params = list(m = m, theta = theta, tau = tau, n = n, B = B,
                      base_p = base_p, alpha_levels = alpha_levels),
        results = all_res,
        comparison = comparison,
        overfit = overfit
    )

    saveRDS(study4, file.path(data_dir, "study4_c3_misspec.rds"))
    cat("\n  Saved: study4_c3_misspec.rds\n")

    study4
}

# =============================================================================
# Study 5: Weibull Series Systems
# =============================================================================

run_study_5 <- function(quick = FALSE) {
    cat("\n=== Study 5: Weibull Series Systems ===\n")

    m <- 2
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    tau <- 8
    n <- 200
    B <- if (quick) 30 else 100
    p <- 0.3
    alpha <- 1
    base_p <- 0.5

    # P matrix for relaxed C2 scenarios
    P <- make_P_matrix(m, "full", values = c(0.3, 0.5))

    study5 <- run_weibull_simulation_study(
        n_sim = B, n_obs = n,
        shapes = shapes, scales = scales,
        tau = tau,
        scenarios = c("W1", "W3", "W4", "W6", "W7"),
        p = p, P = P,
        alpha = alpha, base_p = base_p,
        seed = 42, verbose = TRUE
    )

    study5$params$p <- p
    study5$params$P <- P
    study5$params$alpha <- alpha
    study5$params$base_p <- base_p

    saveRDS(study5, file.path(data_dir, "study5_weibull.rds"))
    cat("\n  Saved: study5_weibull.rds\n")

    study5
}

# =============================================================================
# Figure Generation
# =============================================================================

generate_study4_figures <- function(study4) {
    cat("\n=== Generating Study 4 Figures ===\n")

    df <- study4$comparison
    df$component <- factor(df$component,
                           labels = c("lambda[1]", "lambda[2]", "lambda[3]"))

    # Figure 7: C3 Misspecification Bias Comparison
    df_long <- reshape(df[, c("alpha", "component", "bias_misspec", "bias_correct")],
                       direction = "long",
                       varying = list(c("bias_misspec", "bias_correct")),
                       v.names = "bias",
                       timevar = "model",
                       times = c("Misspecified (C3)", "Correct"),
                       idvar = c("alpha", "component"))

    p7 <- ggplot(df_long, aes(x = alpha, y = bias, color = model, linetype = model)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
        geom_line(linewidth = 1) +
        geom_point(size = 2.5) +
        facet_wrap(~component, labeller = label_parsed) +
        scale_color_manual(values = c("Misspecified (C3)" = "#E41A1C",
                                      "Correct" = "#377EB8")) +
        labs(x = expression(paste("Power Parameter (", alpha, ")")),
             y = "Bias",
             title = "C3 Misspecification: Bias Comparison",
             subtitle = "Misspecified model ignores parameter-dependent masking",
             color = "Model", linetype = "Model") +
        theme_minimal(base_size = 12) +
        theme(
            legend.position = "bottom",
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            strip.text = element_text(size = 12)
        )

    ggsave(file.path(fig_dir, "fig7_c3_misspec_bias.pdf"), p7, width = 9, height = 5)
    ggsave(file.path(fig_dir, "fig7_c3_misspec_bias.png"), p7, width = 9, height = 5, dpi = 300)
    cat("  Saved: fig7_c3_misspec_bias.pdf/png\n")

    # Figure 8: RMSE Ratio by alpha
    p8 <- ggplot(df, aes(x = alpha, y = rmse_ratio, color = component)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                           labels = c(expression(lambda[1]),
                                      expression(lambda[2]),
                                      expression(lambda[3]))) +
        labs(x = expression(paste("Power Parameter (", alpha, ")")),
             y = "RMSE Ratio (Misspecified / Correct)",
             title = "C3 Misspecification: Relative Efficiency Loss",
             subtitle = "Ratio > 1 indicates efficiency loss from ignoring C3 violation",
             color = "Parameter") +
        theme_minimal(base_size = 12) +
        theme(
            legend.position = "right",
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold")
        )

    ggsave(file.path(fig_dir, "fig8_c3_rmse_ratio.pdf"), p8, width = 7, height = 5)
    ggsave(file.path(fig_dir, "fig8_c3_rmse_ratio.png"), p8, width = 7, height = 5, dpi = 300)
    cat("  Saved: fig8_c3_rmse_ratio.pdf/png\n")
}

generate_study5_figures <- function(study5) {
    cat("\n=== Generating Study 5 Figures ===\n")

    df <- study5$summary
    theta_true <- study5$config$theta
    m <- study5$config$m

    # Create parameter labels for Weibull: k1, lambda1, k2, lambda2
    param_labels <- character(2 * m)
    for (j in seq_len(m)) {
        param_labels[2*j - 1] <- paste0("k[", j, "]")
        param_labels[2*j]     <- paste0("lambda[", j, "]")
    }
    df$param_label <- factor(param_labels[df$component],
                             levels = param_labels)

    # Figure 9: Weibull Baseline MLE Performance (W1)
    df_w1 <- df[df$scenario == "W1", ]

    p9 <- ggplot(df_w1, aes(x = param_label, y = bias)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_col(fill = "#377EB8", alpha = 0.8, width = 0.6) +
        geom_errorbar(aes(ymin = bias - rmse, ymax = bias + rmse),
                      width = 0.2, color = "#333333") +
        scale_x_discrete(labels = function(x) parse(text = x)) +
        labs(x = "Parameter", y = "Bias (bars show RMSE range)",
             title = "Weibull Baseline: MLE Performance (W1)",
             subtitle = sprintf("n=%d, B=%d, shapes=(%.1f,%.1f), scales=(%.1f,%.1f)",
                                study5$config$n_obs, study5$config$n_sim,
                                study5$config$shapes[1], study5$config$shapes[2],
                                study5$config$scales[1], study5$config$scales[2])) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold")
        )

    ggsave(file.path(fig_dir, "fig9_weibull_baseline.pdf"), p9, width = 7, height = 5)
    ggsave(file.path(fig_dir, "fig9_weibull_baseline.png"), p9, width = 7, height = 5, dpi = 300)
    cat("  Saved: fig9_weibull_baseline.pdf/png\n")

    # Figure 10: Weibull Misspecification Comparison
    # Compare W3 vs W4 (C2) and W6 vs W7 (C3)
    df_misspec <- df[df$scenario %in% c("W3", "W4", "W6", "W7"), ]

    # Label scenarios
    scenario_labels <- c(
        "W3" = "C2 Misspecified",
        "W4" = "C2 Correct",
        "W6" = "C3 Misspecified",
        "W7" = "C3 Correct"
    )
    df_misspec$scenario_label <- factor(scenario_labels[df_misspec$scenario],
                                        levels = c("C2 Misspecified", "C2 Correct",
                                                   "C3 Misspecified", "C3 Correct"))
    df_misspec$violation <- ifelse(df_misspec$scenario %in% c("W3", "W4"),
                                   "C2 Violation", "C3 Violation")
    df_misspec$model_type <- ifelse(df_misspec$scenario %in% c("W3", "W6"),
                                     "Misspecified", "Correct")

    p10 <- ggplot(df_misspec, aes(x = param_label, y = bias, fill = model_type)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
        facet_wrap(~violation) +
        scale_x_discrete(labels = function(x) parse(text = x)) +
        scale_fill_manual(values = c("Misspecified" = "#E41A1C",
                                     "Correct" = "#377EB8")) +
        labs(x = "Parameter", y = "Bias",
             title = "Weibull: Misspecification Impact",
             subtitle = "Comparing correct vs misspecified models under C2 and C3 violations",
             fill = "Model") +
        theme_minimal(base_size = 12) +
        theme(
            legend.position = "bottom",
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            strip.text = element_text(size = 12)
        )

    ggsave(file.path(fig_dir, "fig10_weibull_misspec.pdf"), p10, width = 9, height = 5)
    ggsave(file.path(fig_dir, "fig10_weibull_misspec.png"), p10, width = 9, height = 5, dpi = 300)
    cat("  Saved: fig10_weibull_misspec.pdf/png\n")
}

# =============================================================================
# Main Execution
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

if (run_study %in% c("all", "4")) {
    study4 <- run_study_4(quick = quick)
    generate_study4_figures(study4)
} else {
    # Try to load existing results for figure generation
    s4_file <- file.path(data_dir, "study4_c3_misspec.rds")
    if (file.exists(s4_file)) study4 <- readRDS(s4_file)
}

if (run_study %in% c("all", "5")) {
    study5 <- run_study_5(quick = quick)
    generate_study5_figures(study5)
} else {
    s5_file <- file.path(data_dir, "study5_weibull.rds")
    if (file.exists(s5_file)) study5 <- readRDS(s5_file)
}

cat("\n=============================================================\n")
cat("All paper simulations complete.\n")
cat("=============================================================\n")
