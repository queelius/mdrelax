# =============================================================================
# Visualization Functions for Simulation Studies
# =============================================================================
#
# Publication-quality plots for:
# - KL-divergence efficiency study
# - Misspecification bias study
# - Identifiability analysis
#
# Dependencies: ggplot2, gridExtra
#
# Author: Alexander Towell
# Package: md_series_system_relaxed_candidate_set_models
# =============================================================================

# Check for ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed. Plotting functions will not work.")
}

# -----------------------------------------------------------------------------
# Theme and Color Settings
# -----------------------------------------------------------------------------

#' Default theme for publication-quality plots
#' @keywords internal
theme_publication <- function(base_size = 12) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

    ggplot2::theme_bw(base_size = base_size) +
        ggplot2::theme(
            # Panel
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "black", linewidth = 0.5),

            # Axes
            axis.title = ggplot2::element_text(size = base_size),
            axis.text = ggplot2::element_text(size = base_size - 1),

            # Legend
            legend.position = "bottom",
            legend.title = ggplot2::element_text(size = base_size - 1),
            legend.text = ggplot2::element_text(size = base_size - 2),
            legend.background = ggplot2::element_rect(fill = "white", color = NA),

            # Strip (for facets)
            strip.background = ggplot2::element_rect(fill = "grey90"),
            strip.text = ggplot2::element_text(size = base_size - 1),

            # Title
            plot.title = ggplot2::element_text(size = base_size + 2, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5)
        )
}

#' Color palette for components
#' @keywords internal
component_colors <- function(m = 3) {
    # Color-blind friendly palette
    colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#000000")
    colors[seq_len(min(m, length(colors)))]
}

#' Color palette for models (correct vs misspecified)
#' @keywords internal
model_colors <- function() {
    c("correct" = "#009E73", "misspecified" = "#D55E00")
}

# -----------------------------------------------------------------------------
# KL-Divergence Efficiency Plots
# -----------------------------------------------------------------------------

#' Plot bias vs KL-divergence
#'
#' @param study_results Results from run_kl_efficiency_study
#' @param title Plot title (optional)
#'
#' @return ggplot object
#' @export
plot_bias_by_kl <- function(study_results, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$summary
    m <- max(df$component)

    if (is.null(title)) {
        title <- "Bias vs. KL-Divergence from Non-Informative Masking"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = kl_target, y = bias,
                                           color = factor(component),
                                           linetype = factor(n))) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::scale_color_manual(
            name = "Component",
            values = component_colors(m),
            labels = paste("lambda", seq_len(m), sep = "_")
        ) +
        ggplot2::scale_linetype_discrete(name = "Sample Size") +
        ggplot2::labs(
            x = "KL-Divergence (d)",
            y = "Bias",
            title = title
        ) +
        theme_publication()

    p
}

#' Plot RMSE vs KL-divergence
#'
#' @param study_results Results from run_kl_efficiency_study
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_rmse_by_kl <- function(study_results, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$summary
    m <- max(df$component)

    if (is.null(title)) {
        title <- "RMSE vs. KL-Divergence"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = kl_target, y = rmse,
                                           color = factor(component))) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::facet_wrap(~ n, labeller = ggplot2::labeller(
            n = function(x) paste("n =", x)
        )) +
        ggplot2::scale_color_manual(
            name = "Component",
            values = component_colors(m),
            labels = paste("lambda", seq_len(m), sep = "_")
        ) +
        ggplot2::labs(
            x = "KL-Divergence (d)",
            y = "RMSE",
            title = title
        ) +
        theme_publication()

    p
}

#' Plot relative efficiency vs KL-divergence
#'
#' @param study_results Results from run_kl_efficiency_study
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_efficiency_by_kl <- function(study_results, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$efficiency
    m <- max(df$component)

    if (is.null(title)) {
        title <- "Relative Efficiency vs. Non-Informative Baseline (d=0)"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = kl_target, y = rel_efficiency,
                                           color = factor(component))) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
        ggplot2::facet_wrap(~ n, labeller = ggplot2::labeller(
            n = function(x) paste("n =", x)
        )) +
        ggplot2::scale_color_manual(
            name = "Component",
            values = component_colors(m),
            labels = paste("lambda", seq_len(m), sep = "_")
        ) +
        ggplot2::labs(
            x = "KL-Divergence (d)",
            y = "Relative Efficiency (MSE_baseline / MSE)",
            title = title,
            subtitle = "Values > 1 indicate improvement over baseline"
        ) +
        theme_publication()

    p
}

#' Plot coverage probability vs KL-divergence
#'
#' @param study_results Results from run_kl_efficiency_study
#' @param target_coverage Target coverage (default 0.95)
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_coverage_by_kl <- function(study_results, target_coverage = 0.95, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$summary
    if (all(is.na(df$coverage))) {
        warning("No coverage data available")
        return(NULL)
    }

    m <- max(df$component)

    if (is.null(title)) {
        title <- sprintf("%.0f%% CI Coverage Probability", target_coverage * 100)
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = kl_target, y = coverage,
                                           color = factor(component))) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_hline(yintercept = target_coverage, linetype = "dashed",
                            color = "grey50") +
        ggplot2::facet_wrap(~ n, labeller = ggplot2::labeller(
            n = function(x) paste("n =", x)
        )) +
        ggplot2::scale_color_manual(
            name = "Component",
            values = component_colors(m),
            labels = paste("lambda", seq_len(m), sep = "_")
        ) +
        ggplot2::scale_y_continuous(limits = c(0.8, 1.0)) +
        ggplot2::labs(
            x = "KL-Divergence (d)",
            y = "Coverage Probability",
            title = title
        ) +
        theme_publication()

    p
}

#' Generate all KL efficiency plots
#'
#' @param study_results Results from run_kl_efficiency_study
#' @param output_dir Directory to save plots (NULL for no saving)
#' @param format Output format ("pdf", "png", "both")
#'
#' @return List of ggplot objects
#' @export
plot_kl_efficiency_all <- function(study_results, output_dir = NULL, format = "pdf") {
    plots <- list(
        bias = plot_bias_by_kl(study_results),
        rmse = plot_rmse_by_kl(study_results),
        efficiency = plot_efficiency_by_kl(study_results),
        coverage = plot_coverage_by_kl(study_results)
    )

    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        for (name in names(plots)) {
            if (!is.null(plots[[name]])) {
                save_plot(plots[[name]],
                          file.path(output_dir, paste0("kl_", name)),
                          format = format)
            }
        }
    }

    plots
}

# -----------------------------------------------------------------------------
# Misspecification Bias Plots
# -----------------------------------------------------------------------------

#' Plot bias comparison between correct and misspecified models
#'
#' @param study_results Results from run_misspecification_study
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_misspec_bias_comparison <- function(study_results, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$summary
    m <- max(df$component)

    if (is.null(title)) {
        title <- "Bias: Correct vs. Misspecified Model"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = alpha, y = bias,
                                           color = model,
                                           linetype = factor(component))) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::scale_color_manual(name = "Model", values = model_colors()) +
        ggplot2::scale_linetype_discrete(
            name = "Component",
            labels = paste("lambda", seq_len(m), sep = "_")
        ) +
        ggplot2::labs(
            x = "Informativeness Parameter (alpha)",
            y = "Bias",
            title = title,
            subtitle = "alpha = 0: non-informative, alpha > 0: informative masking"
        ) +
        theme_publication()

    p
}

#' Plot RMSE comparison between models
#'
#' @param study_results Results from run_misspecification_study
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_misspec_rmse_comparison <- function(study_results, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$summary

    if (is.null(title)) {
        title <- "RMSE: Correct vs. Misspecified Model"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = alpha, y = rmse,
                                           color = model,
                                           group = interaction(model, component))) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::facet_wrap(~ component, labeller = ggplot2::labeller(
            component = function(x) paste("Component", x)
        )) +
        ggplot2::scale_color_manual(name = "Model", values = model_colors()) +
        ggplot2::labs(
            x = "Informativeness Parameter (alpha)",
            y = "RMSE",
            title = title
        ) +
        theme_publication()

    p
}

#' Plot RMSE ratio (misspecified / correct)
#'
#' @param study_results Results from run_misspecification_study
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_misspec_rmse_ratio <- function(study_results, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$comparison
    m <- max(df$component)

    if (is.null(title)) {
        title <- "Misspecification Penalty (RMSE Ratio)"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = alpha, y = rmse_ratio,
                                           color = factor(component))) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
        ggplot2::scale_color_manual(
            name = "Component",
            values = component_colors(m),
            labels = paste("lambda", seq_len(m), sep = "_")
        ) +
        ggplot2::labs(
            x = "Informativeness Parameter (alpha)",
            y = "RMSE Ratio (Misspecified / Correct)",
            title = title,
            subtitle = "Values > 1 indicate misspecification penalty"
        ) +
        theme_publication()

    p
}

#' Generate all misspecification plots
#'
#' @param study_results Results from run_misspecification_study
#' @param output_dir Directory to save plots
#' @param format Output format
#'
#' @return List of ggplot objects
#' @export
plot_misspec_all <- function(study_results, output_dir = NULL, format = "pdf") {
    plots <- list(
        bias = plot_misspec_bias_comparison(study_results),
        rmse = plot_misspec_rmse_comparison(study_results),
        ratio = plot_misspec_rmse_ratio(study_results)
    )

    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        for (name in names(plots)) {
            save_plot(plots[[name]],
                      file.path(output_dir, paste0("misspec_", name)),
                      format = format)
        }
    }

    plots
}

# -----------------------------------------------------------------------------
# Identifiability Analysis Plots
# -----------------------------------------------------------------------------

#' Plot FIM eigenvalues vs correlation
#'
#' @param study_results Results from run_identifiability_study
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_fim_eigenvalues <- function(study_results, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$summary
    df_unique <- df[!duplicated(df$rho), c("rho", "smallest_eigenvalue",
                                            "condition_number", "near_singular")]

    if (is.null(title)) {
        title <- "FIM Smallest Eigenvalue vs. Correlation"
    }

    p <- ggplot2::ggplot(df_unique, ggplot2::aes(x = rho, y = smallest_eigenvalue)) +
        ggplot2::geom_line(linewidth = 0.8, color = "#0072B2") +
        ggplot2::geom_point(ggplot2::aes(color = near_singular), size = 3) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::scale_color_manual(
            name = "Near Singular",
            values = c("FALSE" = "#009E73", "TRUE" = "#D55E00")
        ) +
        ggplot2::scale_y_log10() +
        ggplot2::labs(
            x = "Correlation (rho)",
            y = "Smallest Eigenvalue (log scale)",
            title = title,
            subtitle = "Eigenvalue near 0 indicates non-identifiability"
        ) +
        theme_publication()

    p
}

#' Plot condition number vs correlation
#'
#' @param study_results Results from run_identifiability_study
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_condition_number <- function(study_results, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$summary
    df_unique <- df[!duplicated(df$rho), c("rho", "condition_number")]

    if (is.null(title)) {
        title <- "FIM Condition Number vs. Correlation"
    }

    p <- ggplot2::ggplot(df_unique, ggplot2::aes(x = rho, y = condition_number)) +
        ggplot2::geom_line(linewidth = 0.8, color = "#0072B2") +
        ggplot2::geom_point(size = 3, color = "#0072B2") +
        ggplot2::scale_y_log10() +
        ggplot2::labs(
            x = "Correlation (rho)",
            y = "Condition Number (log scale)",
            title = title,
            subtitle = "High condition number indicates numerical instability"
        ) +
        theme_publication()

    p
}

#' Plot RMSE vs correlation
#'
#' @param study_results Results from run_identifiability_study
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_rmse_by_correlation <- function(study_results, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$summary
    m <- max(df$component)

    if (is.null(title)) {
        title <- "RMSE vs. Candidate Set Correlation"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = rho, y = rmse,
                                           color = factor(component))) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_color_manual(
            name = "Component",
            values = component_colors(m),
            labels = paste("lambda", seq_len(m), sep = "_")
        ) +
        ggplot2::labs(
            x = "Correlation (rho)",
            y = "RMSE",
            title = title
        ) +
        theme_publication()

    p
}

#' Combined identifiability diagnostic plot
#'
#' @param study_results Results from run_identifiability_study
#'
#' @return Combined ggplot (requires gridExtra)
#' @export
plot_identifiability_combined <- function(study_results) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        warning("gridExtra required for combined plot")
        return(NULL)
    }

    p1 <- plot_fim_eigenvalues(study_results) +
        ggplot2::theme(legend.position = "none")
    p2 <- plot_condition_number(study_results)
    p3 <- plot_rmse_by_correlation(study_results)

    gridExtra::grid.arrange(p1, p2, p3, ncol = 2,
                             layout_matrix = rbind(c(1, 2), c(3, 3)))
}

#' Generate all identifiability plots
#'
#' @param study_results Results from run_identifiability_study
#' @param output_dir Directory to save plots
#' @param format Output format
#'
#' @return List of ggplot objects
#' @export
plot_identifiability_all <- function(study_results, output_dir = NULL, format = "pdf") {
    plots <- list(
        eigenvalues = plot_fim_eigenvalues(study_results),
        condition = plot_condition_number(study_results),
        rmse = plot_rmse_by_correlation(study_results)
    )

    # Try to create combined plot
    tryCatch({
        plots$combined <- plot_identifiability_combined(study_results)
    }, error = function(e) {
        warning("Could not create combined plot: ", e$message)
    })

    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        for (name in names(plots)) {
            if (!is.null(plots[[name]])) {
                save_plot(plots[[name]],
                          file.path(output_dir, paste0("ident_", name)),
                          format = format,
                          width = if (name == "combined") 10 else 6,
                          height = if (name == "combined") 8 else 4)
            }
        }
    }

    plots
}

# -----------------------------------------------------------------------------
# MSE Comparison Heatmap
# -----------------------------------------------------------------------------

#' Plot MSE comparison heatmap
#'
#' Creates a heatmap showing MSE across scenarios (e.g., sample size vs KL-divergence)
#'
#' @param study_results Results from run_kl_efficiency_study
#' @param component Which component to plot (default 1)
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_mse_heatmap <- function(study_results, component = 1, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    df <- study_results$summary
    df <- df[df$component == component, ]

    if (is.null(title)) {
        title <- sprintf("MSE Heatmap (Component %d)", component)
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(kl_target),
                                           y = factor(n),
                                           fill = mse)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4f", mse)),
                           color = "white", size = 3) +
        ggplot2::scale_fill_gradient(low = "#009E73", high = "#D55E00",
                                      name = "MSE") +
        ggplot2::labs(
            x = "KL-Divergence (d)",
            y = "Sample Size (n)",
            title = title
        ) +
        theme_publication()

    p
}

# -----------------------------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------------------------

#' Save plot to file
#'
#' @param plot ggplot object
#' @param filename Base filename (without extension)
#' @param format "pdf", "png", or "both"
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi DPI for PNG
#'
#' @keywords internal
save_plot <- function(plot, filename, format = "pdf",
                       width = 6, height = 4, dpi = 300) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

    if (format %in% c("pdf", "both")) {
        ggplot2::ggsave(paste0(filename, ".pdf"), plot,
                        width = width, height = height)
    }

    if (format %in% c("png", "both")) {
        ggplot2::ggsave(paste0(filename, ".png"), plot,
                        width = width, height = height, dpi = dpi)
    }

    invisible(filename)
}

#' Create sampling distribution plot
#'
#' Shows the distribution of MLE estimates across replications
#'
#' @param estimates Matrix of estimates (B x m)
#' @param theta_true True parameter values
#' @param component Which component to plot (NULL for all)
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_sampling_distribution <- function(estimates, theta_true,
                                        component = NULL, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    # Convert to long format
    m <- ncol(estimates)
    df <- data.frame(
        component = rep(seq_len(m), each = nrow(estimates)),
        estimate = as.vector(estimates)
    )
    df <- df[!is.na(df$estimate), ]

    # True values
    truth_df <- data.frame(
        component = seq_len(m),
        theta_true = theta_true
    )

    if (!is.null(component)) {
        df <- df[df$component == component, ]
        truth_df <- truth_df[truth_df$component == component, ]
    }

    if (is.null(title)) {
        title <- "Sampling Distribution of MLE"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = estimate)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                 bins = 30, fill = "steelblue", alpha = 0.7) +
        ggplot2::geom_density(color = "darkblue", linewidth = 0.8) +
        ggplot2::geom_vline(data = truth_df,
                            ggplot2::aes(xintercept = theta_true),
                            color = "red", linetype = "dashed", linewidth = 1) +
        ggplot2::facet_wrap(~ component, scales = "free",
                            labeller = ggplot2::labeller(
                                component = function(x) paste("Component", x)
                            )) +
        ggplot2::labs(
            x = "Estimate",
            y = "Density",
            title = title,
            subtitle = "Red dashed line shows true parameter value"
        ) +
        theme_publication()

    p
}

#' Create QQ plot for normality assessment
#'
#' @param estimates Vector of estimates
#' @param theta_true True parameter value
#' @param title Plot title
#'
#' @return ggplot object
#' @export
plot_qq_normal <- function(estimates, theta_true, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
    }

    estimates <- estimates[!is.na(estimates)]
    n <- length(estimates)

    # Standardize
    z <- (estimates - mean(estimates)) / stats::sd(estimates)

    df <- data.frame(
        theoretical = stats::qnorm(ppoints(n)),
        sample = sort(z)
    )

    if (is.null(title)) {
        title <- "Q-Q Plot (Normality Assessment)"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = theoretical, y = sample)) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::geom_abline(intercept = 0, slope = 1,
                              color = "red", linetype = "dashed") +
        ggplot2::labs(
            x = "Theoretical Quantiles",
            y = "Sample Quantiles",
            title = title
        ) +
        theme_publication()

    p
}
