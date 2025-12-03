# =============================================================================
# GENERATE PUBLICATION FIGURES
# =============================================================================
# Run from the project root directory:
#   Rscript inst/simulations/generate_figures.R

# Set working directory to project root
setwd("/home/spinoza/github/rlang/md-series-systems-relaxed-candidate-set-models")

# Load utilities
source("inst/simulations/sim_utils.R")
source("inst/simulations/sim_plots.R")

# Load results
results_kl <- readRDS("inst/simulations/results/kl_efficiency_results.rds")
results_misspec <- readRDS("inst/simulations/results/misspecification_results.rds")
results_ident <- readRDS("inst/simulations/results/identifiability_results.rds")

cat("\n")
cat("==============================================================\n")
cat("  GENERATING PUBLICATION FIGURES\n")
cat("==============================================================\n\n")

library(ggplot2)
library(dplyr)
library(tidyr)

# Set output directory
fig_dir <- "inst/simulations/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Figure 1: Sample Size Effect on RMSE (from KL study)
# -----------------------------------------------------------------------------
cat("Creating Figure 1: RMSE by Sample Size...\n")

df_kl <- results_kl$summary
df_kl$component <- factor(df_kl$component, labels = c("lambda[1]", "lambda[2]", "lambda[3]"))

p1 <- ggplot(df_kl, aes(x = n, y = rmse, color = component, group = component)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     labels = c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3]))) +
  labs(x = "Sample Size (n)", y = "RMSE",
       title = "MLE Performance by Sample Size",
       subtitle = "Exponential series system with Bernoulli masking (p=0.3)",
       color = "Parameter") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig1_rmse_by_sample_size.pdf"), p1, width = 7, height = 5)
ggsave(file.path(fig_dir, "fig1_rmse_by_sample_size.png"), p1, width = 7, height = 5, dpi = 300)
cat("  Saved: fig1_rmse_by_sample_size.pdf/png\n")

# -----------------------------------------------------------------------------
# Figure 2: Coverage Probability by Sample Size
# -----------------------------------------------------------------------------
cat("Creating Figure 2: Coverage by Sample Size...\n")

p2 <- ggplot(df_kl, aes(x = n, y = coverage, color = component, group = component)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     labels = c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3]))) +
  scale_y_continuous(limits = c(0.85, 1.0)) +
  labs(x = "Sample Size (n)", y = "Coverage Probability",
       title = "95% CI Coverage Probability",
       subtitle = "Dashed line indicates nominal 95% level",
       color = "Parameter") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig2_coverage_by_sample_size.pdf"), p2, width = 7, height = 5)
ggsave(file.path(fig_dir, "fig2_coverage_by_sample_size.png"), p2, width = 7, height = 5, dpi = 300)
cat("  Saved: fig2_coverage_by_sample_size.pdf/png\n")

# -----------------------------------------------------------------------------
# Figure 3: Misspecification Bias Comparison
# -----------------------------------------------------------------------------
cat("Creating Figure 3: Misspecification Bias...\n")

df_comp <- results_misspec$comparison
df_comp$component <- factor(df_comp$component, labels = c("lambda[1]", "lambda[2]", "lambda[3]"))

# Reshape for plotting
df_long <- df_comp %>%
  select(alpha, component, bias_misspec, bias_correct) %>%
  pivot_longer(cols = c(bias_misspec, bias_correct),
               names_to = "model", values_to = "bias") %>%
  mutate(model = ifelse(model == "bias_misspec", "Misspecified (C2)", "Correct"))

p3 <- ggplot(df_long, aes(x = alpha, y = bias, color = model, linetype = model)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~component, labeller = label_parsed) +
  scale_color_manual(values = c("Misspecified (C2)" = "#E41A1C", "Correct" = "#377EB8")) +
  labs(x = expression(paste("Informativeness Parameter (", alpha, ")")),
       y = "Bias",
       title = "Bias Under Correct vs Misspecified Model",
       subtitle = "Misspecified model incorrectly assumes C2 (non-informative masking)",
       color = "Model", linetype = "Model") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 12)
  )

ggsave(file.path(fig_dir, "fig3_misspecification_bias.pdf"), p3, width = 9, height = 5)
ggsave(file.path(fig_dir, "fig3_misspecification_bias.png"), p3, width = 9, height = 5, dpi = 300)
cat("  Saved: fig3_misspecification_bias.pdf/png\n")

# -----------------------------------------------------------------------------
# Figure 4: RMSE Ratio (Misspecified / Correct)
# -----------------------------------------------------------------------------
cat("Creating Figure 4: RMSE Ratio...\n")

p4 <- ggplot(df_comp, aes(x = alpha, y = rmse_ratio, color = component)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     labels = c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3]))) +
  scale_y_continuous(limits = c(0.9, 1.15)) +
  labs(x = expression(paste("Informativeness Parameter (", alpha, ")")),
       y = "RMSE Ratio (Misspecified / Correct)",
       title = "Relative Efficiency Loss from Misspecification",
       subtitle = "Ratio > 1 indicates efficiency loss from assuming C2 when violated",
       color = "Parameter") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig4_rmse_ratio.pdf"), p4, width = 7, height = 5)
ggsave(file.path(fig_dir, "fig4_rmse_ratio.png"), p4, width = 7, height = 5, dpi = 300)
cat("  Saved: fig4_rmse_ratio.pdf/png\n")

# -----------------------------------------------------------------------------
# Figure 5: FIM Eigenvalues by Correlation
# -----------------------------------------------------------------------------
cat("Creating Figure 5: FIM Eigenvalues...\n")

df_ident <- results_ident$summary
df_ident_unique <- unique(df_ident[, c("rho", "smallest_eigenvalue", "condition_number")])

p5 <- ggplot(df_ident_unique, aes(x = rho, y = smallest_eigenvalue)) +
  geom_line(linewidth = 1, color = "#377EB8") +
  geom_point(size = 3, color = "#377EB8") +
  labs(x = expression(paste("Candidate Set Correlation (", rho, ")")),
       y = "Smallest FIM Eigenvalue",
       title = "Identifiability: FIM Smallest Eigenvalue",
       subtitle = "Eigenvalue approaching 0 indicates non-identifiability") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig5_fim_eigenvalue.pdf"), p5, width = 7, height = 5)
ggsave(file.path(fig_dir, "fig5_fim_eigenvalue.png"), p5, width = 7, height = 5, dpi = 300)
cat("  Saved: fig5_fim_eigenvalue.pdf/png\n")

# -----------------------------------------------------------------------------
# Figure 6: RMSE by Correlation
# -----------------------------------------------------------------------------
cat("Creating Figure 6: RMSE by Correlation...\n")

df_ident$component <- factor(df_ident$component, labels = c("lambda[1]", "lambda[2]", "lambda[3]"))

p6 <- ggplot(df_ident, aes(x = rho, y = rmse, color = component)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     labels = c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3]))) +
  labs(x = expression(paste("Candidate Set Correlation (", rho, ")")),
       y = "RMSE",
       title = "MLE Performance vs Candidate Set Correlation",
       color = "Parameter") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig6_rmse_by_correlation.pdf"), p6, width = 7, height = 5)
ggsave(file.path(fig_dir, "fig6_rmse_by_correlation.png"), p6, width = 7, height = 5, dpi = 300)
cat("  Saved: fig6_rmse_by_correlation.pdf/png\n")

cat("\n")
cat("==============================================================\n")
cat("  ALL FIGURES GENERATED SUCCESSFULLY\n")
cat("==============================================================\n")
cat("\nFigures saved to:", fig_dir, "\n")
