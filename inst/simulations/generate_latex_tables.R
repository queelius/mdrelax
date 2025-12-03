# =============================================================================
# GENERATE LATEX TABLES FROM SIMULATION RESULTS
# =============================================================================
# Run from project root: Rscript inst/simulations/generate_latex_tables.R

# Set working directory
setwd("/home/spinoza/github/rlang/md-series-systems-relaxed-candidate-set-models")

# Load utilities
source("inst/simulations/sim_utils.R")

# Load results
results_kl <- readRDS("inst/simulations/results/kl_efficiency_results.rds")
results_misspec <- readRDS("inst/simulations/results/misspecification_results.rds")
results_ident <- readRDS("inst/simulations/results/identifiability_results.rds")

# Create output directory
tables_dir <- "inst/simulations/tables"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("==============================================================\n")
cat("  GENERATING LATEX TABLES\n")
cat("==============================================================\n\n")

# -----------------------------------------------------------------------------
# Table 1: MLE Performance Summary by Sample Size
# -----------------------------------------------------------------------------
cat("Creating Table 1: MLE Performance by Sample Size...\n")

df <- results_kl$summary
df <- df[order(df$n, df$component), ]

# Build LaTeX table
latex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Maximum Likelihood Estimation Performance by Sample Size}",
  "\\label{tab:mle_performance}",
  "\\begin{tabular}{cccccc}",
  "\\toprule",
  "$n$ & Parameter & Bias & RMSE & Coverage & Mean CI Width \\\\",
  "\\midrule"
)

sample_sizes <- sort(unique(df$n))
for (n in sample_sizes) {
  sub <- df[df$n == n, ]
  first_n <- TRUE
  for (i in 1:nrow(sub)) {
    row <- sub[i, ]
    n_str <- if (first_n) as.character(n) else ""
    param <- sprintf("$\\lambda_%d$", row$component)
    latex <- c(latex, sprintf("%s & %s & %.4f & %.4f & %.3f & %.3f \\\\",
                              n_str, param, row$bias, row$rmse,
                              row$coverage, row$mean_ci_width))
    first_n <- FALSE
  }
  if (n != tail(sample_sizes, 1)) {
    latex <- c(latex, "\\midrule")
  }
}

latex <- c(latex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}",
  "\\small",
  "\\item Note: Results based on 200 Monte Carlo replications. True parameters:",
  "\\item $\\lambda_1 = 1.0$, $\\lambda_2 = 1.5$, $\\lambda_3 = 2.0$.",
  "\\item Bernoulli masking with $p = 0.3$, censoring proportion $\\approx 20\\%$.",
  "\\end{tablenotes}",
  "\\end{table}"
)

writeLines(latex, file.path(tables_dir, "table1_mle_performance.tex"))
cat("  Saved: table1_mle_performance.tex\n")

# Print to console
cat("\n--- Table 1: MLE Performance ---\n")
cat(paste(latex, collapse = "\n"))
cat("\n\n")

# -----------------------------------------------------------------------------
# Table 2: Misspecification Bias Comparison
# -----------------------------------------------------------------------------
cat("Creating Table 2: Misspecification Bias...\n")

df_comp <- results_misspec$comparison
df_comp <- df_comp[order(df_comp$alpha, df_comp$component), ]

latex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Bias Comparison: Correct vs Misspecified Model}",
  "\\label{tab:misspecification_bias}",
  "\\begin{tabular}{ccccc}",
  "\\toprule",
  "$\\alpha$ & Parameter & Bias (Correct) & Bias (Misspec.) & RMSE Ratio \\\\",
  "\\midrule"
)

alpha_vals <- sort(unique(df_comp$alpha))
for (alpha in alpha_vals) {
  sub <- df_comp[df_comp$alpha == alpha, ]
  first_alpha <- TRUE
  for (i in 1:nrow(sub)) {
    row <- sub[i, ]
    alpha_str <- if (first_alpha) sprintf("%.0f", alpha) else ""
    param <- sprintf("$\\lambda_%d$", row$component)
    latex <- c(latex, sprintf("%s & %s & %.4f & %.4f & %.3f \\\\",
                              alpha_str, param, row$bias_correct, row$bias_misspec,
                              row$rmse_ratio))
    first_alpha <- FALSE
  }
  if (alpha != tail(alpha_vals, 1)) {
    latex <- c(latex, "\\midrule")
  }
}

latex <- c(latex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}",
  "\\small",
  "\\item Note: $\\alpha = 0$ corresponds to non-informative masking (C2 satisfied).",
  "\\item As $\\alpha$ increases, masking becomes more informative.",
  "\\item RMSE Ratio = RMSE(Misspecified) / RMSE(Correct); values $> 1$ indicate efficiency loss.",
  "\\end{tablenotes}",
  "\\end{table}"
)

writeLines(latex, file.path(tables_dir, "table2_misspecification.tex"))
cat("  Saved: table2_misspecification.tex\n")

# Print to console
cat("\n--- Table 2: Misspecification Bias ---\n")
cat(paste(latex, collapse = "\n"))
cat("\n\n")

# -----------------------------------------------------------------------------
# Table 3: Identifiability Analysis
# -----------------------------------------------------------------------------
cat("Creating Table 3: Identifiability Analysis...\n")

df_ident <- results_ident$summary
df_ident_unique <- unique(df_ident[, c("rho", "smallest_eigenvalue", "condition_number")])
df_ident_unique <- df_ident_unique[order(df_ident_unique$rho), ]

latex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Fisher Information Matrix Analysis by Candidate Set Correlation}",
  "\\label{tab:identifiability}",
  "\\begin{tabular}{ccc}",
  "\\toprule",
  "$\\rho$ & Smallest Eigenvalue & Condition Number \\\\",
  "\\midrule"
)

for (i in 1:nrow(df_ident_unique)) {
  row <- df_ident_unique[i, ]
  latex <- c(latex, sprintf("%.1f & %.4f & %.2f \\\\",
                            row$rho, row$smallest_eigenvalue, row$condition_number))
}

latex <- c(latex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}",
  "\\small",
  "\\item Note: $\\rho$ measures correlation between candidate set indicators.",
  "\\item As $\\rho \\to 1$, components always co-occur in candidate sets,",
  "\\item leading to non-identifiability (eigenvalue $\\to 0$, condition number $\\to \\infty$).",
  "\\end{tablenotes}",
  "\\end{table}"
)

writeLines(latex, file.path(tables_dir, "table3_identifiability.tex"))
cat("  Saved: table3_identifiability.tex\n")

# Print to console
cat("\n--- Table 3: Identifiability ---\n")
cat(paste(latex, collapse = "\n"))
cat("\n\n")

# -----------------------------------------------------------------------------
# Table 4: RMSE by Correlation
# -----------------------------------------------------------------------------
cat("Creating Table 4: RMSE by Correlation...\n")

df_rmse <- df_ident[order(df_ident$rho, df_ident$component), ]

latex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{MLE Performance by Candidate Set Correlation}",
  "\\label{tab:rmse_correlation}",
  "\\begin{tabular}{ccccc}",
  "\\toprule",
  "$\\rho$ & Parameter & Bias & RMSE & Relative Bias \\\\",
  "\\midrule"
)

rho_vals <- sort(unique(df_rmse$rho))
for (rho in rho_vals) {
  sub <- df_rmse[df_rmse$rho == rho, ]
  first_rho <- TRUE
  for (i in 1:nrow(sub)) {
    row <- sub[i, ]
    rho_str <- if (first_rho) sprintf("%.1f", rho) else ""
    param <- sprintf("$\\lambda_%d$", row$component)
    latex <- c(latex, sprintf("%s & %s & %.4f & %.4f & %.4f \\\\",
                              rho_str, param, row$bias, row$rmse, row$rel_bias))
    first_rho <- FALSE
  }
  if (rho != tail(rho_vals, 1)) {
    latex <- c(latex, "\\midrule")
  }
}

latex <- c(latex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}",
  "\\small",
  "\\item Note: Relative bias = Bias / True parameter value.",
  "\\end{tablenotes}",
  "\\end{table}"
)

writeLines(latex, file.path(tables_dir, "table4_rmse_correlation.tex"))
cat("  Saved: table4_rmse_correlation.tex\n")

# Print to console
cat("\n--- Table 4: RMSE by Correlation ---\n")
cat(paste(latex, collapse = "\n"))
cat("\n\n")

# -----------------------------------------------------------------------------
# Summary statistics for paper text
# -----------------------------------------------------------------------------
cat("Creating summary statistics for paper narrative...\n")

summary_stats <- list()

# KL efficiency study
summary_stats$kl <- list(
  n_values = unique(df$n),
  B = results_kl$params$B,
  rates = results_kl$params$rates,
  coverage_range = range(df$coverage),
  rmse_range = range(df$rmse),
  bias_range = range(df$bias)
)

# Misspecification study
summary_stats$misspec <- list(
  alpha_values = unique(df_comp$alpha),
  B = results_misspec$params$B,
  max_rmse_ratio = max(df_comp$rmse_ratio),
  bias_at_alpha_10 = df_comp$bias_misspec[df_comp$alpha == 10]
)

# Identifiability study
summary_stats$ident <- list(
  rho_values = unique(df_ident_unique$rho),
  B = results_ident$params$B,
  eigenvalue_at_rho_0 = df_ident_unique$smallest_eigenvalue[df_ident_unique$rho == 0],
  eigenvalue_at_rho_09 = df_ident_unique$smallest_eigenvalue[df_ident_unique$rho == 0.9],
  condition_range = range(df_ident_unique$condition_number)
)

saveRDS(summary_stats, file.path(tables_dir, "summary_statistics.rds"))
cat("  Saved: summary_statistics.rds\n")

# Print summary
cat("\n--- Summary Statistics ---\n")
cat(sprintf("KL Efficiency Study: n = %s, B = %d\n",
            paste(summary_stats$kl$n_values, collapse = ", "),
            summary_stats$kl$B))
cat(sprintf("  Coverage range: [%.3f, %.3f]\n",
            summary_stats$kl$coverage_range[1], summary_stats$kl$coverage_range[2]))
cat(sprintf("  RMSE range: [%.4f, %.4f]\n",
            summary_stats$kl$rmse_range[1], summary_stats$kl$rmse_range[2]))

cat(sprintf("\nMisspecification Study: alpha = %s, B = %d\n",
            paste(summary_stats$misspec$alpha_values, collapse = ", "),
            summary_stats$misspec$B))
cat(sprintf("  Max RMSE ratio: %.3f\n", summary_stats$misspec$max_rmse_ratio))

cat(sprintf("\nIdentifiability Study: rho = %s, B = %d\n",
            paste(round(summary_stats$ident$rho_values, 1), collapse = ", "),
            summary_stats$ident$B))
cat(sprintf("  Smallest eigenvalue at rho=0: %.4f\n", summary_stats$ident$eigenvalue_at_rho_0))
cat(sprintf("  Smallest eigenvalue at rho=0.9: %.4f\n", summary_stats$ident$eigenvalue_at_rho_09))
cat(sprintf("  Condition number range: [%.2f, %.2f]\n",
            summary_stats$ident$condition_range[1], summary_stats$ident$condition_range[2]))

cat("\n")
cat("==============================================================\n")
cat("  ALL LATEX TABLES GENERATED SUCCESSFULLY\n")
cat("==============================================================\n")
cat("\nTables saved to:", tables_dir, "\n")
