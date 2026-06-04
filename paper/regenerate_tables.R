#!/usr/bin/env Rscript
# =============================================================================
# Regenerate Tables 1, 2, and 3 of the paper from the canonical sweep RDS
# files. All three tables come from a single code path so the numbers
# in the manuscript stay consistent with paper/data/*.rds.
# =============================================================================
#
# Inputs:
#   paper/data/sensitivity_sweep_c1.rds      (C1 sweep at n=500, B=200)
#   paper/data/sensitivity_sweep.rds         (C2 sweep at n=500, B=200)
#   paper/data/sensitivity_sweep_c3.rds      (C3 sweep at n=500, B=200)
#   paper/data/sensitivity_sweep_c1_n2000.rds (C1 ablation at n=2000)
#   paper/data/sensitivity_sweep_n2000.rds    (C2 ablation at n=2000)
#   paper/data/sensitivity_sweep_c3_n2000.rds (C3 ablation at n=2000)
#
# Outputs:
#   Prints the three tables as plain text and as LaTeX tabular bodies
#   to stdout. Optionally writes the LaTeX bodies to a file when --write
#   is supplied (so they can be diffed against the manuscript).
#
# Usage:
#   Rscript paper/regenerate_tables.R [--write]
#
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
write_out <- "--write" %in% args

script_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error = function(e) normalizePath("paper")
)
data_dir <- file.path(script_dir, "data")

# -----------------------------------------------------------------------------
# Load all six sweep results
# -----------------------------------------------------------------------------

df_c1    <- readRDS(file.path(data_dir, "sensitivity_sweep_c1.rds"))
df_c2    <- readRDS(file.path(data_dir, "sensitivity_sweep.rds"))
df_c3    <- readRDS(file.path(data_dir, "sensitivity_sweep_c3.rds"))
df_c1_2k <- readRDS(file.path(data_dir, "sensitivity_sweep_c1_n2000.rds"))
df_c2_2k <- readRDS(file.path(data_dir, "sensitivity_sweep_n2000.rds"))
df_c3_2k <- readRDS(file.path(data_dir, "sensitivity_sweep_c3_n2000.rds"))

# Each frame has the same schema: component, bias, rmse, coverage,
# median_bias, mad, mean_est, true_value, n_converged, n_reps, plus the
# severity column named "alpha" (C1, C3) or "severity" (C2).

# -----------------------------------------------------------------------------
# Table 1: Worst component-level relative bias per violation
# -----------------------------------------------------------------------------

#' For each component (excluding sys_hazard), find the severity at which the
#' absolute relative bias |bias|/|true_value| is largest. Return the worst
#' three components, ranked by that maximum.
worst_component_fragility <- function(df, sev_col, top_n = 3L) {
    df$rel_bias <- abs(df$bias) / abs(df$true_value)
    comps <- setdiff(unique(df$component), "sys_hazard")

    rows <- lapply(comps, function(c) {
        sub <- df[df$component == c, ]
        idx <- which.max(sub$rel_bias)
        data.frame(
            component        = c,
            pct_bias         = 100 * sub$rel_bias[idx],
            coverage         = sub$coverage[idx],
            severity         = sub[[sev_col]][idx],
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out <- out[order(-out$pct_bias), ]
    head(out, top_n)
}

#' Pretty-print a component name (k1 -> $k_1$, lambda3 -> $\lambda_3$) for
#' LaTeX output.
tex_comp <- function(name) {
    if (grepl("^k([0-9]+)$", name)) {
        idx <- sub("^k", "", name)
        sprintf("$k_%s$", idx)
    } else if (grepl("^lambda([0-9]+)$", name)) {
        idx <- sub("^lambda", "", name)
        sprintf("$\\lambda_%s$", idx)
    } else {
        name
    }
}

cov_str <- function(c) {
    sprintf("%.2f", c)
}

sev_str <- function(label, sev) {
    if (label == "C1") {
        sprintf("$\\alpha_{C_1} = %.2f$", sev)
    } else if (label == "C2") {
        sprintf("$s = %.1f$", sev)
    } else {
        sprintf("$\\alpha_{C_3} = %.2f$", sev)
    }
}

emit_table1 <- function(out_path = NULL) {
    fragility_block <- function(label, df_sub) {
        lines <- character()
        first <- TRUE
        for (i in seq_len(nrow(df_sub))) {
            row <- df_sub[i, ]
            cell_violation <- if (first) label else ""
            lines <- c(lines, sprintf("%s & %s & %.0f\\%% & %s & %s \\\\",
                cell_violation,
                tex_comp(row$component),
                row$pct_bias,
                cov_str(row$coverage),
                sev_str(label, row$severity)
            ))
            first <- FALSE
        }
        lines
    }

    c1_block <- fragility_block("C1", worst_component_fragility(df_c1, "alpha"))
    c2_block <- fragility_block("C2", worst_component_fragility(df_c2, "severity"))
    c3_block <- fragility_block("C3", worst_component_fragility(df_c3, "alpha"))

    body <- c(
        "\\toprule",
        "Violation & Parameter & Max $|\\text{bias}|/|\\text{truth}|$ & Coverage at worst & Severity \\\\",
        "\\midrule",
        c1_block,
        "\\midrule",
        c2_block,
        "\\midrule",
        c3_block,
        "\\bottomrule"
    )

    cat("\n%% ---- Table 1: Component-level fragility (regenerated) ----\n")
    cat(paste(body, collapse = "\n"), "\n", sep = "")

    if (!is.null(out_path)) {
        writeLines(body, out_path)
    }
}

# -----------------------------------------------------------------------------
# Tables 2 and 3: First-order behavior test (quadratic linear term)
# -----------------------------------------------------------------------------

#' Fit bias ~ alpha + I(alpha^2) on the sys_hazard subset and return the
#' linear coefficient and its p-value. This is the test that produced the
#' verdicts in Tables 2 and 3 of the manuscript.
quadratic_linear_term <- function(df, sev_col) {
    sys <- df[df$component == "sys_hazard", ]
    sys$sev <- sys[[sev_col]]
    fit <- lm(bias ~ sev + I(sev^2), data = sys)
    cs <- summary(fit)$coefficients
    list(linear = cs[2, 1], pval = cs[2, 4])
}

#' Format coefficient with two decimals normally, but switch to three
#' decimals when the magnitude is below 0.01 so very small slopes do not
#' collapse to "+0.00" and lose information.
fmt_coef <- function(x) {
    if (abs(x) < 0.01) sprintf("%+.3f", x) else sprintf("%+.2f", x)
}

fmt_p <- function(p) {
    if (p < 0.001) "$< 0.001$" else sprintf("%.3f", p)
}

verdict_label <- function(linear, pval) {
    if (pval < 0.05) "broken" else "preserved"
}

theory_match <- function(label, verdict) {
    if (label == "C1") {
        "matches \\Cref{cor:hierarchy} item~3"
    } else if (label == "C2") {
        "matches \\Cref{prop:hazard-robust}"
    } else {
        "matches \\Cref{cor:hierarchy} item~3"
    }
}

emit_table2 <- function(out_path = NULL) {
    t_c1 <- quadratic_linear_term(df_c1, "alpha")
    t_c2 <- quadratic_linear_term(df_c2, "severity")
    t_c3 <- quadratic_linear_term(df_c3, "alpha")

    fmt_row <- function(label, t) {
        v <- verdict_label(t$linear, t$pval)
        sprintf("%s & $%s$  & %s & %s          & %s \\\\",
                label, fmt_coef(t$linear), fmt_p(t$pval), v, theory_match(label, v))
    }

    body <- c(
        "\\toprule",
        "Violation & Quadratic linear term & $p$-value & First-order verdict & Theory match \\\\",
        "\\midrule",
        fmt_row("C1", t_c1),
        fmt_row("C2", t_c2),
        fmt_row("C3", t_c3),
        "\\bottomrule"
    )

    cat("\n%% ---- Table 2: Slope tests at n=500 (regenerated) ----\n")
    cat(paste(body, collapse = "\n"), "\n", sep = "")

    if (!is.null(out_path)) {
        writeLines(body, out_path)
    }
}

emit_table3 <- function(out_path = NULL) {
    t1 <- quadratic_linear_term(df_c1,    "alpha")
    t2 <- quadratic_linear_term(df_c2,    "severity")
    t3 <- quadratic_linear_term(df_c3,    "alpha")
    u1 <- quadratic_linear_term(df_c1_2k, "alpha")
    u2 <- quadratic_linear_term(df_c2_2k, "severity")
    u3 <- quadratic_linear_term(df_c3_2k, "alpha")

    fmt_row <- function(label, lo, hi) {
        sprintf("%-13s & $%s$      & %s  & $%s$      & %s  \\\\",
                label, fmt_coef(lo$linear), fmt_p(lo$pval),
                fmt_coef(hi$linear), fmt_p(hi$pval))
    }

    body <- c(
        "\\toprule",
        "              & \\multicolumn{2}{c}{$n = 500$, $B = 200$} & \\multicolumn{2}{c}{$n = 2000$, $B = 100$} \\\\",
        "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
        "Violation     & Linear coef. & $p$        & Linear coef. & $p$        \\\\",
        "\\midrule",
        fmt_row("C1", t1, u1),
        fmt_row("C2", t2, u2),
        fmt_row("C3", t3, u3),
        "\\bottomrule"
    )

    cat("\n%% ---- Table 3: Slope tests at two sample sizes (regenerated) ----\n")
    cat(paste(body, collapse = "\n"), "\n", sep = "")

    if (!is.null(out_path)) {
        writeLines(body, out_path)
    }
}

# -----------------------------------------------------------------------------
# Diagnostic print: every row that appears in Table 1, including the
# components that fall below the top-3 cutoff. Useful for re-checking by
# hand.
# -----------------------------------------------------------------------------

cat("=== Component fragility (top 5 per violation, full numeric) ===\n\n")
cat("--- C1 ---\n")
print(worst_component_fragility(df_c1, "alpha", top_n = 5L), row.names = FALSE)
cat("\n--- C2 ---\n")
print(worst_component_fragility(df_c2, "severity", top_n = 5L), row.names = FALSE)
cat("\n--- C3 ---\n")
print(worst_component_fragility(df_c3, "alpha", top_n = 5L), row.names = FALSE)

cat("\n=== Quadratic linear-term test (sys_hazard ~ severity) ===\n")
for (lbl_df in list(
    list("C1 n=500",  df_c1,    "alpha"),
    list("C2 n=500",  df_c2,    "severity"),
    list("C3 n=500",  df_c3,    "alpha"),
    list("C1 n=2000", df_c1_2k, "alpha"),
    list("C2 n=2000", df_c2_2k, "severity"),
    list("C3 n=2000", df_c3_2k, "alpha")
)) {
    t <- quadratic_linear_term(lbl_df[[2]], lbl_df[[3]])
    cat(sprintf("  %-10s  linear = %+.5f  p = %.4f\n",
                lbl_df[[1]], t$linear, t$pval))
}

# -----------------------------------------------------------------------------
# Emit LaTeX bodies
# -----------------------------------------------------------------------------

out_dir <- file.path(script_dir, "tables")
if (write_out) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

emit_table1(if (write_out) file.path(out_dir, "table1_body.tex") else NULL)
emit_table2(if (write_out) file.path(out_dir, "table2_body.tex") else NULL)
emit_table3(if (write_out) file.path(out_dir, "table3_body.tex") else NULL)

cat("\n=== Done. Re-run with --write to dump LaTeX bodies into paper/tables/. ===\n")
