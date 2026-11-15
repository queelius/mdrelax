#!/usr/bin/env Rscript
# Run with an installed mdrelax package; no devtools or workspace loading.
# GPL-3.0-or-later. Alexander Richard Towell; substantial AI assistance.
suppressPackageStartupMessages(library(mdrelax))

near <- function(actual, expected, tolerance = 2e-5) {
    stopifnot(all(is.finite(actual)), all(is.finite(expected)),
              max(abs(actual - expected)) <= tolerance)
}

t <- seq_len(12) / 10
delta <- c(rep(1, 10), 0, 0)
causes <- c(rep(1, 5), rep(2, 3), rep(3, 2))
C <- matrix(FALSE, length(t), 3)
C[cbind(seq_along(causes), causes)] <- TRUE
full <- matrix(TRUE, length(t), 3)
full[delta == 0, ] <- FALSE
cases <- list(singletons = C, relabelled = C[, 3:1], fully_masked = full)
expected_total <- sum(delta) / sum(t)
phi <- c(0.2, 0.3, 0.5)

rows <- lapply(names(cases), function(label) {
    candidates <- cases[[label]]
    ll <- loglik_exp_series(t, candidates, delta)
    composition <- sum(log(drop(candidates[delta == 1, , drop = FALSE] %*% phi)))
    for (total in c(0.1, expected_total, 4)) {
        near(ll(total * phi), sum(delta) * log(total) - total * sum(t) + composition, 1e-10)
        near(sum(phi * score_exp_series(t, candidates, delta)(total * phi)),
             sum(delta) / total - sum(t), 1e-10)
    }
    fit <- mle_exp_series(t, candidates, delta)
    stopifnot(fit$converged)
    near(sum(fit$theta), expected_total)
    data.frame(observation = label, rate_1 = fit$theta[1],
               rate_2 = fit$theta[2], rate_3 = fit$theta[3],
               total = sum(fit$theta), failures = sum(delta), exposure = sum(t))
})
summary <- do.call(rbind, rows)
near(as.numeric(summary[1, c("rate_1", "rate_2", "rate_3")]), c(5, 3, 2) / sum(t))
near(as.numeric(summary[2, c("rate_1", "rate_2", "rate_3")]), c(2, 3, 5) / sum(t))
stopifnot(qr(fim_exp_series(t, full, delta)(expected_total * phi))$rank == 1)
scalar_fit <- optimize(function(total) sum(delta) * log(total) - total * sum(t),
                       c(0.01, 10), maximum = TRUE, tol = 1e-10)
near(scalar_fit$maximum, expected_total)

# Reject a common false inference: ignoring censored exposure changes the rate.
stopifnot(abs(sum(delta) / sum(t[delta == 1]) - expected_total) > 0.5)

# Deterministic quantile grid: retaining only late failures changes the target.
true_rate <- 3
times <- qexp((seq_len(1000) - 0.5) / 1000, rate = true_rate)
retained <- times[times > qexp(0.5, rate = true_rate)]
selected <- length(retained) / sum(retained)
near(length(times) / sum(times), true_rate, 0.01)
near(selected, true_rate / (1 + log(2)), 0.01)
stopifnot(abs(selected - true_rate) > 1)
selected_fit <- mle_exp_series(retained, matrix(TRUE, length(retained), 2))
stopifnot(selected_fit$converged)
near(sum(selected_fit$theta), selected)

# No failures gives a boundary supremum; upper bounds can bind away from d/E.
no_failures <- loglik_exp_series(t, matrix(FALSE, length(t), 3), rep(0, length(t)))
stopifnot(no_failures(rep(0.01, 3)) > no_failures(rep(0.1, 3)))
bounded <- mle_exp_series(t, full, delta, theta0 = rep(0.05, 3), upper = 0.1)
stopifnot(bounded$converged, sum(bounded$theta) < expected_total)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L)
if (length(args) == 1L) {
    dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
    write.csv(summary, file.path(args[1], "fits.csv"), row.names = FALSE)
    write.csv(data.frame(population_rate = true_rate,
                         full_sample_fit = length(times) / sum(times),
                         selected_sample_fit = selected,
                         selected_population_limit = true_rate / (1 + log(2))),
              file.path(args[1], "selection.csv"), row.names = FALSE)
}

write.table(summary, stdout(), sep = ",", row.names = FALSE, quote = FALSE)
cat(sprintf("\nSelective-retention control: full %.6f; selected %.6f; population rate %.1f\n",
            length(times) / sum(times), selected, true_rate))
cat("All factorization, fitting, identification and scope controls passed.\n")
print(sessionInfo())
