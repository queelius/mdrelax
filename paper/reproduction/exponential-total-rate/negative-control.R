#!/usr/bin/env Rscript
# Deliberately wrong likelihood: omit censored exposure. The verifier must fail.
suppressPackageStartupMessages(library(mdrelax))
sandbox <- new.env(parent = globalenv())
sys.source("R/exp_series_c1_c2_c3.R", envir = sandbox)
sandbox$correct_loglik <- sandbox$loglik_exp_series
evalq(loglik_exp_series <- function(t, C, delta = NULL) {
    if (is.null(delta)) delta <- rep(1, length(t))
    t[delta == 0] <- 0
    correct_loglik(t, C, delta)
}, envir = sandbox)
result <- tryCatch({
    sys.source("paper/reproduction/exponential-total-rate/verify.R", envir = sandbox)
    NULL
}, error = identity)
stopifnot(inherits(result, "simpleError"),
          identical(conditionMessage(result),
                    "max(abs(actual - expected)) <= tolerance is not TRUE"))
cat("Rejected the omitted-censored-exposure mutation at the independent likelihood check.\n")
