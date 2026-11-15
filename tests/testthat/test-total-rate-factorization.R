# Independent controls for the existing paper's exponential factorization.

test_that("candidate sets affect composition but not the total-rate likelihood", {
    t <- c(0.2, 0.4, 0.7, 1.1, 1.6)
    delta <- c(1, 1, 1, 0, 0)
    C <- rbind(c(TRUE, FALSE, FALSE), c(TRUE, TRUE, FALSE),
               c(FALSE, TRUE, TRUE), c(FALSE, FALSE, FALSE),
               c(FALSE, FALSE, FALSE))
    phi <- c(0.2, 0.3, 0.5)
    d <- sum(delta)
    exposure <- sum(t)
    composition <- sum(log(drop(C[delta == 1, , drop = FALSE] %*% phi)))
    for (S in c(0.1, d / exposure, 3)) {
        expect_equal(loglik_exp_series(t, C, delta)(S * phi),
                     d * log(S) - S * exposure + composition,
                     tolerance = 1e-12)
        expect_equal(sum(phi * score_exp_series(t, C, delta)(S * phi)),
                     d / S - exposure, tolerance = 1e-12)
        expect_equal(drop(t(phi) %*% fim_exp_series(t, C, delta)(S * phi) %*% phi),
                     d / S^2, tolerance = 1e-12)
    }
})

test_that("relabelling singleton causes changes rates while preserving the total", {
    t <- seq_len(12) / 10
    delta <- c(rep(1, 10), 0, 0)
    causes <- c(rep(1, 5), rep(2, 3), rep(3, 2))
    C <- matrix(FALSE, 12, 3)
    C[cbind(seq_along(causes), causes)] <- TRUE
    fit <- mle_exp_series(t, C, delta)
    relabelled <- mle_exp_series(t, C[, 3:1], delta)
    expect_true(fit$converged)
    expect_true(relabelled$converged)
    expect_equal(fit$theta, c(5, 3, 2) / sum(t), tolerance = 2e-5)
    expect_equal(relabelled$theta, rev(fit$theta), tolerance = 2e-5)
    expect_equal(sum(fit$theta), sum(delta) / sum(t), tolerance = 2e-5)
    expect_false(isTRUE(all.equal(sum(fit$theta), sum(delta) / sum(t[delta == 1]))))
})

test_that("full masking identifies the total without identifying composition", {
    t <- c(0.1, 0.3, 0.7, 1.0)
    C <- matrix(TRUE, length(t), 3)
    ll <- loglik_exp_series(t, C)
    S <- length(t) / sum(t)
    expect_equal(ll(S * c(0.1, 0.2, 0.7)), ll(S * c(0.6, 0.3, 0.1)))
    expect_equal(qr(fim_exp_series(t, C)(rep(S / 3, 3)))$rank, 1L)
    fit <- mle_exp_series(t, C)
    expect_true(fit$converged)
    expect_equal(sum(fit$theta), S, tolerance = 2e-5)
})

test_that("zero failures and binding rate bounds are not interior optima", {
    t <- c(0.2, 0.5, 1)
    C <- matrix(FALSE, length(t), 2)
    ll <- loglik_exp_series(t, C, rep(0, length(t)))
    expect_gt(ll(c(0.01, 0.01)), ll(c(0.1, 0.1)))
    C[,] <- TRUE
    fit <- mle_exp_series(t, C, theta0 = c(0.05, 0.05), upper = 0.1)
    expect_true(fit$converged)
    expect_equal(fit$theta, c(0.1, 0.1))
    expect_lt(sum(fit$theta), length(t) / sum(t))
})

test_that("the same fitted identity does not remove selective-retention bias", {
    rate <- 3
    t <- qexp((seq_len(1000) - 0.5) / 1000, rate = rate)
    retained <- t[t > qexp(0.5, rate = rate)]
    full_rate <- length(t) / sum(t)
    selected_rate <- length(retained) / sum(retained)
    expect_lt(abs(full_rate - rate), 0.01)
    expect_gt(abs(selected_rate - rate), 1)
    expect_lt(abs(selected_rate - rate / (1 + log(2))), 0.01)
    C <- matrix(TRUE, length(retained), 2)
    fit <- mle_exp_series(retained, C)
    expect_true(fit$converged)
    expect_equal(sum(fit$theta), selected_rate, tolerance = 2e-5)
})
