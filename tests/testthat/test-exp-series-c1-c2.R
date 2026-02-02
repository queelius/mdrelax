# Tests for exp_series_c1_c2.R - Relaxed C3 (Power-Weighted Masking)

# -----------------------------------------------------------------------------
# Data generation tests
# -----------------------------------------------------------------------------

test_that("rexp_series_md_c1_c2 generates correct structure", {
    set.seed(42)
    sim <- rexp_series_md_c1_c2(n = 50, theta = c(1, 2), alpha = 1,
                                 base_p = 0.4, tau = 5)
    expect_equal(length(sim$t), 50)
    expect_equal(length(sim$delta), 50)
    expect_equal(dim(sim$C), c(50, 2))
    expect_equal(length(sim$k), 50)
    expect_true(all(sim$delta %in% c(0, 1)))
})

test_that("rexp_series_md_c1_c2 satisfies C1", {
    set.seed(42)
    sim <- rexp_series_md_c1_c2(n = 100, theta = c(1, 2, 3), alpha = 1,
                                 base_p = 0.5, tau = 10)
    uncensored <- sim$delta == 1
    for (i in which(uncensored)) {
        expect_true(sim$C[i, sim$k[i]])
    }
})

test_that("rexp_series_md_c1_c2 with alpha=0 gives uniform masking", {
    set.seed(42)
    n <- 1000
    theta <- c(1, 5)  # Very different rates
    sim <- rexp_series_md_c1_c2(n = n, theta = theta, alpha = 0,
                                 base_p = 0.5, tau = 10)

    # With alpha=0, both components should have same inclusion rate
    # (for non-failed component)
    uncensored <- sim$delta == 1
    # Inclusion rates for component 1 when component 2 failed
    k2 <- uncensored & sim$k == 2
    if (sum(k2) > 50) {
        rate1_given_k2 <- mean(sim$C[k2, 1])
        expect_equal(rate1_given_k2, 0.5, tolerance = 0.1)
    }
})

test_that("rexp_series_md_c1_c2 with large alpha gives informative masking", {
    set.seed(42)
    n <- 1000
    theta <- c(1, 4)  # Component 2 has higher rate
    sim <- rexp_series_md_c1_c2(n = n, theta = theta, alpha = 5,
                                 base_p = 0.5, tau = 10)

    # With large alpha, comp 2 (higher rate) should be more often in C
    uncensored <- sim$delta == 1
    k1 <- uncensored & sim$k == 1  # When comp 1 failed

    if (sum(k1) > 50) {
        rate2_given_k1 <- mean(sim$C[k1, 2])
        # With alpha=5, theta=(1,4): weights proportional to 1^5=1, 4^5=1024
        # w_norm = c(1/1024, 1) * 0.5 ≈ c(0.0005, 0.5)
        # So component 2 should be included ~50% when comp 1 failed
        # Component 1 (low rate) should rarely be included when comp 2 failed
        k2 <- uncensored & sim$k == 2
        if (sum(k2) > 50) {
            rate1_given_k2 <- mean(sim$C[k2, 1])
            expect_true(rate2_given_k1 > rate1_given_k2)
        }
    }
})

test_that("rexp_series_md_c1_c2 censored obs have empty candidate sets", {
    set.seed(42)
    sim <- rexp_series_md_c1_c2(n = 100, theta = c(1, 2), alpha = 1,
                                 base_p = 0.5, tau = 0.2)
    censored <- sim$delta == 0
    if (any(censored)) {
        expect_true(all(rowSums(sim$C[censored, , drop = FALSE]) == 0))
    }
})

# -----------------------------------------------------------------------------
# Log-likelihood tests
# -----------------------------------------------------------------------------

test_that("loglik_exp_series_c1_c2 returns a function (fixed alpha)", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    ll <- loglik_exp_series_c1_c2(t, C, alpha = 0, base_p = 0.5)
    expect_type(ll, "closure")
})

test_that("loglik_exp_series_c1_c2 returns -Inf for non-positive theta", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    ll <- loglik_exp_series_c1_c2(t, C, alpha = 0, base_p = 0.5)
    expect_equal(ll(c(0, 1)), -Inf)
    expect_equal(ll(c(-1, 1)), -Inf)
})

test_that("loglik_exp_series_c1_c2 with alpha=0 matches C1-C2-C3 (up to constant)", {
    set.seed(42)
    sim <- rexp_series_md(n = 50, theta = c(1, 2), p = 0.3, tau = 5)

    ll_c123 <- loglik_exp_series(sim$t, sim$C, sim$delta)
    ll_c12 <- loglik_exp_series_c1_c2(sim$t, sim$C, sim$delta, alpha = 0,
                                       base_p = 0.3)

    theta <- c(1, 2)
    # With alpha=0, the power weights are uniform, so w_norm = base_p for all
    # The masking contribution is constant (doesn't depend on theta)
    # So the functions should differ by a constant
    val_c123 <- ll_c123(theta)
    val_c12 <- ll_c12(theta)

    # Check that they move in the same direction
    theta2 <- c(1.5, 2.5)
    diff1 <- ll_c123(theta2) - val_c123
    diff2 <- ll_c12(theta2) - val_c12
    expect_equal(diff1, diff2, tolerance = 1e-6)
})

test_that("loglik_exp_series_c1_c2 handles censored observations", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    delta <- c(1, 0)
    ll <- loglik_exp_series_c1_c2(t, C, delta, alpha = 0, base_p = 0.5)
    val <- ll(c(1, 1))
    expect_true(is.finite(val))
})

# -----------------------------------------------------------------------------
# Score tests
# -----------------------------------------------------------------------------

test_that("score_exp_series_c1_c2 matches numerical gradient", {
    set.seed(42)
    sim <- rexp_series_md_c1_c2(n = 30, theta = c(1, 2), alpha = 1,
                                 base_p = 0.4, tau = 5)

    ll <- loglik_exp_series_c1_c2(sim$t, sim$C, sim$delta, alpha = 1,
                                   base_p = 0.4)
    sc <- score_exp_series_c1_c2(sim$t, sim$C, sim$delta, alpha = 1,
                                  base_p = 0.4)

    theta <- c(1, 2)
    h <- 1e-6
    num_grad <- numeric(2)
    for (j in 1:2) {
        theta_plus <- theta
        theta_plus[j] <- theta_plus[j] + h
        num_grad[j] <- (ll(theta_plus) - ll(theta)) / h
    }

    expect_equal(sc(theta), num_grad, tolerance = 1e-3)
})

test_that("score with alpha=0 matches C1-C2-C3 score (hazard part)", {
    t <- c(0.5, 1.0, 0.8)
    C <- matrix(c(TRUE, TRUE, FALSE,
                  TRUE, FALSE, TRUE), nrow = 3, ncol = 2)
    delta <- c(1, 1, 1)

    sc_c123 <- score_exp_series(t, C, delta)
    sc_c12 <- score_exp_series_c1_c2(t, C, delta, alpha = 0, base_p = 0.5)

    theta <- c(1, 2)
    # With alpha=0 the masking gradient should be zero
    # So the scores should match
    expect_equal(sc_c12(theta), sc_c123(theta), tolerance = 1e-4)
})

# -----------------------------------------------------------------------------
# FIM tests
# -----------------------------------------------------------------------------

test_that("fim_exp_series_c1_c2 returns correct dimensions", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    fim <- fim_exp_series_c1_c2(t, C, alpha = 1, base_p = 0.5)
    I <- fim(c(1, 2))
    expect_equal(dim(I), c(2, 2))
})

test_that("fim_exp_series_c1_c2 is positive definite at MLE", {
    set.seed(42)
    sim <- rexp_series_md_c1_c2(n = 100, theta = c(1, 2), alpha = 1,
                                 base_p = 0.4, tau = 5)
    fit <- mle_exp_series_c1_c2(sim$t, sim$C, sim$delta, alpha = 1,
                                 base_p = 0.4)
    if (fit$converged) {
        eigs <- eigen(fit$fim)$values
        expect_true(all(eigs > 0))
    }
})

# -----------------------------------------------------------------------------
# MLE tests
# -----------------------------------------------------------------------------

test_that("mle_exp_series_c1_c2 converges with fixed alpha", {
    set.seed(42)
    theta_true <- c(1, 2)
    sim <- rexp_series_md_c1_c2(n = 200, theta = theta_true, alpha = 1,
                                 base_p = 0.4, tau = 5)
    fit <- mle_exp_series_c1_c2(sim$t, sim$C, sim$delta, alpha = 1,
                                 base_p = 0.4)
    expect_true(fit$converged)
    expect_equal(length(fit$theta), 2)
    expect_true(all(fit$theta > 0))
})

test_that("mle_exp_series_c1_c2 with alpha=0 matches C1-C2-C3 MLE", {
    set.seed(42)
    theta_true <- c(1, 2)
    sim <- rexp_series_md(n = 200, theta = theta_true, p = 0.3, tau = 5)

    fit_c123 <- mle_exp_series(sim$t, sim$C, sim$delta)
    fit_c12 <- mle_exp_series_c1_c2(sim$t, sim$C, sim$delta, alpha = 0,
                                     base_p = 0.3)

    if (fit_c123$converged && fit_c12$converged) {
        expect_equal(fit_c12$theta, fit_c123$theta, tolerance = 0.05)
    }
})

test_that("mle_exp_series_c1_c2 returns correct structure", {
    set.seed(42)
    sim <- rexp_series_md_c1_c2(n = 100, theta = c(1, 2), alpha = 1,
                                 base_p = 0.4, tau = 5)
    fit <- mle_exp_series_c1_c2(sim$t, sim$C, sim$delta, alpha = 1,
                                 base_p = 0.4)

    expect_true(is.list(fit))
    expect_true("theta" %in% names(fit))
    expect_true("alpha" %in% names(fit))
    expect_true("se" %in% names(fit))
    expect_true("vcov" %in% names(fit))
    expect_true("loglik" %in% names(fit))
    expect_true("converged" %in% names(fit))
})

test_that("mle_exp_series_c1_c2 with joint alpha estimation", {
    set.seed(42)
    theta_true <- c(1, 2)
    sim <- rexp_series_md_c1_c2(n = 200, theta = theta_true, alpha = 1,
                                 base_p = 0.4, tau = 5)
    fit <- mle_exp_series_c1_c2(sim$t, sim$C, sim$delta, alpha = NULL,
                                 base_p = 0.4)
    expect_true(fit$converged)
    expect_true(fit$alpha >= 0)
})

# -----------------------------------------------------------------------------
# Monte Carlo tests
# -----------------------------------------------------------------------------

test_that("MLE is approximately unbiased with known alpha (Monte Carlo)", {
    skip_on_cran()

    set.seed(42)
    theta_true <- c(1, 2)
    alpha <- 1
    base_p <- 0.4
    n_sim <- 30
    n_obs <- 200

    estimates <- matrix(nrow = n_sim, ncol = 2)
    for (i in seq_len(n_sim)) {
        sim <- rexp_series_md_c1_c2(n = n_obs, theta = theta_true,
                                     alpha = alpha, base_p = base_p, tau = 5)
        fit <- mle_exp_series_c1_c2(sim$t, sim$C, sim$delta,
                                     alpha = alpha, base_p = base_p)
        if (fit$converged) {
            estimates[i, ] <- fit$theta
        }
    }

    mean_est <- colMeans(estimates, na.rm = TRUE)
    expect_equal(mean_est, theta_true, tolerance = 0.3)
})

# -----------------------------------------------------------------------------
# Data frame interface tests
# -----------------------------------------------------------------------------

test_that("loglik_exp_series_c1_c2_df matches direct call", {
    set.seed(42)
    sim <- rexp_series_md_c1_c2(n = 20, theta = c(1, 2), alpha = 1,
                                 base_p = 0.4, tau = 5)
    df <- as_dataframe(sim)

    ll_direct <- loglik_exp_series_c1_c2(sim$t, sim$C, sim$delta,
                                          alpha = 1, base_p = 0.4)
    ll_df <- loglik_exp_series_c1_c2_df(df, alpha = 1, base_p = 0.4)

    theta <- c(1, 2)
    expect_equal(ll_direct(theta), ll_df(theta))
})
