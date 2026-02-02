# Tests for exp_series_relaxed_c1.R - Relaxed C1 (P(K∈C) < 1)

# -----------------------------------------------------------------------------
# Data generation tests
# -----------------------------------------------------------------------------

test_that("rexp_series_md_relaxed_c1 generates correct structure", {
    set.seed(42)
    P <- matrix(c(0.9, 0.3, 0.4, 0.85), nrow = 2)  # Diagonal < 1
    sim <- rexp_series_md_relaxed_c1(n = 50, theta = c(1, 2), P = P, tau = 5)

    expect_equal(length(sim$t), 50)
    expect_equal(length(sim$delta), 50)
    expect_equal(dim(sim$C), c(50, 2))
    expect_equal(length(sim$k), 50)
})

test_that("rexp_series_md_relaxed_c1 can miss the failed component", {
    set.seed(42)
    P <- matrix(c(0.5, 0.3, 0.3, 0.5), nrow = 2)  # Low diagonal
    sim <- rexp_series_md_relaxed_c1(n = 200, theta = c(1, 2), P = P, tau = 10)

    # Some uncensored observations should have k NOT in C
    uncensored <- sim$delta == 1
    k_in_C <- sapply(which(uncensored), function(i) sim$C[i, sim$k[i]])
    # With p_k(k) = 0.5, about 50% should have k not in C
    expect_true(mean(!k_in_C) > 0.2)
})

test_that("rexp_series_md_relaxed_c1 with diag=1 matches C1-C3", {
    set.seed(42)
    P <- matrix(c(1, 0.3, 0.4, 1), nrow = 2)  # Diagonal = 1 (C1 holds)
    sim_r <- rexp_series_md_relaxed_c1(n = 200, theta = c(1, 2), P = P, tau = 5)

    # All uncensored observations should have k in C
    uncensored <- sim_r$delta == 1
    for (i in which(uncensored)) {
        expect_true(sim_r$C[i, sim_r$k[i]])
    }
})

# -----------------------------------------------------------------------------
# Log-likelihood tests
# -----------------------------------------------------------------------------

test_that("loglik_exp_series_relaxed_c1 returns a function", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    P <- matrix(c(0.9, 0.3, 0.3, 0.9), nrow = 2)
    ll <- loglik_exp_series_relaxed_c1(t, C, delta = c(1, 1), P)
    expect_type(ll, "closure")
})

test_that("loglik_exp_series_relaxed_c1 returns -Inf for non-positive theta", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    P <- matrix(c(0.9, 0.3, 0.3, 0.9), nrow = 2)
    ll <- loglik_exp_series_relaxed_c1(t, C, delta = c(1, 1), P)
    expect_equal(ll(c(0, 1)), -Inf)
    expect_equal(ll(c(-1, 1)), -Inf)
})

test_that("loglik_exp_series_relaxed_c1 with diag(P)=1 matches C1-C3", {
    set.seed(42)
    P <- matrix(c(1, 0.3, 0.4, 1), nrow = 2)
    sim <- rexp_series_md_c1_c3(n = 50, theta = c(1, 2), P = P, tau = 5)

    ll_c13 <- loglik_exp_series_c1_c3(sim$t, sim$C, sim$delta, P = P)
    ll_rc1 <- loglik_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P = P)

    theta <- c(1, 2)
    expect_equal(ll_c13(theta), ll_rc1(theta), tolerance = 1e-8)

    theta2 <- c(1.5, 2.5)
    expect_equal(ll_c13(theta2), ll_rc1(theta2), tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
# Score tests
# -----------------------------------------------------------------------------

test_that("score_exp_series_relaxed_c1 matches numerical gradient", {
    set.seed(42)
    P <- matrix(c(0.8, 0.3, 0.4, 0.7), nrow = 2)
    sim <- rexp_series_md_relaxed_c1(n = 30, theta = c(1, 2), P = P, tau = 5)

    ll <- loglik_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)
    sc <- score_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)

    theta <- c(1, 2)
    h <- 1e-6
    num_grad <- numeric(2)
    for (j in 1:2) {
        theta_plus <- theta
        theta_plus[j] <- theta_plus[j] + h
        num_grad[j] <- (ll(theta_plus) - ll(theta)) / h
    }

    expect_equal(sc(theta), num_grad, tolerance = 1e-4)
})

test_that("score_exp_series_relaxed_c1 with diag=1 matches C1-C3 score", {
    set.seed(42)
    P <- matrix(c(1, 0.3, 0.4, 1), nrow = 2)
    sim <- rexp_series_md_c1_c3(n = 30, theta = c(1, 2), P = P, tau = 5)

    sc_c13 <- score_exp_series_c1_c3(sim$t, sim$C, sim$delta, P)
    sc_rc1 <- score_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)

    theta <- c(1, 2)
    expect_equal(sc_c13(theta), sc_rc1(theta), tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
# FIM tests
# -----------------------------------------------------------------------------

test_that("fim_exp_series_relaxed_c1 equals negative Hessian", {
    set.seed(42)
    P <- matrix(c(0.8, 0.3, 0.4, 0.7), nrow = 2)
    sim <- rexp_series_md_relaxed_c1(n = 30, theta = c(1, 2), P = P, tau = 5)

    ll <- loglik_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)
    fim <- fim_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)

    theta <- c(1, 2)
    h <- 1e-5

    # Numerical Hessian
    H <- matrix(0, 2, 2)
    for (j in 1:2) {
        for (k in 1:2) {
            theta_jk <- theta
            theta_jk[j] <- theta_jk[j] + h
            theta_jk[k] <- theta_jk[k] + h

            theta_j <- theta
            theta_j[j] <- theta_j[j] + h

            theta_k <- theta
            theta_k[k] <- theta_k[k] + h

            H[j, k] <- (ll(theta_jk) - ll(theta_j) - ll(theta_k) + ll(theta)) / h^2
        }
    }

    expect_equal(fim(theta), -H, tolerance = 1e-3)
})

test_that("fim_exp_series_relaxed_c1 with diag=1 matches C1-C3 FIM", {
    set.seed(42)
    P <- matrix(c(1, 0.3, 0.4, 1), nrow = 2)
    sim <- rexp_series_md_c1_c3(n = 30, theta = c(1, 2), P = P, tau = 5)

    fim_c13 <- fim_exp_series_c1_c3(sim$t, sim$C, sim$delta, P)
    fim_rc1 <- fim_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)

    theta <- c(1, 2)
    expect_equal(fim_c13(theta), fim_rc1(theta), tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
# MLE tests
# -----------------------------------------------------------------------------

test_that("mle_exp_series_relaxed_c1 converges", {
    set.seed(42)
    P <- matrix(c(0.9, 0.3, 0.4, 0.85), nrow = 2)
    sim <- rexp_series_md_relaxed_c1(n = 200, theta = c(1, 2), P = P, tau = 5)

    fit <- mle_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)
    expect_true(fit$converged)
    expect_equal(length(fit$theta), 2)
})

test_that("mle_exp_series_relaxed_c1 returns correct structure", {
    set.seed(42)
    P <- matrix(c(0.9, 0.3, 0.4, 0.85), nrow = 2)
    sim <- rexp_series_md_relaxed_c1(n = 100, theta = c(1, 2), P = P, tau = 5)

    fit <- mle_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)
    expect_true(is.list(fit))
    expect_true("theta" %in% names(fit))
    expect_true("se" %in% names(fit))
    expect_true("vcov" %in% names(fit))
    expect_true("loglik" %in% names(fit))
    expect_true("converged" %in% names(fit))
    expect_true("fim" %in% names(fit))
})

test_that("mle_exp_series_relaxed_c1 with diag=1 matches C1-C3 MLE", {
    set.seed(42)
    theta_true <- c(1, 2)
    P <- matrix(c(1, 0.3, 0.4, 1), nrow = 2)
    sim <- rexp_series_md_c1_c3(n = 200, theta = theta_true, P = P, tau = 5)

    fit_c13 <- mle_exp_series_c1_c3(sim$t, sim$C, sim$delta, fixed_P = P)
    fit_rc1 <- mle_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)

    if (fit_c13$converged && fit_rc1$converged) {
        expect_equal(fit_rc1$theta, fit_c13$theta, tolerance = 0.01)
    }
})

# -----------------------------------------------------------------------------
# Monte Carlo tests
# -----------------------------------------------------------------------------

test_that("MLE with high P diagonal is approximately unbiased (Monte Carlo)", {
    skip_on_cran()

    set.seed(42)
    theta_true <- c(1, 2)
    P <- matrix(c(0.95, 0.3, 0.4, 0.9), nrow = 2)
    n_sim <- 30
    n_obs <- 200

    estimates <- matrix(nrow = n_sim, ncol = 2)
    for (i in seq_len(n_sim)) {
        sim <- rexp_series_md_relaxed_c1(n = n_obs, theta = theta_true,
                                          P = P, tau = 5)
        fit <- mle_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)
        if (fit$converged) {
            estimates[i, ] <- fit$theta
        }
    }

    mean_est <- colMeans(estimates, na.rm = TRUE)
    # Relaxed C1 with high diagonal should still give reasonable estimates
    expect_equal(mean_est, theta_true, tolerance = 0.4)
})

# -----------------------------------------------------------------------------
# Edge cases
# -----------------------------------------------------------------------------

test_that("relaxed C1 with very low diagonal still works", {
    set.seed(42)
    P <- matrix(c(0.6, 0.3, 0.3, 0.6), nrow = 2)
    sim <- rexp_series_md_relaxed_c1(n = 200, theta = c(1, 2), P = P, tau = 5)

    fit <- mle_exp_series_relaxed_c1(sim$t, sim$C, sim$delta, P)
    expect_true(fit$converged)
    # With low diagonal, estimates may be biased but should be finite
    expect_true(all(is.finite(fit$theta)))
})

test_that("relaxed C1 handles all observations censored", {
    t <- c(2, 3)
    C <- matrix(FALSE, nrow = 2, ncol = 2)
    delta <- c(0, 0)
    P <- matrix(c(0.9, 0.3, 0.3, 0.9), nrow = 2)

    ll <- loglik_exp_series_relaxed_c1(t, C, delta, P)
    # With all censored, ll is just -sum(t) * sum(theta)
    expect_equal(ll(c(1, 1)), -5 * 2)
})
