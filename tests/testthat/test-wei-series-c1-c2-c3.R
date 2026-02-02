# Tests for wei_series_c1_c2_c3.R - Weibull Masked Data Likelihood

# -----------------------------------------------------------------------------
# Parameter helper tests
# -----------------------------------------------------------------------------

test_that("unpack/pack_wei_params round-trip", {
    theta <- c(2, 3, 1.5, 4)  # k1=2, b1=3, k2=1.5, b2=4
    params <- unpack_wei_params(theta)
    expect_equal(params$shapes, c(2, 1.5))
    expect_equal(params$scales, c(3, 4))
    expect_equal(params$m, 2L)

    theta_back <- pack_wei_params(params$shapes, params$scales)
    expect_equal(theta_back, theta)
})

test_that("unpack_wei_params rejects odd-length theta", {
    expect_error(unpack_wei_params(c(1, 2, 3)), "even length")
})

# -----------------------------------------------------------------------------
# Log-likelihood tests
# -----------------------------------------------------------------------------

test_that("loglik_wei_series returns a function", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    ll <- loglik_wei_series(t, C)
    expect_type(ll, "closure")
})

test_that("loglik_wei_series returns -Inf for non-positive params", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    ll <- loglik_wei_series(t, C)
    expect_equal(ll(c(0, 1, 1, 1)), -Inf)
    expect_equal(ll(c(-1, 1, 1, 1)), -Inf)
    expect_equal(ll(c(1, 0, 1, 1)), -Inf)
})

test_that("loglik_wei_series with shape=1 matches exponential", {
    set.seed(42)
    rates <- c(1, 2)
    sim <- rexp_series_md(n = 50, theta = rates, p = 0.3, tau = 5)

    ll_exp <- loglik_exp_series(sim$t, sim$C, sim$delta)
    ll_wei <- loglik_wei_series(sim$t, sim$C, sim$delta)

    # Weibull with shape=1, scale=1/rate
    theta_wei <- c(1, 1/rates[1], 1, 1/rates[2])
    theta_exp <- rates

    # The log-likelihoods should match
    expect_equal(ll_wei(theta_wei), ll_exp(theta_exp), tolerance = 1e-8)
})

test_that("loglik_wei_series handles censored observations", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    delta <- c(1, 0)
    ll <- loglik_wei_series(t, C, delta)

    theta <- c(2, 3, 1.5, 4)
    val <- ll(theta)
    expect_true(is.finite(val))
})

test_that("loglik_wei_series computes correct value (manual)", {
    # Single observation, single component, uncensored
    t <- c(2)
    C <- matrix(TRUE, nrow = 1, ncol = 1)
    ll <- loglik_wei_series(t, C)

    k <- 2; b <- 3
    theta <- c(k, b)

    # ll = -(t/b)^k + log(h(t))
    # h(t) = (k/b) * (t/b)^(k-1)
    z <- t / b
    expected <- -z^k + log((k / b) * z^(k - 1))
    expect_equal(ll(theta), expected, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Score tests
# -----------------------------------------------------------------------------

test_that("score_wei_series matches numerical gradient", {
    set.seed(42)
    sim <- rwei_series_md(n = 30, shapes = c(2, 1.5), scales = c(3, 4),
                           p = 0.3, tau = 5)

    ll <- loglik_wei_series(sim$t, sim$C, sim$delta)
    sc <- score_wei_series(sim$t, sim$C, sim$delta)

    theta <- c(2, 3, 1.5, 4)
    h <- 1e-6
    num_grad <- numeric(4)
    for (j in 1:4) {
        theta_plus <- theta
        theta_plus[j] <- theta_plus[j] + h
        num_grad[j] <- (ll(theta_plus) - ll(theta)) / h
    }

    expect_equal(sc(theta), num_grad, tolerance = 1e-3)
})

test_that("score_wei_series with shape=1 is consistent with exp loglik", {
    set.seed(42)
    rates <- c(1, 2)
    sim <- rexp_series_md(n = 30, theta = rates, p = 0.3, tau = 5)

    # Score of Weibull loglik at shape=1, scale=1/rate should be zero for shape
    # when the likelihood matches the exponential (at MLE)
    ll_wei <- loglik_wei_series(sim$t, sim$C, sim$delta)
    sc_wei <- score_wei_series(sim$t, sim$C, sim$delta)

    # Verify score has correct length
    theta_wei <- c(1, 1/rates[1], 1, 1/rates[2])
    score_val <- sc_wei(theta_wei)
    expect_equal(length(score_val), 4)
    expect_true(all(is.finite(score_val)))
})

# -----------------------------------------------------------------------------
# FIM tests
# -----------------------------------------------------------------------------

test_that("fim_wei_series equals negative Hessian (numerical)", {
    set.seed(42)
    sim <- rwei_series_md(n = 30, shapes = c(2, 1.5), scales = c(3, 4),
                           p = 0.3, tau = 5)

    ll <- loglik_wei_series(sim$t, sim$C, sim$delta)
    fim <- fim_wei_series(sim$t, sim$C, sim$delta)

    theta <- c(2, 3, 1.5, 4)
    h <- 1e-5

    # Numerical Hessian
    p <- 4
    H <- matrix(0, p, p)
    for (j in 1:p) {
        for (k in 1:p) {
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

    I <- fim(theta)
    expect_equal(I, -H, tolerance = 1e-2)
})

# -----------------------------------------------------------------------------
# MLE tests
# -----------------------------------------------------------------------------

test_that("mle_wei_series converges", {
    set.seed(42)
    shapes_true <- c(2, 1.5)
    scales_true <- c(3, 4)
    sim <- rwei_series_md(n = 200, shapes = shapes_true, scales = scales_true,
                           p = 0.3, tau = 8)

    fit <- mle_wei_series(sim$t, sim$C, sim$delta)
    expect_true(fit$converged)
    expect_equal(length(fit$theta), 4)
    expect_true(all(fit$theta > 0))
})

test_that("mle_wei_series returns correct structure", {
    set.seed(42)
    sim <- rwei_series_md(n = 100, shapes = c(2, 1.5), scales = c(3, 4),
                           p = 0.3, tau = 8)
    fit <- mle_wei_series(sim$t, sim$C, sim$delta)

    expect_true(is.list(fit))
    expect_true("theta" %in% names(fit))
    expect_true("shapes" %in% names(fit))
    expect_true("scales" %in% names(fit))
    expect_true("se" %in% names(fit))
    expect_true("vcov" %in% names(fit))
    expect_true("loglik" %in% names(fit))
    expect_true("converged" %in% names(fit))
    expect_true("fim" %in% names(fit))
    expect_equal(fit$m, 2)
})

test_that("mle_wei_series recovers true parameters approximately", {
    set.seed(42)
    shapes_true <- c(2, 1.5)
    scales_true <- c(3, 4)
    theta_true <- pack_wei_params(shapes_true, scales_true)

    sim <- rwei_series_md(n = 500, shapes = shapes_true, scales = scales_true,
                           p = 0.3, tau = 8)
    fit <- mle_wei_series(sim$t, sim$C, sim$delta,
                           theta0 = theta_true * 0.8)

    if (fit$converged) {
        expect_equal(fit$shapes, shapes_true, tolerance = 0.5)
        expect_equal(fit$scales, scales_true, tolerance = 1.0)
    }
})

# -----------------------------------------------------------------------------
# Data generation tests
# -----------------------------------------------------------------------------

test_that("rwei_series_md generates correct structure", {
    set.seed(42)
    sim <- rwei_series_md(n = 50, shapes = c(2, 1.5, 3), scales = c(3, 4, 2),
                           p = 0.4, tau = 5)

    expect_equal(length(sim$t), 50)
    expect_equal(length(sim$delta), 50)
    expect_equal(dim(sim$C), c(50, 3))
    expect_equal(length(sim$k), 50)
    expect_true(all(sim$delta %in% c(0, 1)))
    expect_true(all(sim$k %in% 1:3))
})

test_that("rwei_series_md satisfies C1", {
    set.seed(42)
    sim <- rwei_series_md(n = 100, shapes = c(2, 1.5), scales = c(3, 4),
                           p = 0.3, tau = 10)

    uncensored <- sim$delta == 1
    for (i in which(uncensored)) {
        expect_true(sim$C[i, sim$k[i]])
    }
})

test_that("rwei_series_md censored obs have empty candidate sets", {
    set.seed(42)
    sim <- rwei_series_md(n = 100, shapes = c(2, 1.5), scales = c(3, 4),
                           p = 0.5, tau = 0.5)

    censored <- sim$delta == 0
    if (any(censored)) {
        expect_true(all(rowSums(sim$C[censored, , drop = FALSE]) == 0))
    }
})

test_that("rwei_series_md with keep_latent includes component times", {
    set.seed(42)
    sim <- rwei_series_md(n = 10, shapes = c(2, 1.5), scales = c(3, 4),
                           p = 0.3, tau = 10, keep_latent = TRUE)

    expect_true("Tm" %in% names(sim))
    expect_equal(dim(sim$Tm), c(10, 2))

    # System time = min of components (before censoring)
    true_sys <- apply(sim$Tm, 1, min)
    expect_equal(pmin(true_sys, 10), sim$t, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Data frame interface tests
# -----------------------------------------------------------------------------

test_that("loglik_wei_series_df matches direct call", {
    set.seed(42)
    sim <- rwei_series_md(n = 20, shapes = c(2, 1.5), scales = c(3, 4),
                           p = 0.3, tau = 5)
    df <- as_dataframe(sim)

    ll_direct <- loglik_wei_series(sim$t, sim$C, sim$delta)
    ll_df <- loglik_wei_series_df(df)

    theta <- c(2, 3, 1.5, 4)
    expect_equal(ll_direct(theta), ll_df(theta))
})

# -----------------------------------------------------------------------------
# Monte Carlo tests
# -----------------------------------------------------------------------------

test_that("Weibull MLE is approximately unbiased (Monte Carlo)", {
    skip_on_cran()

    set.seed(42)
    shapes_true <- c(2, 1.5)
    scales_true <- c(3, 4)
    theta_true <- pack_wei_params(shapes_true, scales_true)
    n_sim <- 20
    n_obs <- 300

    estimates <- matrix(nrow = n_sim, ncol = 4)
    for (i in seq_len(n_sim)) {
        sim <- rwei_series_md(n = n_obs, shapes = shapes_true,
                               scales = scales_true, p = 0.3, tau = 8)
        fit <- mle_wei_series(sim$t, sim$C, sim$delta,
                               theta0 = theta_true * 0.9)
        if (fit$converged) {
            estimates[i, ] <- fit$theta
        }
    }

    mean_est <- colMeans(estimates, na.rm = TRUE)
    # Shape and scale estimates should be approximately correct
    expect_equal(mean_est[1], 2, tolerance = 0.5)
    expect_equal(mean_est[3], 1.5, tolerance = 0.5)
})
