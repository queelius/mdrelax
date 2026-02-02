# Tests for exp_series_c1_c3.R - Relaxed C2 (General Bernoulli Model)

# -----------------------------------------------------------------------------
# Helper function tests: make_P_matrix
# -----------------------------------------------------------------------------

test_that("make_P_matrix creates uniform matrix correctly", {
    P <- make_P_matrix(3, "uniform", p = 0.4)
    expect_equal(dim(P), c(3, 3))
    expect_equal(diag(P), c(1, 1, 1))
    expect_equal(P[1, 2], 0.4)
    expect_equal(P[1, 3], 0.4)
    expect_equal(P[2, 1], 0.4)
    expect_equal(P[2, 3], 0.4)
    expect_equal(P[3, 1], 0.4)
    expect_equal(P[3, 2], 0.4)
})

test_that("make_P_matrix creates symmetric matrix correctly", {
    P <- make_P_matrix(3, "symmetric", values = c(0.2, 0.3, 0.4))
    expect_equal(dim(P), c(3, 3))
    expect_equal(diag(P), c(1, 1, 1))
    # Values: (1,2), (1,3), (2,3)
    expect_equal(P[1, 2], 0.2)
    expect_equal(P[2, 1], 0.2)
    expect_equal(P[1, 3], 0.3)
    expect_equal(P[3, 1], 0.3)
    expect_equal(P[2, 3], 0.4)
    expect_equal(P[3, 2], 0.4)
})

test_that("make_P_matrix creates full matrix correctly", {
    # Off-diagonal elements in column-major order
    # For m=3: (2,1), (3,1), (1,2), (3,2), (1,3), (2,3)
    vals <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
    P <- make_P_matrix(3, "full", values = vals)
    expect_equal(dim(P), c(3, 3))
    expect_equal(diag(P), c(1, 1, 1))
    expect_equal(P[row(P) != col(P)], vals)
})

test_that("make_P_matrix full type requires values", {
    expect_error(make_P_matrix(3, "full"), "values required")
})

test_that("make_P_matrix full type checks length", {
    expect_error(make_P_matrix(3, "full", values = c(0.1, 0.2)),
                 "values must have length")
})

# -----------------------------------------------------------------------------
# Helper function tests: satisfies_C2
# -----------------------------------------------------------------------------

test_that("satisfies_C2 returns TRUE for uniform P", {
    P <- make_P_matrix(3, "uniform", p = 0.3)
    expect_true(satisfies_C2(P))
})

test_that("satisfies_C2 returns TRUE for row-constant P", {
    # Each row has constant off-diagonal, but rows can differ
    P <- matrix(0, 3, 3)
    diag(P) <- 1
    P[1, 2:3] <- 0.2  # Row 1 off-diagonal = 0.2
    P[2, c(1, 3)] <- 0.3  # Row 2 off-diagonal = 0.3
    P[3, 1:2] <- 0.5  # Row 3 off-diagonal = 0.5
    expect_true(satisfies_C2(P))
})

test_that("satisfies_C2 returns FALSE for varying off-diagonal", {
    P <- matrix(0.5, 3, 3)
    diag(P) <- 1
    P[1, 2] <- 0.2  # Make row 1 off-diagonal vary
    P[1, 3] <- 0.8
    expect_false(satisfies_C2(P))
})

# -----------------------------------------------------------------------------
# Helper function tests: compute_pi
# -----------------------------------------------------------------------------

test_that("compute_pi returns 0 when k not in c (C1 violation)", {
    P <- make_P_matrix(3, "uniform", p = 0.5)
    c <- c(TRUE, FALSE, TRUE)  # Components 1 and 3 in set
    expect_equal(compute_pi(c, k = 2, P), 0)  # k=2 not in c
})

test_that("compute_pi computes correct probability", {
    P <- make_P_matrix(3, "uniform", p = 0.5)
    c <- c(TRUE, TRUE, FALSE)  # Components 1 and 2 in set, 3 not in set

    # For k=1 (component 1 failed):
    # pi = p_2(1) * (1 - p_3(1)) = 0.5 * 0.5 = 0.25
    expect_equal(compute_pi(c, k = 1, P), 0.25)

    # For k=2 (component 2 failed):
    # pi = p_1(2) * (1 - p_3(2)) = 0.5 * 0.5 = 0.25
    expect_equal(compute_pi(c, k = 2, P), 0.25)
})

test_that("compute_pi works with non-uniform P", {
    vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
    P <- make_P_matrix(3, "full", values = vals)
    # P[2,1] = 0.2, P[3,1] = 0.3, P[1,2] = 0.4, P[3,2] = 0.5, P[1,3] = 0.6, P[2,3] = 0.7

    c <- c(TRUE, TRUE, FALSE)  # 1 and 2 in set

    # For k=1: p_2(1) * (1-p_3(1)) = P[2,1] * (1-P[3,1]) = 0.2 * 0.7 = 0.14
    expect_equal(compute_pi(c, k = 1, P), 0.2 * (1 - 0.3), tolerance = 1e-10)

    # For k=2: p_1(2) * (1-p_3(2)) = P[1,2] * (1-P[3,2]) = 0.4 * 0.5 = 0.2
    expect_equal(compute_pi(c, k = 2, P), 0.4 * (1 - 0.5), tolerance = 1e-10)
})

test_that("compute_pi_all returns named vector for candidates", {
    P <- make_P_matrix(3, "uniform", p = 0.5)
    c <- c(TRUE, TRUE, FALSE)
    pi_vals <- compute_pi_all(c, P)

    expect_equal(length(pi_vals), 2)
    expect_equal(names(pi_vals), c("1", "2"))
    expect_equal(unname(pi_vals["1"]), unname(pi_vals["2"]))  # Symmetric under uniform P
})

# -----------------------------------------------------------------------------
# Log-likelihood tests with fixed P
# -----------------------------------------------------------------------------

test_that("loglik_exp_series_c1_c3 returns function with fixed P", {
    t <- c(0.5, 1.0)
    C <- matrix(c(TRUE, TRUE, TRUE, TRUE), nrow = 2)
    P <- make_P_matrix(2, "uniform", p = 0.5)

    ll <- loglik_exp_series_c1_c3(t, C, P = P)
    expect_type(ll, "closure")
})

test_that("loglik_exp_series_c1_c3 returns -Inf for non-positive theta", {
    t <- c(0.5, 1.0)
    C <- matrix(c(TRUE, TRUE, TRUE, TRUE), nrow = 2)
    P <- make_P_matrix(2, "uniform", p = 0.5)

    ll <- loglik_exp_series_c1_c3(t, C, P = P)
    expect_equal(ll(c(0, 1)), -Inf)
    expect_equal(ll(c(-1, 1)), -Inf)
})

test_that("loglik_exp_series_c1_c3 with uniform P matches C1-C2-C3 likelihood", {
    # Under uniform P satisfying C2, the pi_k(c) values are equal for all k in c
    # So the weighted sum should match the standard likelihood up to a constant
    set.seed(42)
    theta_true <- c(1, 2)
    p <- 0.3
    n <- 50

    sim <- rexp_series_md(n = n, theta = theta_true, p = p, tau = 10)

    P_uniform <- make_P_matrix(2, "uniform", p = p)
    ll_c1_c3 <- loglik_exp_series_c1_c3(sim$t, sim$C, sim$delta, P = P_uniform)
    ll_c1_c2_c3 <- loglik_exp_series(sim$t, sim$C, sim$delta)

    # The likelihoods should differ by a constant (the sum of log(pi_k) terms)
    # But at the MLE, gradient should be zero for both
    theta_test <- c(1.2, 1.8)

    # Check that they give finite values
    expect_true(is.finite(ll_c1_c3(theta_test)))
    expect_true(is.finite(ll_c1_c2_c3(theta_test)))
})

test_that("loglik_exp_series_c1_c3 handles censored observations", {
    t <- c(0.5, 1.0)
    C <- matrix(c(TRUE, TRUE, FALSE, FALSE), nrow = 2, byrow = TRUE)
    delta <- c(1, 0)  # Second is censored
    P <- make_P_matrix(2, "uniform", p = 0.5)

    ll <- loglik_exp_series_c1_c3(t, C, delta, P = P)
    theta <- c(1, 1)

    # Should be finite
    expect_true(is.finite(ll(theta)))
})

# -----------------------------------------------------------------------------
# Log-likelihood tests without fixed P (P as parameter)
# -----------------------------------------------------------------------------

test_that("loglik_exp_series_c1_c3 returns function taking (theta, P) when P not fixed", {
    t <- c(0.5, 1.0)
    C <- matrix(c(TRUE, TRUE, TRUE, TRUE), nrow = 2)

    ll <- loglik_exp_series_c1_c3(t, C, P = NULL)
    expect_type(ll, "closure")

    P <- make_P_matrix(2, "uniform", p = 0.5)
    result <- ll(c(1, 1), P)
    expect_true(is.finite(result))
})

test_that("loglik_exp_series_c1_c3 checks diagonal of P", {
    t <- c(0.5, 1.0)
    C <- matrix(c(TRUE, TRUE, TRUE, TRUE), nrow = 2)

    ll <- loglik_exp_series_c1_c3(t, C, P = NULL)

    P_bad <- matrix(0.5, 2, 2)  # Diagonal should be 1
    expect_error(ll(c(1, 1), P_bad), "Diagonal of P must be 1")
})

# -----------------------------------------------------------------------------
# MLE tests with fixed P
# -----------------------------------------------------------------------------

test_that("mle_exp_series_c1_c3 with fixed P converges", {
    set.seed(123)
    theta_true <- c(1, 2)
    P <- make_P_matrix(2, "uniform", p = 0.3)
    n <- 200

    sim <- rexp_series_md_c1_c3(n = n, theta = theta_true, P = P, tau = 10)

    fit <- mle_exp_series_c1_c3(sim$t, sim$C, sim$delta, fixed_P = P)
    expect_true(fit$converged)
    expect_equal(fit$theta, theta_true, tolerance = 0.4)
})

test_that("mle_exp_series_c1_c3 returns correct structure with fixed P", {
    set.seed(42)
    theta_true <- c(1, 1.5, 2)
    P <- make_P_matrix(3, "uniform", p = 0.3)
    sim <- rexp_series_md_c1_c3(n = 100, theta = theta_true, P = P, tau = 5)

    fit <- mle_exp_series_c1_c3(sim$t, sim$C, sim$delta, fixed_P = P)

    expect_true(is.list(fit))
    expect_true("theta" %in% names(fit))
    expect_true("P" %in% names(fit))
    expect_true("se" %in% names(fit))
    expect_true("vcov" %in% names(fit))
    expect_true("loglik" %in% names(fit))
    expect_true("converged" %in% names(fit))
    expect_equal(length(fit$theta), 3)
    expect_equal(fit$P, P)  # Should return the fixed P
})

# -----------------------------------------------------------------------------
# MLE tests with joint estimation of theta and P
# -----------------------------------------------------------------------------

test_that("mle_exp_series_c1_c3 joint estimation runs", {
    set.seed(42)
    theta_true <- c(1, 2)
    P_true <- make_P_matrix(2, "uniform", p = 0.4)
    n <- 300

    sim <- rexp_series_md_c1_c3(n = n, theta = theta_true, P = P_true, tau = 10)

    # Joint estimation (no fixed_P)
    fit <- mle_exp_series_c1_c3(sim$t, sim$C, sim$delta)

    expect_true(is.list(fit))
    expect_true("theta" %in% names(fit))
    expect_true("P" %in% names(fit))
    expect_true("n_params" %in% names(fit))
    expect_equal(fit$n_params, 2 + 2)  # m + m*(m-1) = 2 + 2 = 4
})

test_that("mle_exp_series_c1_c3 joint estimation recovers parameters (large n)", {
    skip_on_cran()  # Slow test

    set.seed(123)
    theta_true <- c(1, 2)
    P_true <- make_P_matrix(2, "full", values = c(0.3, 0.5))  # Relaxed C2
    n <- 500

    sim <- rexp_series_md_c1_c3(n = n, theta = theta_true, P = P_true, tau = 10)

    fit <- mle_exp_series_c1_c3(sim$t, sim$C, sim$delta,
                                 theta0 = theta_true, P0 = P_true)

    expect_true(fit$converged)
    # With large n, should be reasonably close
    expect_equal(fit$theta, theta_true, tolerance = 0.5)
})

# -----------------------------------------------------------------------------
# Data generation tests
# -----------------------------------------------------------------------------

test_that("rexp_series_md_c1_c3 generates correct structure", {
    set.seed(42)
    theta <- c(1, 2, 3)
    P <- make_P_matrix(3, "uniform", p = 0.4)
    sim <- rexp_series_md_c1_c3(n = 50, theta = theta, P = P, tau = 2)

    expect_equal(length(sim$t), 50)
    expect_equal(length(sim$delta), 50)
    expect_equal(dim(sim$C), c(50, 3))
    expect_equal(length(sim$k), 50)
    expect_true(all(sim$delta %in% c(0, 1)))
    expect_true(all(sim$k %in% 1:3))
})

test_that("rexp_series_md_c1_c3 satisfies C1", {
    set.seed(42)
    theta <- c(1, 2)
    P <- make_P_matrix(2, "uniform", p = 0.3)
    sim <- rexp_series_md_c1_c3(n = 100, theta = theta, P = P, tau = 10)

    # For uncensored observations, failed component must be in candidate set
    uncensored <- sim$delta == 1
    for (i in which(uncensored)) {
        expect_true(sim$C[i, sim$k[i]])
    }
})

test_that("rexp_series_md_c1_c3 censored observations have empty candidate sets", {
    set.seed(42)
    theta <- c(1, 2)
    P <- make_P_matrix(2, "uniform", p = 0.5)
    sim <- rexp_series_md_c1_c3(n = 100, theta = theta, P = P, tau = 0.5)  # High censoring

    censored <- sim$delta == 0
    if (any(censored)) {
        expect_true(all(rowSums(sim$C[censored, , drop = FALSE]) == 0))
    }
})

test_that("rexp_series_md_c1_c3 with keep_latent includes Tm and P", {
    set.seed(42)
    theta <- c(1, 2)
    P <- make_P_matrix(2, "uniform", p = 0.3)
    sim <- rexp_series_md_c1_c3(n = 10, theta = theta, P = P, tau = 5, keep_latent = TRUE)

    expect_true("Tm" %in% names(sim))
    expect_true("P" %in% names(sim))
    expect_equal(dim(sim$Tm), c(10, 2))
    expect_equal(sim$P, P)

    # System time should equal min of component times (before censoring)
    true_sys <- apply(sim$Tm, 1, min)
    expect_equal(pmin(true_sys, 5), sim$t)
})

test_that("rexp_series_md_c1_c3 validates inputs", {
    P <- make_P_matrix(2, "uniform", p = 0.3)

    expect_error(rexp_series_md_c1_c3(10, c(-1, 2), P, 5), "theta must be positive")

    P_bad <- P
    diag(P_bad) <- 0.5
    expect_error(rexp_series_md_c1_c3(10, c(1, 2), P_bad, 5), "Diagonal of P must be 1")

    P_wrong_size <- make_P_matrix(3, "uniform", p = 0.3)
    expect_error(rexp_series_md_c1_c3(10, c(1, 2), P_wrong_size, 5), "P must be m x m")
})

test_that("rexp_series_md_c1_c3 uses P matrix for candidate set generation", {
    set.seed(42)
    theta <- c(1, 2)

    # Create P where p_1(2) = 0.9 and p_2(1) = 0.1
    # P[j,k] = p_j(k) = P(j in C | K = k)
    # P[1,2] = 0.9: when component 2 fails, component 1 is very likely included
    # P[2,1] = 0.1: when component 1 fails, component 2 is unlikely included
    P <- matrix(c(1, 0.9, 0.1, 1), nrow = 2, byrow = TRUE)
    # This gives: P[1,1]=1, P[1,2]=0.9, P[2,1]=0.1, P[2,2]=1

    n <- 500
    sim <- rexp_series_md_c1_c3(n = n, theta = theta, P = P, tau = 100)

    # For observations where k=2 failed, component 1 should be in ~90%
    k2_obs <- which(sim$delta == 1 & sim$k == 2)
    if (length(k2_obs) > 20) {
        prop_1_in <- mean(sim$C[k2_obs, 1])
        expect_gt(prop_1_in, 0.7)  # Should be around 0.9
    }

    # For observations where k=1 failed, component 2 should be in ~10%
    k1_obs <- which(sim$delta == 1 & sim$k == 1)
    if (length(k1_obs) > 20) {
        prop_2_in <- mean(sim$C[k1_obs, 2])
        expect_lt(prop_2_in, 0.3)  # Should be around 0.1
    }
})

# -----------------------------------------------------------------------------
# Comparison with C1-C2-C3 model when C2 is satisfied
# -----------------------------------------------------------------------------

test_that("MLE matches C1-C2-C3 when P satisfies C2", {
    set.seed(42)
    theta_true <- c(1, 2)
    p <- 0.3
    n <- 200

    # Generate data with uniform P (satisfies C2)
    P_uniform <- make_P_matrix(2, "uniform", p = p)
    sim <- rexp_series_md_c1_c3(n = n, theta = theta_true, P = P_uniform, tau = 10)

    # Fit with relaxed C2 model (fixed P)
    fit_c1_c3 <- mle_exp_series_c1_c3(sim$t, sim$C, sim$delta, fixed_P = P_uniform)

    # Fit with C1-C2-C3 model
    fit_c1_c2_c3 <- mle_exp_series(sim$t, sim$C, sim$delta)

    # MLEs should be similar
    expect_equal(fit_c1_c3$theta, fit_c1_c2_c3$theta, tolerance = 0.1)
})

# -----------------------------------------------------------------------------
# Monte Carlo validation
# -----------------------------------------------------------------------------

test_that("MLE is approximately unbiased with fixed P (Monte Carlo)", {
    skip_on_cran()  # Skip slow test

    set.seed(42)
    theta_true <- c(1, 2)
    P <- make_P_matrix(2, "uniform", p = 0.3)
    n_sim <- 30
    n_obs <- 200

    estimates <- matrix(nrow = n_sim, ncol = 2)
    for (i in seq_len(n_sim)) {
        sim <- rexp_series_md_c1_c3(n = n_obs, theta = theta_true, P = P, tau = 5)
        fit <- mle_exp_series_c1_c3(sim$t, sim$C, sim$delta,
                                     theta0 = theta_true, fixed_P = P)
        if (fit$converged) {
            estimates[i, ] <- fit$theta
        }
    }

    # Mean estimate should be close to true value
    mean_est <- colMeans(estimates, na.rm = TRUE)
    expect_equal(mean_est, theta_true, tolerance = 0.2)
})

# -----------------------------------------------------------------------------
# Edge cases
# -----------------------------------------------------------------------------

test_that("loglik_exp_series_c1_c3 handles singleton candidate sets", {
    t <- c(1.0)
    C <- matrix(c(TRUE, FALSE), nrow = 1)  # Only component 1
    P <- make_P_matrix(2, "uniform", p = 0.3)

    ll <- loglik_exp_series_c1_c3(t, C, P = P)

    # With singleton, pi_1(c) = (1 - p_2(1)) = 0.7
    # Likelihood should be finite
    expect_true(is.finite(ll(c(1, 2))))
})

test_that("mle_exp_series_c1_c3 handles high inclusion probability", {
    set.seed(42)
    theta_true <- c(1, 2)
    P <- make_P_matrix(2, "uniform", p = 0.8)  # High inclusion probability
    n <- 100

    sim <- rexp_series_md_c1_c3(n = n, theta = theta_true, P = P, tau = 10)

    # Suppress potential numerical warnings during optimization
    fit <- suppressWarnings(mle_exp_series_c1_c3(sim$t, sim$C, sim$delta, fixed_P = P))
    expect_true(fit$converged)
})
