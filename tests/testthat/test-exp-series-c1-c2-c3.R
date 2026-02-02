# Tests for exp_series_c1_c2_c3.R - minimal dependency-free implementation

# -----------------------------------------------------------------------------
# Helper function tests
# -----------------------------------------------------------------------------

test_that("decode_matrix extracts prefixed columns correctly", {
    df <- data.frame(x1 = c(T, F), x2 = c(F, T), x3 = c(T, T), other = c(1, 2))
    mat <- decode_matrix(df, "x")
    expect_equal(ncol(mat), 3)
    expect_equal(nrow(mat), 2)
    expect_equal(colnames(mat), c("x1", "x2", "x3"))
    expect_equal(mat[1, ], c(x1 = TRUE, x2 = FALSE, x3 = TRUE))
})

test_that("decode_matrix returns NULL for no matches", {
    df <- data.frame(a = 1, b = 2)
    expect_null(decode_matrix(df, "x"))
})

test_that("encode_matrix creates prefixed columns", {
    mat <- matrix(c(1, 2, 3, 4), nrow = 2)
    df <- encode_matrix(mat, "y")
    expect_equal(names(df), c("y1", "y2"))
    expect_equal(nrow(df), 2)
})

# -----------------------------------------------------------------------------
# Log-likelihood tests
# -----------------------------------------------------------------------------

test_that("loglik_exp_series returns a function", {
    t <- c(0.5, 1.0)
    C <- matrix(c(TRUE, TRUE, TRUE, TRUE), nrow = 2)
    ll <- loglik_exp_series(t, C)
    expect_type(ll, "closure")
})
test_that("loglik_exp_series returns -Inf for non-positive theta", {
    t <- c(0.5, 1.0)
    C <- matrix(c(TRUE, TRUE, TRUE, TRUE), nrow = 2)
    ll <- loglik_exp_series(t, C)
    expect_equal(ll(c(0, 1)), -Inf)
    expect_equal(ll(c(-1, 1)), -Inf)
    expect_equal(ll(c(1, -1)), -Inf)
})

test_that("loglik_exp_series computes correct value", {
    # Simple case: two obs, two components, all in candidate set
    t <- c(1.0, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    ll <- loglik_exp_series(t, C)

    # At theta = (1, 1):
    # -sum(t) * sum(theta) = -2 * 2 = -4
    # + 2 * log(sum(theta)) = 2 * log(2)
    expected <- -4 + 2 * log(2)
    expect_equal(ll(c(1, 1)), expected, tolerance = 1e-10)
})

test_that("loglik_exp_series handles censored observations", {
    t <- c(0.5, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    delta <- c(1, 0)  # Second observation is censored
    ll <- loglik_exp_series(t, C, delta)

    # At theta = (1, 1):
    # Survival: -sum(t) * sum(theta) = -1.5 * 2 = -3
    # Hazard from first only: log(2)
    expected <- -3 + log(2)
    expect_equal(ll(c(1, 1)), expected, tolerance = 1e-10)
})

test_that("loglik_exp_series handles partial candidate sets", {
    t <- c(1.0)
    C <- matrix(c(TRUE, FALSE), nrow = 1)  # Only component 1 in set
    ll <- loglik_exp_series(t, C)

    # At theta = (1, 2):
    # Survival: -1 * 3 = -3
    # Hazard: log(1) = 0 (only theta[1] = 1 in candidate set)
    expected <- -3 + log(1)
    expect_equal(ll(c(1, 2)), expected, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Score tests
# -----------------------------------------------------------------------------

test_that("score_exp_series returns correct gradient", {
    t <- c(1.0, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    sc <- score_exp_series(t, C)

    # At theta = (1, 1):
    # For each j: -sum(t) + n * (1/sum(theta)) = -2 + 2 * (1/2) = -1
    expected <- c(-1, -1)
    expect_equal(sc(c(1, 1)), expected, tolerance = 1e-10)
})

test_that("score matches numerical gradient", {
    t <- c(0.5, 1.2, 0.8)
    C <- matrix(c(TRUE, TRUE, FALSE,
                  TRUE, FALSE, TRUE,
                  FALSE, TRUE, TRUE), nrow = 3, byrow = TRUE)
    ll <- loglik_exp_series(t, C)
    sc <- score_exp_series(t, C)

    theta <- c(1, 1.5, 2)
    h <- 1e-6

    # Numerical gradient
    num_grad <- numeric(3)
    for (j in 1:3) {
        theta_plus <- theta
        theta_plus[j] <- theta_plus[j] + h
        num_grad[j] <- (ll(theta_plus) - ll(theta)) / h
    }

    expect_equal(sc(theta), num_grad, tolerance = 1e-4)
})

# -----------------------------------------------------------------------------
# FIM tests
# -----------------------------------------------------------------------------

test_that("fim_exp_series returns correct FIM", {
    t <- c(1.0, 1.0)
    C <- matrix(TRUE, nrow = 2, ncol = 2)
    fim <- fim_exp_series(t, C)

    # At theta = (1, 1):
    # I[j,k] = n * (1/sum(theta)^2) = 2 * (1/4) = 0.5
    expected <- matrix(0.5, nrow = 2, ncol = 2)
    expect_equal(fim(c(1, 1)), expected, tolerance = 1e-10)
})

test_that("fim_exp_series handles partial candidate sets", {
    t <- c(1.0, 1.0)
    C <- matrix(c(TRUE, FALSE,
                  FALSE, TRUE), nrow = 2, byrow = TRUE)
    fim <- fim_exp_series(t, C)

    # At theta = (1, 1):
    # First obs: C1 = {1}, contributes 1 to I[1,1]
    # Second obs: C2 = {2}, contributes 1 to I[2,2]
    # Off-diagonal: 0
    expected <- diag(c(1, 1))
    expect_equal(fim(c(1, 1)), expected, tolerance = 1e-10)
})

test_that("FIM equals negative Hessian (numerical check)", {
    t <- c(0.5, 1.2, 0.8)
    C <- matrix(c(TRUE, TRUE, FALSE,
                  TRUE, FALSE, TRUE,
                  FALSE, TRUE, TRUE), nrow = 3, byrow = TRUE)
    ll <- loglik_exp_series(t, C)
    fim <- fim_exp_series(t, C)

    theta <- c(1, 1.5, 2)
    h <- 1e-5

    # Numerical Hessian
    H <- matrix(0, 3, 3)
    for (j in 1:3) {
        for (k in 1:3) {
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

    # FIM should equal -H
    expect_equal(fim(theta), -H, tolerance = 1e-3)
})

# -----------------------------------------------------------------------------
# MLE tests
# -----------------------------------------------------------------------------

test_that("mle_exp_series converges on simple data", {
    set.seed(123)
    theta_true <- c(1, 2)
    n <- 200

    # Use rexp_series_md which creates proper partial candidate sets
    # Full candidate sets (C all TRUE) cause non-identifiability
    sim <- rexp_series_md(n = n, theta = theta_true, p = 0.3, tau = 10)

    fit <- mle_exp_series(sim$t, sim$C, sim$delta)
    expect_true(fit$converged)
    expect_equal(fit$theta, theta_true, tolerance = 0.4)
})

test_that("full candidate sets cause non-identifiability", {
    # When all candidate sets are full (all TRUE), only sum(theta) is identifiable
    set.seed(42)
    theta_true <- c(1, 2)
    n <- 200

    Tm <- cbind(rexp(n, rate = theta_true[1]), rexp(n, rate = theta_true[2]))
    t <- apply(Tm, 1, min)
    C <- matrix(TRUE, nrow = n, ncol = 2)  # Full candidate sets

    fit <- mle_exp_series(t, C)
    expect_true(fit$converged)

    # Individual components NOT identifiable, but sum should be correct
    expect_equal(sum(fit$theta), sum(theta_true), tolerance = 0.3)

    # MLE will converge to something like (1.5, 1.5) not (1, 2)
    # because likelihood only depends on sum(theta) when C is all TRUE
})

test_that("mle_exp_series returns correct structure", {
    set.seed(42)
    sim <- rexp_series_md(n = 100, theta = c(1, 1.5, 2), p = 0.3, tau = 5)
    fit <- mle_exp_series(sim$t, sim$C, sim$delta)

    expect_true(is.list(fit))
    expect_true("theta" %in% names(fit))
    expect_true("se" %in% names(fit))
    expect_true("vcov" %in% names(fit))
    expect_true("loglik" %in% names(fit))
    expect_true("converged" %in% names(fit))
    expect_true("fim" %in% names(fit))
    expect_equal(length(fit$theta), 3)
    expect_equal(length(fit$se), 3)
    expect_equal(dim(fit$vcov), c(3, 3))
})

# -----------------------------------------------------------------------------
# Data generation tests
# -----------------------------------------------------------------------------

test_that("rexp_series_md generates correct structure", {
    set.seed(42)
    sim <- rexp_series_md(n = 50, theta = c(1, 2, 3), p = 0.4, tau = 2)

    expect_equal(length(sim$t), 50)
    expect_equal(length(sim$delta), 50)
    expect_equal(dim(sim$C), c(50, 3))
    expect_equal(length(sim$k), 50)
    expect_true(all(sim$delta %in% c(0, 1)))
    expect_true(all(sim$k %in% 1:3))
})

test_that("rexp_series_md satisfies C1 (failed component in candidate set)", {
    set.seed(42)
    sim <- rexp_series_md(n = 100, theta = c(1, 2), p = 0.3, tau = 10)

    # For uncensored observations, failed component must be in candidate set
    uncensored <- sim$delta == 1
    for (i in which(uncensored)) {
        expect_true(sim$C[i, sim$k[i]])
    }
})

test_that("rexp_series_md censored observations have empty candidate sets", {
    set.seed(42)
    sim <- rexp_series_md(n = 100, theta = c(1, 2), p = 0.5, tau = 0.5)  # High censoring

    censored <- sim$delta == 0
    if (any(censored)) {
        expect_true(all(rowSums(sim$C[censored, , drop = FALSE]) == 0))
    }
})

test_that("rexp_series_md with keep_latent includes component times", {
    set.seed(42)
    sim <- rexp_series_md(n = 10, theta = c(1, 2), p = 0.3, tau = 5, keep_latent = TRUE)

    expect_true("Tm" %in% names(sim))
    expect_equal(dim(sim$Tm), c(10, 2))

    # System time should equal min of component times (before censoring)
    true_sys <- apply(sim$Tm, 1, min)
    expect_equal(pmin(true_sys, 5), sim$t)
})

# -----------------------------------------------------------------------------
# Data frame interface tests
# -----------------------------------------------------------------------------

test_that("as_dataframe converts simulation output", {
    set.seed(42)
    sim <- rexp_series_md(n = 5, theta = c(1, 2), p = 0.3, tau = 2)
    df <- as_dataframe(sim)

    expect_true(is.data.frame(df))
    expect_true("t" %in% names(df))
    expect_true("delta" %in% names(df))
    expect_true("x1" %in% names(df))
    expect_true("x2" %in% names(df))
    expect_equal(nrow(df), 5)
})

test_that("loglik_exp_series_df matches direct call", {
    set.seed(42)
    sim <- rexp_series_md(n = 20, theta = c(1, 1.5), p = 0.3, tau = 2)
    df <- as_dataframe(sim)

    ll_direct <- loglik_exp_series(sim$t, sim$C, sim$delta)
    ll_df <- loglik_exp_series_df(df)

    theta <- c(1, 1.5)
    expect_equal(ll_direct(theta), ll_df(theta))
})

test_that("mle_exp_series_df matches direct call", {
    set.seed(42)
    sim <- rexp_series_md(n = 100, theta = c(1, 2), p = 0.3, tau = 3)
    df <- as_dataframe(sim)

    fit_direct <- mle_exp_series(sim$t, sim$C, sim$delta, theta0 = c(1, 1))
    fit_df <- mle_exp_series_df(df, theta0 = c(1, 1))

    expect_equal(fit_direct$theta, fit_df$theta, tolerance = 1e-8)
    expect_equal(fit_direct$loglik, fit_df$loglik, tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
# Monte Carlo validation
# -----------------------------------------------------------------------------

test_that("MLE is approximately unbiased (Monte Carlo)", {
    skip_on_cran()  # Skip slow test on CRAN

    set.seed(42)
    theta_true <- c(1, 2)
    n_sim <- 50
    n_obs <- 200

    estimates <- matrix(nrow = n_sim, ncol = 2)
    for (i in seq_len(n_sim)) {
        sim <- rexp_series_md(n = n_obs, theta = theta_true, p = 0.3, tau = 5)
        fit <- mle_exp_series(sim$t, sim$C, sim$delta, theta0 = theta_true)
        if (fit$converged) {
            estimates[i, ] <- fit$theta
        }
    }

    # Mean estimate should be close to true value
    mean_est <- colMeans(estimates, na.rm = TRUE)
    expect_equal(mean_est, theta_true, tolerance = 0.15)
})
