# Tests for wei_series_c1_c3.R - Weibull Masked Data (Relaxed C2)
# Informative masking where P[j,k] = P(j in C | K = k) can vary with k

# -----------------------------------------------------------------------------
# Function existence tests (RED phase - will fail until implemented)
# -----------------------------------------------------------------------------

test_that("wei_series_c1_c3 functions exist", {
  expect_true(exists("loglik_wei_series_c1_c3"))
  expect_true(exists("score_wei_series_c1_c3"))
  expect_true(exists("fim_wei_series_c1_c3"))
  expect_true(exists("mle_wei_series_c1_c3"))
  expect_true(exists("rwei_series_md_c1_c3"))
})

# -----------------------------------------------------------------------------
# Log-likelihood tests
# -----------------------------------------------------------------------------

test_that("loglik_wei_series_c1_c3 returns a function", {
  t <- c(0.5, 1.0)
  C <- matrix(TRUE, nrow = 2, ncol = 2)
  P <- make_P_matrix(2, "uniform", p = 0.5)
  ll <- loglik_wei_series_c1_c3(t, C, P = P)
  expect_type(ll, "closure")
})

test_that("loglik_wei_series_c1_c3 returns -Inf for non-positive params", {
  t <- c(0.5, 1.0)
  C <- matrix(TRUE, nrow = 2, ncol = 2)
  P <- make_P_matrix(2, "uniform", p = 0.5)
  ll <- loglik_wei_series_c1_c3(t, C, P = P)
  expect_equal(ll(c(0, 1, 1, 1)), -Inf)
  expect_equal(ll(c(-1, 1, 1, 1)), -Inf)
  expect_equal(ll(c(1, 0, 1, 1)), -Inf)
})

test_that("loglik_wei_series_c1_c3 with uniform P matches wei_series", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  p <- 0.3
  P_uniform <- make_P_matrix(2, "uniform", p = p)

  sim <- rwei_series_md(n = 50, shapes = shapes, scales = scales, p = p, tau = 5)

  ll_baseline <- loglik_wei_series(sim$t, sim$C, sim$delta)
  ll_c13 <- loglik_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P_uniform)

  theta <- pack_wei_params(shapes, scales)

  # Both likelihoods should be finite at the same theta
  expect_true(is.finite(ll_baseline(theta)))
  expect_true(is.finite(ll_c13(theta)))

  # Gradients at MLE should both be near zero
})

test_that("loglik_wei_series_c1_c3 handles non-uniform P", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  # For m=2, "full" needs m*(m-1) = 2 off-diagonal values
  P <- make_P_matrix(2, "full", values = c(0.3, 0.7))

  sim <- rwei_series_md_c1_c3(n = 50, shapes = shapes, scales = scales, P = P, tau = 5)

  ll <- loglik_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)
  theta <- pack_wei_params(shapes, scales)

  val <- ll(theta)
  expect_true(is.finite(val))
  expect_true(val < 0)  # Log-likelihood should be negative
})

test_that("loglik_wei_series_c1_c3 handles censored observations", {
  t <- c(0.5, 1.0)
  C <- matrix(TRUE, nrow = 2, ncol = 2)
  C[2, ] <- FALSE  # Censored observation has empty candidate set
  delta <- c(1, 0)
  P <- make_P_matrix(2, "uniform", p = 0.5)
  ll <- loglik_wei_series_c1_c3(t, C, delta, P = P)

  theta <- c(2, 3, 1.5, 4)
  val <- ll(theta)
  expect_true(is.finite(val))
})

test_that("loglik_wei_series_c1_c3 computes correct value (manual)", {
  # Single observation, single component, uncensored
  t <- c(2)
  C <- matrix(TRUE, nrow = 1, ncol = 1)
  P <- matrix(1, nrow = 1, ncol = 1)  # Trivial case
  ll <- loglik_wei_series_c1_c3(t, C, P = P)

  k <- 2; b <- 3
  theta <- c(k, b)

  # ll = -(t/b)^k + log(h(t) * pi_1(c))
  # For single component, pi_1(c) = 1, h(t) = (k/b) * (t/b)^(k-1)
  z <- t / b
  expected <- -z^k + log((k / b) * z^(k - 1))
  # Use unname to strip names attribute
  expect_equal(unname(ll(theta)), expected[[1]], tolerance = 1e-10)
})

test_that("loglik_wei_series_c1_c3 with 3 components", {
  set.seed(42)
  shapes <- c(2, 1.5, 2.5)
  scales <- c(3, 4, 2)
  P <- make_P_matrix(3, "full", values = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7))

  sim <- rwei_series_md_c1_c3(n = 100, shapes = shapes, scales = scales, P = P, tau = 5)

  ll <- loglik_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)
  theta <- pack_wei_params(shapes, scales)

  val <- ll(theta)
  expect_true(is.finite(val))
})

# -----------------------------------------------------------------------------
# Score tests
# -----------------------------------------------------------------------------

test_that("score_wei_series_c1_c3 matches numerical gradient", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  P <- make_P_matrix(2, "full", values = c(0.3, 0.7))

  sim <- rwei_series_md_c1_c3(n = 30, shapes = shapes, scales = scales, P = P, tau = 5)

  ll <- loglik_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)
  sc <- score_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)

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

test_that("score_wei_series_c1_c3 handles non-true parameter values", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  P <- make_P_matrix(2, "uniform", p = 0.4)

  sim <- rwei_series_md_c1_c3(n = 30, shapes = shapes, scales = scales, P = P, tau = 5)
  sc <- score_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)

  # Test at various parameter values
  theta1 <- c(1.5, 2.5, 2, 5)
  theta2 <- c(3, 4, 1, 3)

  g1 <- sc(theta1)
  g2 <- sc(theta2)

  expect_equal(length(g1), 4)
  expect_equal(length(g2), 4)
  expect_true(all(is.finite(g1)))
  expect_true(all(is.finite(g2)))
})

test_that("score_wei_series_c1_c3 with 3 components matches numerical", {
  set.seed(42)
  shapes <- c(2, 1.5, 2.5)
  scales <- c(3, 4, 2)
  P <- make_P_matrix(3, "full", values = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7))

  sim <- rwei_series_md_c1_c3(n = 30, shapes = shapes, scales = scales, P = P, tau = 5)

  ll <- loglik_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)
  sc <- score_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)

  theta <- pack_wei_params(shapes, scales)
  h <- 1e-6
  num_grad <- sapply(seq_along(theta), function(j) {
    theta_plus <- theta
    theta_plus[j] <- theta_plus[j] + h
    (ll(theta_plus) - ll(theta)) / h
  })

  expect_equal(sc(theta), unname(num_grad), tolerance = 1e-3)
})

# -----------------------------------------------------------------------------
# FIM tests
# -----------------------------------------------------------------------------

test_that("fim_wei_series_c1_c3 equals negative Hessian", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  P <- make_P_matrix(2, "uniform", p = 0.4)

  sim <- rwei_series_md_c1_c3(n = 30, shapes = shapes, scales = scales, P = P, tau = 5)

  ll <- loglik_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)
  fim <- fim_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)

  theta <- c(2, 3, 1.5, 4)
  h <- 1e-5

  # Numerical Hessian (4-point method)
  p <- 4
  H <- matrix(0, p, p)
  for (j in 1:p) {
    for (k in 1:p) {
      theta_pp <- theta_pm <- theta_mp <- theta_mm <- theta
      theta_pp[j] <- theta_pp[j] + h; theta_pp[k] <- theta_pp[k] + h
      theta_pm[j] <- theta_pm[j] + h; theta_pm[k] <- theta_pm[k] - h
      theta_mp[j] <- theta_mp[j] - h; theta_mp[k] <- theta_mp[k] + h
      theta_mm[j] <- theta_mm[j] - h; theta_mm[k] <- theta_mm[k] - h
      H[j, k] <- (ll(theta_pp) - ll(theta_pm) -
                 ll(theta_mp) + ll(theta_mm)) / (4 * h^2)
    }
  }

  I <- fim(theta)
  expect_equal(I, -H, tolerance = 1e-2)
})

test_that("fim_wei_series_c1_c3 is positive definite at reasonable theta", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  P <- make_P_matrix(2, "uniform", p = 0.4)

  sim <- rwei_series_md_c1_c3(n = 100, shapes = shapes, scales = scales, P = P, tau = 5)
  fim <- fim_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)

  theta <- pack_wei_params(shapes, scales)
  I <- fim(theta)

  # Check positive definiteness via eigenvalues
  eig <- eigen(I, symmetric = TRUE)$values
  expect_true(all(eig > 0),
              info = paste("Eigenvalues:", paste(round(eig, 4), collapse=", ")))
})

# -----------------------------------------------------------------------------
# MLE tests
# -----------------------------------------------------------------------------

test_that("mle_wei_series_c1_c3 converges", {
  set.seed(42)
  shapes_true <- c(2, 1.5)
  scales_true <- c(3, 4)
  P <- make_P_matrix(2, "uniform", p = 0.3)

  sim <- rwei_series_md_c1_c3(n = 200, shapes = shapes_true, scales = scales_true,
                               P = P, tau = 8)

  fit <- mle_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)
  expect_true(fit$converged)
  expect_equal(length(fit$theta), 4)
  expect_true(all(fit$theta > 0))
})

test_that("mle_wei_series_c1_c3 returns correct structure", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  P <- make_P_matrix(2, "uniform", p = 0.3)

  sim <- rwei_series_md_c1_c3(n = 100, shapes = shapes, scales = scales,
                               P = P, tau = 8)
  fit <- mle_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)

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

test_that("mle_wei_series_c1_c3 recovers true parameters approximately", {
  set.seed(42)
  shapes_true <- c(2, 1.5)
  scales_true <- c(3, 4)
  theta_true <- pack_wei_params(shapes_true, scales_true)
  P <- make_P_matrix(2, "uniform", p = 0.3)

  sim <- rwei_series_md_c1_c3(n = 500, shapes = shapes_true, scales = scales_true,
                               P = P, tau = 8)
  fit <- mle_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P,
                               theta0 = theta_true * 0.8)

  if (fit$converged) {
    expect_equal(fit$shapes, shapes_true, tolerance = 0.5)
    expect_equal(fit$scales, scales_true, tolerance = 1.0)
  }
})

test_that("mle_wei_series_c1_c3 with non-uniform P", {
  set.seed(42)
  shapes_true <- c(2, 1.5)
  scales_true <- c(3, 4)
  theta_true <- pack_wei_params(shapes_true, scales_true)
  P <- make_P_matrix(2, "full", values = c(0.3, 0.7))

  sim <- rwei_series_md_c1_c3(n = 500, shapes = shapes_true, scales = scales_true,
                               P = P, tau = 8)
  fit <- mle_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P,
                               theta0 = theta_true * 0.8)

  expect_true(fit$converged)
  if (fit$converged) {
    expect_equal(fit$shapes, shapes_true, tolerance = 0.5)
    expect_equal(fit$scales, scales_true, tolerance = 1.0)
  }
})

test_that("mle_wei_series_c1_c3 with 3 components", {
  set.seed(42)
  shapes_true <- c(2, 1.5, 2.5)
  scales_true <- c(3, 4, 2)
  theta_true <- pack_wei_params(shapes_true, scales_true)
  P <- make_P_matrix(3, "uniform", p = 0.4)

  sim <- rwei_series_md_c1_c3(n = 500, shapes = shapes_true, scales = scales_true,
                               P = P, tau = 5)
  fit <- mle_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P,
                               theta0 = theta_true * 0.8)

  expect_true(fit$converged)
  expect_equal(fit$m, 3)
})

# -----------------------------------------------------------------------------
# Data generation tests
# -----------------------------------------------------------------------------

test_that("rwei_series_md_c1_c3 generates correct structure", {
  set.seed(42)
  shapes <- c(2, 1.5, 3)
  scales <- c(3, 4, 2)
  P <- make_P_matrix(3, "uniform", p = 0.4)

  sim <- rwei_series_md_c1_c3(n = 50, shapes = shapes, scales = scales, P = P, tau = 5)

  expect_equal(length(sim$t), 50)
  expect_equal(length(sim$delta), 50)
  expect_equal(dim(sim$C), c(50, 3))
  expect_equal(length(sim$k), 50)
  expect_true(all(sim$delta %in% c(0, 1)))
  expect_true(all(sim$k %in% 1:3))
})

test_that("rwei_series_md_c1_c3 satisfies C1", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  P <- make_P_matrix(2, "full", values = c(0.3, 0.7))

  sim <- rwei_series_md_c1_c3(n = 100, shapes = shapes, scales = scales, P = P, tau = 10)

  uncensored <- sim$delta == 1
  for (i in which(uncensored)) {
    expect_true(sim$C[i, sim$k[i]])
  }
})

test_that("rwei_series_md_c1_c3 censored obs have empty candidate sets", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  P <- make_P_matrix(2, "uniform", p = 0.5)

  sim <- rwei_series_md_c1_c3(n = 100, shapes = shapes, scales = scales, P = P, tau = 0.5)

  censored <- sim$delta == 0
  if (any(censored)) {
    expect_true(all(rowSums(sim$C[censored, , drop = FALSE]) == 0))
  }
})

test_that("rwei_series_md_c1_c3 with keep_latent includes component times", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  P <- make_P_matrix(2, "uniform", p = 0.3)

  sim <- rwei_series_md_c1_c3(n = 10, shapes = shapes, scales = scales,
                               P = P, tau = 10, keep_latent = TRUE)

  expect_true("Tm" %in% names(sim))
  expect_equal(dim(sim$Tm), c(10, 2))

  # System time = min of components (before censoring)
  true_sys <- apply(sim$Tm, 1, min)
  expect_equal(pmin(true_sys, 10), sim$t, tolerance = 1e-10)
})

test_that("rwei_series_md_c1_c3 masking varies with k for non-uniform P", {
  set.seed(123)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  # P[j,k] = P(j in C | K = k), stored by column
  # We want: P[2,1] = 0.2 (low chance comp 2 in C when comp 1 fails)
  #          P[1,2] = 0.8 (high chance comp 1 in C when comp 2 fails)
  # The matrix is:
  # [1, P[1,2]]   [1,   0.8]
  # [P[2,1], 1] = [0.2, 1  ]
  P <- matrix(c(1, 0.2, 0.8, 1), nrow = 2, byrow = FALSE)

  sim <- rwei_series_md_c1_c3(n = 1000, shapes = shapes, scales = scales, P = P, tau = 10)

  # When component 1 fails (K=1), component 2's inclusion depends on P[2,1] = 0.2
  # When component 2 fails (K=2), component 1's inclusion depends on P[1,2] = 0.8
  k1_obs <- which(sim$k == 1 & sim$delta == 1)
  k2_obs <- which(sim$k == 2 & sim$delta == 1)

  if (length(k1_obs) > 50 && length(k2_obs) > 50) {
    prop_2_in_c_given_k1 <- mean(sim$C[k1_obs, 2])  # Should be ~0.2
    prop_1_in_c_given_k2 <- mean(sim$C[k2_obs, 1])  # Should be ~0.8

    expect_true(prop_2_in_c_given_k1 < 0.35,
                info = paste("Prop 2 in C given K=1:", round(prop_2_in_c_given_k1, 3)))
    expect_true(prop_1_in_c_given_k2 > 0.65,
                info = paste("Prop 1 in C given K=2:", round(prop_1_in_c_given_k2, 3)))
  }
})

# -----------------------------------------------------------------------------
# Data frame interface tests
# -----------------------------------------------------------------------------

test_that("loglik_wei_series_c1_c3_df matches direct call", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  P <- make_P_matrix(2, "uniform", p = 0.3)

  sim <- rwei_series_md_c1_c3(n = 20, shapes = shapes, scales = scales, P = P, tau = 5)
  df <- as_dataframe(sim)

  ll_direct <- loglik_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P)
  ll_df <- loglik_wei_series_c1_c3_df(df, P = P)

  theta <- pack_wei_params(shapes, scales)
  expect_equal(ll_direct(theta), ll_df(theta))
})

test_that("mle_wei_series_c1_c3_df matches direct call", {
  set.seed(42)
  shapes <- c(2, 1.5)
  scales <- c(3, 4)
  theta_true <- pack_wei_params(shapes, scales)
  P <- make_P_matrix(2, "uniform", p = 0.3)

  sim <- rwei_series_md_c1_c3(n = 100, shapes = shapes, scales = scales, P = P, tau = 5)
  df <- as_dataframe(sim)

  fit_direct <- mle_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P, theta0 = theta_true * 0.9)
  fit_df <- mle_wei_series_c1_c3_df(df, P = P, theta0 = theta_true * 0.9)

  expect_equal(fit_direct$theta, fit_df$theta)
  expect_equal(fit_direct$loglik, fit_df$loglik)
})

# -----------------------------------------------------------------------------
# Monte Carlo tests
# -----------------------------------------------------------------------------

test_that("Weibull C1-C3 MLE is approximately unbiased (Monte Carlo)", {
  skip_on_cran()

  set.seed(42)
  shapes_true <- c(2, 1.5)
  scales_true <- c(3, 4)
  theta_true <- pack_wei_params(shapes_true, scales_true)
  P <- make_P_matrix(2, "uniform", p = 0.3)
  n_sim <- 20
  n_obs <- 300

  estimates <- matrix(nrow = n_sim, ncol = 4)
  for (i in seq_len(n_sim)) {
    sim <- rwei_series_md_c1_c3(n = n_obs, shapes = shapes_true,
                                 scales = scales_true, P = P, tau = 8)
    fit <- mle_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P,
                                 theta0 = theta_true * 0.9)
    if (fit$converged) {
      estimates[i, ] <- fit$theta
    }
  }

  mean_est <- colMeans(estimates, na.rm = TRUE)
  expect_equal(mean_est[1], shapes_true[1], tolerance = 0.5)
  expect_equal(mean_est[3], shapes_true[2], tolerance = 0.5)
})

test_that("Weibull C1-C3 MLE with non-uniform P is unbiased (Monte Carlo)", {
  skip_on_cran()

  set.seed(42)
  shapes_true <- c(2, 1.5)
  scales_true <- c(3, 4)
  theta_true <- pack_wei_params(shapes_true, scales_true)
  P <- make_P_matrix(2, "full", values = c(0.3, 0.7))
  n_sim <- 20
  n_obs <- 300

  estimates <- matrix(nrow = n_sim, ncol = 4)
  for (i in seq_len(n_sim)) {
    sim <- rwei_series_md_c1_c3(n = n_obs, shapes = shapes_true,
                                 scales = scales_true, P = P, tau = 8)
    fit <- mle_wei_series_c1_c3(sim$t, sim$C, sim$delta, P = P,
                                 theta0 = theta_true * 0.9)
    if (fit$converged) {
      estimates[i, ] <- fit$theta
    }
  }

  mean_est <- colMeans(estimates, na.rm = TRUE)
  bias <- mean_est - theta_true
  # Bias should be small for correctly specified model
  expect_true(max(abs(bias)) < 0.5,
              info = paste("Bias:", paste(round(bias, 3), collapse=", ")))
})
