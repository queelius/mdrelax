# Tests for exponential series system likelihood functions under C1, C2, C3

test_that("md_loglike_exp_series_C1_C2_C3 returns a function", {
  md <- tibble::tibble(
    t = c(0.5, 1.2, 0.8),
    x1 = c(TRUE, TRUE, FALSE),
    x2 = c(TRUE, FALSE, TRUE),
    x3 = c(FALSE, TRUE, TRUE)
  )
  ll <- md_loglike_exp_series_C1_C2_C3(md)
  expect_type(ll, "closure")
})

test_that("md_loglike_exp_series_C1_C2_C3 returns -Inf for non-positive theta", {
  md <- tibble::tibble(
    t = c(0.5, 1.0),
    x1 = c(TRUE, TRUE),
    x2 = c(TRUE, FALSE)
  )
  ll <- md_loglike_exp_series_C1_C2_C3(md)
  expect_equal(ll(c(0, 1)), -Inf)
  expect_equal(ll(c(-1, 1)), -Inf)
  expect_equal(ll(c(1, -1)), -Inf)
})

test_that("md_loglike_exp_series_C1_C2_C3 computes correct log-likelihood", {
  # Simple case: two observations, two components, all in candidate set
  md <- tibble::tibble(
    t = c(1.0, 1.0),  # Both have lifetime 1
    x1 = c(TRUE, TRUE),
    x2 = c(TRUE, TRUE)
  )
  ll <- md_loglike_exp_series_C1_C2_C3(md)

  # At theta = (1, 1):
  # -sum(t) * sum(theta) = -2 * 2 = -4
  # + 2 * log(sum(theta)) = 2 * log(2)
  expected <- -4 + 2 * log(2)
  expect_equal(ll(c(1, 1)), expected, tolerance = 1e-10)
})

test_that("md_loglike_exp_series_C1_C2_C3 handles censored observations", {
  md <- tibble::tibble(
    t = c(0.5, 1.0),
    x1 = c(TRUE, TRUE),
    x2 = c(TRUE, TRUE),
    delta = c(FALSE, TRUE)  # Second observation is censored
  )
  ll <- md_loglike_exp_series_C1_C2_C3(md)

  # At theta = (1, 1):
  # Survival: -sum(t) * sum(theta) = -1.5 * 2 = -3
  # Hazard contribution only from first (uncensored): log(sum(theta)) = log(2)
  expected <- -3 + log(2)
  expect_equal(ll(c(1, 1)), expected, tolerance = 1e-10)
})

test_that("md_score_exp_series_C1_C2_C3 returns a function", {
  md <- tibble::tibble(
    t = c(0.5, 1.2, 0.8),
    x1 = c(TRUE, TRUE, FALSE),
    x2 = c(TRUE, FALSE, TRUE),
    x3 = c(FALSE, TRUE, TRUE)
  )
  score <- md_score_exp_series_C1_C2_C3(md)
  expect_type(score, "closure")
})

test_that("md_score_exp_series_C1_C2_C3 computes correct gradient", {
  md <- tibble::tibble(
    t = c(1.0, 1.0),
    x1 = c(TRUE, TRUE),
    x2 = c(TRUE, TRUE)
  )
  score <- md_score_exp_series_C1_C2_C3(md)

  # At theta = (1, 1):
  # For each j: -sum(t) + n * (1/sum(theta)) = -2 + 2 * (1/2) = -1
  expected <- c(-1, -1)
  expect_equal(score(c(1, 1)), expected, tolerance = 1e-10)
})

test_that("md_fim_exp_series_C1_C2_C3 returns correct FIM", {
  md <- tibble::tibble(
    t = c(1.0, 1.0),
    x1 = c(TRUE, TRUE),
    x2 = c(TRUE, TRUE)
  )
  fim <- md_fim_exp_series_C1_C2_C3(md)

  # At theta = (1, 1):
  # I[j,k] = n * (1/sum(theta)^2) = 2 * (1/4) = 0.5
  expected <- matrix(0.5, nrow = 2, ncol = 2)
  result <- fim(c(1, 1))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("md_fim_exp_series_C1_C2_C3 handles partial candidate sets", {
  md <- tibble::tibble(
    t = c(1.0, 1.0),
    x1 = c(TRUE, FALSE),  # First obs: only comp 1
    x2 = c(FALSE, TRUE)   # Second obs: only comp 2
  )
  fim <- md_fim_exp_series_C1_C2_C3(md)

  # At theta = (1, 1):
  # First obs: C1 = {1}, contributes (1/1)^2 = 1 to I[1,1]
  # Second obs: C2 = {2}, contributes (1/1)^2 = 1 to I[2,2]
  # Off-diagonal: 0 (no overlap)
  expected <- diag(c(1, 1))
  result <- fim(c(1, 1))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("score is gradient of log-likelihood (numerical check)", {
  md <- tibble::tibble(
    t = c(0.5, 1.2, 0.8),
    x1 = c(TRUE, TRUE, FALSE),
    x2 = c(TRUE, FALSE, TRUE),
    x3 = c(FALSE, TRUE, TRUE)
  )
  ll <- md_loglike_exp_series_C1_C2_C3(md)
  score <- md_score_exp_series_C1_C2_C3(md)

  theta <- c(1, 1.5, 2)
  h <- 1e-6

  # Numerical gradient
  num_grad <- numeric(3)
  for (j in 1:3) {
    theta_plus <- theta
    theta_plus[j] <- theta_plus[j] + h
    num_grad[j] <- (ll(theta_plus) - ll(theta)) / h
  }

  analytic_grad <- score(theta)
  expect_equal(analytic_grad, num_grad, tolerance = 1e-4)
})

test_that("grad_descent converges on simple quadratic", {
  # f(x) = (x - 2)^2, minimum at x = 2
  f <- function(x) (x - 2)^2
  result <- grad_descent(f, x0 = 0, lr = 0.1, eps = 1e-6, max_iter = 1000)
  expect_true(result$converged)
  expect_equal(result$param, 2, tolerance = 1e-3)
})

test_that("grad_descent respects support constraint", {
  # f(x) = x^2, but constrained to x >= 1
  f <- function(x) x^2
  sup <- function(x) x >= 1
  result <- grad_descent(f, x0 = 5, sup = sup, lr = 0.1, eps = 1e-6, max_iter = 1000)
  # Minimum in support is at x = 1
  expect_true(result$param >= 1 - 1e-6)
})

test_that("informative_masking_by_rank uses correct parameter names", {
  # Test that the function works with alpha and beta (not alpha0/beta0)
  ts <- c(0.5, 1.0, 1.5)  # component 1 fails first
  probs <- informative_masking_by_rank(ts, alpha = 0, beta = 0.5)

  # With alpha = 0, all non-failed components have same probability beta
  # Failed component (rank 1) has prob 1
  # probs should be c(1, 0.5, 0.5) permuted by order of ts
  expect_equal(probs[1], 1)  # Component 1 failed first
  expect_equal(probs[2], 0.5, tolerance = 1e-10)
  expect_equal(probs[3], 0.5, tolerance = 1e-10)
})

test_that("informative_masking_by_rank validates parameters",
{
  ts <- c(0.5, 1.0, 1.5)
  expect_error(informative_masking_by_rank(ts, alpha = -1, beta = 0.5))
  expect_error(informative_masking_by_rank(ts, alpha = 0, beta = -0.1))
  expect_error(informative_masking_by_rank(ts, alpha = 0, beta = 1.1))
})
