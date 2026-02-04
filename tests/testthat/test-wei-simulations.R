# =============================================================================
# Tests for Weibull Simulation Framework
# =============================================================================

# -----------------------------------------------------------------------------
# Configuration tests
# -----------------------------------------------------------------------------

test_that("wei_sim_config creates valid configuration", {
  config <- wei_sim_config(n_sim = 50, n_obs = 200,
                            shapes = c(2, 1.5), scales = c(3, 4),
                            tau = 8, seed = 42)

  expect_equal(config$n_sim, 50)
  expect_equal(config$n_obs, 200)
  expect_equal(config$shapes, c(2, 1.5))
  expect_equal(config$scales, c(3, 4))
  expect_equal(config$theta, pack_wei_params(c(2, 1.5), c(3, 4)))
  expect_equal(config$m, 2)
  expect_equal(config$tau, 8)
})

test_that("wei_sim_config validates inputs", {
  expect_error(wei_sim_config(shapes = c(2, 1.5), scales = c(3)),
               "same length")
})

# -----------------------------------------------------------------------------
# Single replication tests
# -----------------------------------------------------------------------------

test_that("run_wei_replication W1 works (baseline)", {
  config <- wei_sim_config(n_sim = 1, n_obs = 100,
                            shapes = c(2, 1.5), scales = c(3, 4))
  set.seed(42)
  result <- run_wei_replication(config, "W1")

  expect_equal(result$scenario, "W1")
  expect_true(is.logical(result$converged))
  expect_true(is.finite(result$loglik) || is.na(result$loglik))
})

test_that("run_wei_replication W2 works (overfit C2)", {
  config <- wei_sim_config(n_sim = 1, n_obs = 100,
                            shapes = c(2, 1.5), scales = c(3, 4))
  set.seed(42)
  result <- run_wei_replication(config, "W2")

  expect_equal(result$scenario, "W2")
  expect_true(is.logical(result$converged))
})

test_that("run_wei_replication W3 works (misspec C2)", {
  config <- wei_sim_config(n_sim = 1, n_obs = 100,
                            shapes = c(2, 1.5), scales = c(3, 4))
  set.seed(42)
  P <- make_P_matrix(2, "full", values = c(0.3, 0.7))
  result <- run_wei_replication(config, "W3", P = P)

  expect_equal(result$scenario, "W3")
  expect_true(is.logical(result$converged))
})

test_that("run_wei_replication W4 works (correct C2)", {
  config <- wei_sim_config(n_sim = 1, n_obs = 100,
                            shapes = c(2, 1.5), scales = c(3, 4))
  set.seed(42)
  P <- make_P_matrix(2, "full", values = c(0.3, 0.7))
  result <- run_wei_replication(config, "W4", P = P)

  expect_equal(result$scenario, "W4")
  expect_true(is.logical(result$converged))
})

test_that("run_wei_replication W5 works (overfit C3)", {
  config <- wei_sim_config(n_sim = 1, n_obs = 100,
                            shapes = c(2, 1.5), scales = c(3, 4))
  set.seed(42)
  result <- run_wei_replication(config, "W5")

  expect_equal(result$scenario, "W5")
  expect_true(is.logical(result$converged))
})

test_that("run_wei_replication W6 works (misspec C3)", {
  config <- wei_sim_config(n_sim = 1, n_obs = 100,
                            shapes = c(2, 1.5), scales = c(3, 4))
  set.seed(42)
  result <- run_wei_replication(config, "W6")

  expect_equal(result$scenario, "W6")
  expect_true(is.logical(result$converged))
})

test_that("run_wei_replication W7 works (correct C3)", {
  config <- wei_sim_config(n_sim = 1, n_obs = 100,
                            shapes = c(2, 1.5), scales = c(3, 4))
  set.seed(42)
  result <- run_wei_replication(config, "W7")

  expect_equal(result$scenario, "W7")
  expect_true(is.logical(result$converged))
})

test_that("run_wei_replication rejects unknown scenario", {
  config <- wei_sim_config()
  expect_error(run_wei_replication(config, "W99"), "Unknown scenario")
})

# -----------------------------------------------------------------------------
# Full simulation study tests
# -----------------------------------------------------------------------------

test_that("run_weibull_simulation_study runs W1 baseline", {
  skip_on_cran()

  study <- run_weibull_simulation_study(
    n_sim = 5, n_obs = 100,
    shapes = c(2, 1.5), scales = c(3, 4),
    scenarios = "W1", verbose = FALSE
  )

  expect_true(is.list(study))
  expect_true("config" %in% names(study))
  expect_true("results" %in% names(study))
  expect_true("summary" %in% names(study))
  expect_equal(nrow(study$results), 5)
})

test_that("run_weibull_simulation_study runs multiple scenarios", {
  skip_on_cran()

  study <- run_weibull_simulation_study(
    n_sim = 3, n_obs = 100,
    shapes = c(2, 1.5), scales = c(3, 4),
    scenarios = c("W1", "W4"), verbose = FALSE
  )

  expect_equal(length(unique(study$results$scenario)), 2)
  expect_true("W1" %in% study$results$scenario)
  expect_true("W4" %in% study$results$scenario)
})

test_that("run_weibull_simulation_study with 3 components", {
  skip_on_cran()

  study <- run_weibull_simulation_study(
    n_sim = 3, n_obs = 100,
    shapes = c(2, 1.5, 2.5), scales = c(3, 4, 2),
    scenarios = "W1", verbose = FALSE
  )

  expect_equal(study$config$m, 3)
  expect_equal(length(study$config$theta), 6)
})

# -----------------------------------------------------------------------------
# Simulation result validation tests
# -----------------------------------------------------------------------------

test_that("W1 baseline produces reasonable estimates", {
  skip_on_cran()

  study <- run_weibull_simulation_study(
    n_sim = 10, n_obs = 300,
    shapes = c(2, 1.5), scales = c(3, 4),
    scenarios = "W1", verbose = FALSE
  )

  # Check convergence rate
  conv_rate <- mean(study$results$converged)
  expect_true(conv_rate > 0.5, info = paste("Convergence rate:", conv_rate))

  # Check summary statistics
  summary_df <- study$summary
  converged <- summary_df$n_converged[1]
  expect_true(converged > 0)
})

test_that("W4 correctly specified relaxed C2 produces reasonable estimates", {
  skip_on_cran()

  P <- make_P_matrix(2, "full", values = c(0.3, 0.7))

  study <- run_weibull_simulation_study(
    n_sim = 10, n_obs = 300,
    shapes = c(2, 1.5), scales = c(3, 4),
    scenarios = "W4", P = P, verbose = FALSE
  )

  conv_rate <- mean(study$results$converged)
  expect_true(conv_rate > 0.5, info = paste("Convergence rate:", conv_rate))
})

# -----------------------------------------------------------------------------
# Summary statistics tests
# -----------------------------------------------------------------------------

test_that("summarize_simulation works with Weibull results", {
  skip_on_cran()

  study <- run_weibull_simulation_study(
    n_sim = 5, n_obs = 100,
    shapes = c(2, 1.5), scales = c(3, 4),
    scenarios = "W1", verbose = FALSE
  )

  summary_df <- study$summary
  expect_true(is.data.frame(summary_df))
  expect_true("bias" %in% names(summary_df))
  expect_true("variance" %in% names(summary_df))
  expect_true("mse" %in% names(summary_df))
})
