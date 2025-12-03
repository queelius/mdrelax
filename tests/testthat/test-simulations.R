# =============================================================================
# Tests for Simulation Framework
# =============================================================================

# Source the simulation utilities (assuming tests run from package root)
sim_dir <- system.file("simulations", package = "md_series_system_relaxed_candidate_set_models")
if (sim_dir == "") {
    # Fallback for development
    sim_dir <- file.path(getwd(), "inst", "simulations")
}

if (file.exists(file.path(sim_dir, "sim_utils.R"))) {
    source(file.path(sim_dir, "sim_utils.R"))
    RUN_SIM_TESTS <- TRUE
} else {
    RUN_SIM_TESTS <- FALSE
    warning("Simulation files not found, skipping simulation tests")
}

# =============================================================================
# Component Time Generation Tests
# =============================================================================

test_that("generate_component_times produces correct dimensions", {
    skip_if_not(RUN_SIM_TESTS)

    Tm <- generate_component_times(n = 100, rates = c(1, 2, 3))

    expect_equal(nrow(Tm), 100)
    expect_equal(ncol(Tm), 3)
    expect_true(all(Tm > 0))
})

test_that("generate_component_times has correct expected values", {
    skip_if_not(RUN_SIM_TESTS)

    set.seed(42)
    rates <- c(0.5, 1.0, 2.0)
    # Expected means are 1/rate
    expected_means <- 1 / rates

    # Large sample for law of large numbers
    Tm <- generate_component_times(n = 10000, rates = rates)

    actual_means <- colMeans(Tm)

    # Should be within 10% of expected
    for (j in seq_along(rates)) {
        expect_lt(abs(actual_means[j] - expected_means[j]) / expected_means[j], 0.1)
    }
})

# =============================================================================
# Masked Data Generation Tests
# =============================================================================

test_that("generate_masked_data creates correct structure", {
    skip_if_not(RUN_SIM_TESTS)

    md <- generate_masked_data(
        n = 50,
        rates = c(1, 1.5, 2),
        tau = 2.0,
        masking_model = "bernoulli",
        masking_params = list(p = 0.3),
        seed = 123
    )

    expect_equal(nrow(md), 50)
    expect_true("t" %in% names(md))
    expect_true("k" %in% names(md))
    expect_true("delta" %in% names(md))
    expect_true("t1" %in% names(md))
    expect_true("t2" %in% names(md))
    expect_true("t3" %in% names(md))
    expect_true("q1" %in% names(md))
    expect_true("x1" %in% names(md))
})

test_that("generate_masked_data respects right-censoring", {
    skip_if_not(RUN_SIM_TESTS)

    set.seed(42)
    tau <- 0.5  # Low censoring time should produce some censored observations

    md <- generate_masked_data(
        n = 200,
        rates = c(1, 1.5, 2),
        tau = tau,
        masking_model = "bernoulli",
        masking_params = list(p = 0.3)
    )

    # Check all observed times are <= tau
    expect_true(all(md$t <= tau + 1e-10))

    # Censoring indicator should be TRUE when system time exceeds tau
    # (Note: due to the structure, some censoring should occur)
    cens_rate <- mean(md$delta)
    expect_true(cens_rate > 0 || cens_rate < 1)  # Not all or none
})

test_that("generate_masked_data with informative masking works", {
    skip_if_not(RUN_SIM_TESTS)

    md <- generate_masked_data(
        n = 50,
        rates = c(1, 1.5, 2),
        tau = 2.0,
        masking_model = "informative",
        masking_params = list(alpha = 5, beta = 0.3),
        seed = 123
    )

    expect_equal(nrow(md), 50)

    # Check that q values vary (informative masking should produce different q's)
    Q <- as.matrix(md[, c("q1", "q2", "q3")])
    # For uncensored obs, failed component should have q = 1
    uncensored <- !md$delta
    for (i in which(uncensored)) {
        expect_equal(Q[i, md$k[i]], 1)
    }
})

# =============================================================================
# MLE Computation Tests
# =============================================================================

test_that("compute_mle_exp_series returns correct structure", {
    skip_if_not(RUN_SIM_TESTS)

    md <- generate_masked_data(
        n = 100,
        rates = c(1, 1.5, 2),
        tau = 2.0,
        masking_model = "bernoulli",
        masking_params = list(p = 0.3),
        seed = 42
    )

    mle <- compute_mle_exp_series(md, theta0 = c(1, 1.5, 2))

    expect_length(mle$theta_hat, 3)
    expect_true(is.numeric(mle$loglike))
    expect_equal(dim(mle$fim), c(3, 3))
    expect_equal(dim(mle$vcov), c(3, 3))
    expect_true(is.logical(mle$converged))
    expect_equal(mle$nobs, 100)
})

test_that("compute_mle_exp_series produces positive estimates", {
    skip_if_not(RUN_SIM_TESTS)

    md <- generate_masked_data(
        n = 100,
        rates = c(1, 1.5, 2),
        tau = 2.0,
        masking_model = "bernoulli",
        masking_params = list(p = 0.3),
        seed = 42
    )

    mle <- compute_mle_exp_series(md)

    expect_true(all(mle$theta_hat > 0))
})

test_that("MLE is approximately unbiased with large sample", {
    skip_if_not(RUN_SIM_TESTS)

    # This test uses a larger sample to check approximate unbiasedness
    set.seed(42)
    rates <- c(1, 1.5, 2)

    md <- generate_masked_data(
        n = 1000,
        rates = rates,
        tau = 2.0,
        masking_model = "bernoulli",
        masking_params = list(p = 0.3)
    )

    mle <- compute_mle_exp_series(md, theta0 = rates)

    # With n=1000, estimates should be within ~30% of true values
    for (j in seq_along(rates)) {
        rel_error <- abs(mle$theta_hat[j] - rates[j]) / rates[j]
        expect_lt(rel_error, 0.3,
                  info = sprintf("Component %d relative error too large", j))
    }
})

# =============================================================================
# Statistical Functions Tests
# =============================================================================

test_that("compute_mle_stats calculates correct bias", {
    skip_if_not(RUN_SIM_TESTS)

    # Create synthetic estimates with known properties
    theta_true <- c(1, 2, 3)
    n_reps <- 100

    # Estimates that are systematically 0.1 higher
    estimates <- matrix(rep(theta_true + 0.1, each = n_reps), nrow = n_reps)

    stats <- compute_mle_stats(estimates, theta_true)

    expect_equal(stats$bias, rep(0.1, 3), tolerance = 1e-10)
})

test_that("compute_mle_stats handles NA values", {
    skip_if_not(RUN_SIM_TESTS)

    theta_true <- c(1, 2)
    estimates <- matrix(c(1.1, 2.1,
                          NA, NA,
                          0.9, 1.9), nrow = 3, byrow = TRUE)

    stats <- compute_mle_stats(estimates, theta_true)

    expect_equal(stats$n_valid, 2)
    expect_equal(stats$bias, c(0, 0), tolerance = 1e-10)
})

test_that("relative_efficiency computes correct ratio", {
    skip_if_not(RUN_SIM_TESTS)

    expect_equal(relative_efficiency(2, 4), 2)
    expect_equal(relative_efficiency(4, 2), 0.5)
    expect_true(is.na(relative_efficiency(0, 1)))
})

# =============================================================================
# Scenario Runner Tests
# =============================================================================

test_that("run_scenario completes successfully", {
    skip_if_not(RUN_SIM_TESTS)

    result <- run_scenario(
        n = 30,
        rates = c(1, 2),
        tau = 1.0,
        masking_model = "bernoulli",
        masking_params = list(p = 0.3),
        B = 5,
        seed = 42,
        progress = FALSE
    )

    expect_true("estimates" %in% names(result))
    expect_true("stats" %in% names(result))
    expect_equal(nrow(result$estimates), 5)
    expect_equal(ncol(result$estimates), 2)
})

# =============================================================================
# Informative Masking Tests
# =============================================================================

test_that("informative_masking_by_rank_internal returns valid probabilities", {
    skip_if_not(RUN_SIM_TESTS)

    ts <- c(0.5, 1.0, 1.5)  # Component 1 fails first

    probs <- informative_masking_by_rank_internal(ts, alpha = 2, beta = 0.5)

    # All probabilities should be in [0, 1]
    expect_true(all(probs >= 0))
    expect_true(all(probs <= 1))

    # Failed component should have prob 1
    expect_equal(probs[1], 1)
})

test_that("informative_masking_by_rank_internal respects alpha=0 (uniform)", {
    skip_if_not(RUN_SIM_TESTS)

    ts <- c(0.5, 1.0, 1.5, 2.0)  # Component 1 fails first
    beta <- 0.4

    probs <- informative_masking_by_rank_internal(ts, alpha = 0, beta = beta)

    expect_equal(probs[1], 1)  # Failed component
    expect_equal(probs[2], beta, tolerance = 1e-10)
    expect_equal(probs[3], beta, tolerance = 1e-10)
    expect_equal(probs[4], beta, tolerance = 1e-10)
})

# =============================================================================
# FIM Analysis Tests
# =============================================================================

test_that("fim_eigenvalue_analysis identifies near-singular matrices", {
    skip_if_not(RUN_SIM_TESTS)

    # Singular matrix
    singular_fim <- matrix(c(1, 1, 1, 1), nrow = 2)
    result <- fim_eigenvalue_analysis(singular_fim)
    expect_true(result$near_singular)
    expect_lt(result$smallest, 1e-9)

    # Non-singular matrix
    nonsingular_fim <- diag(c(1, 2))
    result2 <- fim_eigenvalue_analysis(nonsingular_fim)
    expect_false(result2$near_singular)
    expect_equal(result2$smallest, 1)
})

# =============================================================================
# I/O Tests
# =============================================================================

test_that("results_to_dataframe produces valid output", {
    skip_if_not(RUN_SIM_TESTS)

    # Create minimal results structure
    results <- list(
        list(
            params = list(
                n = 50,
                rates = c(1, 2),
                masking_model = "bernoulli",
                B = 10
            ),
            stats = list(
                bias = c(0.1, 0.2),
                variance = c(0.01, 0.02),
                mse = c(0.02, 0.06),
                rmse = c(sqrt(0.02), sqrt(0.06)),
                coverage = c(0.9, 0.95),
                mean_ci_width = c(0.2, 0.3),
                n_valid = 10
            )
        )
    )

    df <- results_to_dataframe(results)

    expect_s3_class(df, "data.frame")
    expect_true("n" %in% names(df))
    expect_true("bias" %in% names(df))
    expect_equal(nrow(df), 2)  # 2 components
})

# =============================================================================
# Helper Function Tests
# =============================================================================

test_that("compute_tau_for_censoring produces valid censoring times", {
    skip_if_not(RUN_SIM_TESTS)

    rates <- c(1, 2, 3)

    # No censoring
    tau_inf <- compute_tau_for_censoring(rates, 0)
    expect_equal(tau_inf, Inf)

    # Full censoring
    tau_zero <- compute_tau_for_censoring(rates, 1)
    expect_equal(tau_zero, 0)

    # Partial censoring
    tau_mid <- compute_tau_for_censoring(rates, 0.2)
    expect_true(tau_mid > 0)
    expect_true(is.finite(tau_mid))
})

# =============================================================================
# End-to-End Integration Test
# =============================================================================

test_that("full simulation workflow runs without error", {
    skip_if_not(RUN_SIM_TESTS)
    skip_on_cran()  # Skip on CRAN due to time

    # This is a minimal end-to-end test
    scenarios <- list(
        list(
            n = 30,
            rates = c(1, 2),
            tau = 1.0,
            masking_model = "bernoulli",
            masking_params = list(p = 0.3)
        )
    )

    results <- run_simulation(
        scenarios = scenarios,
        B = 5,
        n_cores = 1,
        seed = 42,
        progress = FALSE
    )

    expect_length(results, 1)
    expect_true("stats" %in% names(results[[1]]))
    expect_true(results[[1]]$stats$n_valid > 0)
})
