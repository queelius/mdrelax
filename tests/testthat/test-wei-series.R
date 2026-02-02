# Tests for wei_series.R - Weibull Series System Distribution Functions

# -----------------------------------------------------------------------------
# PDF tests
# -----------------------------------------------------------------------------

test_that("dwei_series returns correct density", {
    # Single component Weibull is just dweibull
    shapes <- c(2)
    scales <- c(3)
    t <- c(0.5, 1, 2, 3, 5)

    for (ti in t) {
        expect_equal(dwei_series(ti, shapes, scales),
                     dweibull(ti, shape = 2, scale = 3),
                     tolerance = 1e-10)
    }
})

test_that("dwei_series with log=TRUE", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    t <- 1.5

    expect_equal(dwei_series(t, shapes, scales, log = TRUE),
                 log(dwei_series(t, shapes, scales)),
                 tolerance = 1e-10)
})

test_that("dwei_series integrates to 1", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)

    result <- integrate(function(t) dwei_series(t, shapes, scales),
                        lower = 0, upper = 50)
    expect_equal(result$value, 1, tolerance = 1e-4)
})

test_that("dwei_series is non-negative", {
    shapes <- c(2, 1.5, 3)
    scales <- c(3, 4, 2)
    t <- seq(0.01, 10, by = 0.1)

    d <- sapply(t, function(ti) dwei_series(ti, shapes, scales))
    expect_true(all(d >= 0))
})

# -----------------------------------------------------------------------------
# CDF tests
# -----------------------------------------------------------------------------

test_that("pwei_series single component matches pweibull", {
    shapes <- c(2)
    scales <- c(3)
    t <- c(0.5, 1, 2, 3, 5)

    expect_equal(pwei_series(t, shapes, scales),
                 pweibull(t, shape = 2, scale = 3),
                 tolerance = 1e-10)
})

test_that("pwei_series is monotonically increasing", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    t <- seq(0, 10, by = 0.5)

    p <- pwei_series(t, shapes, scales)
    expect_true(all(diff(p) >= 0))
})

test_that("pwei_series at 0 is 0, at Inf approaches 1", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)

    expect_equal(pwei_series(0, shapes, scales), 0)
    expect_equal(pwei_series(100, shapes, scales), 1, tolerance = 1e-6)
})

test_that("pwei_series lower.tail=FALSE gives survival", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    t <- 2

    p <- pwei_series(t, shapes, scales)
    s <- pwei_series(t, shapes, scales, lower.tail = FALSE)
    expect_equal(p + s, 1, tolerance = 1e-10)
})

test_that("pwei_series matches surv_wei_series", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    t <- c(0.5, 1, 2, 3)

    p <- pwei_series(t, shapes, scales, lower.tail = FALSE)
    s <- surv_wei_series(t, shapes, scales)
    expect_equal(p, s, tolerance = 1e-10)
})

test_that("pwei_series consistent with dwei_series", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    t <- 2

    # P(T <= t) = integral from 0 to t of f(x)dx
    p_cdf <- pwei_series(t, shapes, scales)
    p_int <- integrate(function(x) sapply(x, function(xi) dwei_series(xi, shapes, scales)),
                        lower = 0, upper = t)$value
    expect_equal(p_cdf, p_int, tolerance = 1e-6)
})

# -----------------------------------------------------------------------------
# Quantile tests
# -----------------------------------------------------------------------------

test_that("qwei_series inverts pwei_series", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    p <- c(0.1, 0.25, 0.5, 0.75, 0.9)

    q <- qwei_series(p, shapes, scales)
    p_back <- pwei_series(q, shapes, scales)
    expect_equal(p_back, p, tolerance = 1e-6)
})

test_that("qwei_series handles boundary values", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)

    expect_equal(qwei_series(0, shapes, scales), 0)
    expect_equal(qwei_series(1, shapes, scales), Inf)
})

# -----------------------------------------------------------------------------
# Random generation tests
# -----------------------------------------------------------------------------

test_that("rwei_series generates correct distribution", {
    set.seed(42)
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    n <- 10000

    x <- rwei_series(n, shapes, scales)

    # Check mean is approximately correct
    mttf <- wei_series_mttf(shapes, scales)
    expect_equal(mean(x), mttf, tolerance = 0.05)

    # Check CDF matches empirical
    t <- 2
    p_theory <- pwei_series(t, shapes, scales)
    p_emp <- mean(x <= t)
    expect_equal(p_emp, p_theory, tolerance = 0.02)
})

# -----------------------------------------------------------------------------
# Hazard function tests
# -----------------------------------------------------------------------------

test_that("hazard_wei_series is sum of component hazards", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    t <- 2

    h <- hazard_wei_series(t, shapes, scales)

    # Manual computation
    h1 <- (shapes[1] / scales[1]) * (t / scales[1])^(shapes[1] - 1)
    h2 <- (shapes[2] / scales[2]) * (t / scales[2])^(shapes[2] - 1)
    expect_equal(h, h1 + h2, tolerance = 1e-10)
})

test_that("hazard_wei_series = f(t)/R(t)", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    t <- 2

    h <- hazard_wei_series(t, shapes, scales)
    f <- dwei_series(t, shapes, scales)
    R <- surv_wei_series(t, shapes, scales)
    expect_equal(h, f / R, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Survival function tests
# -----------------------------------------------------------------------------

test_that("surv_wei_series at 0 is 1", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    expect_equal(surv_wei_series(0, shapes, scales), 1)
})

test_that("surv_wei_series log.p", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    t <- 2

    log_s <- surv_wei_series(t, shapes, scales, log.p = TRUE)
    s <- surv_wei_series(t, shapes, scales)
    expect_equal(log_s, log(s), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# MTTF tests
# -----------------------------------------------------------------------------

test_that("wei_series_mttf is positive", {
    shapes <- c(2, 1.5)
    scales <- c(3, 4)
    mttf <- wei_series_mttf(shapes, scales)
    expect_true(mttf > 0)
})

test_that("wei_series_mttf for single exponential", {
    # Exponential: shape=1, MTTF = scale
    shapes <- c(1)
    scales <- c(5)
    mttf <- wei_series_mttf(shapes, scales)
    expect_equal(mttf, 5, tolerance = 0.01)
})

test_that("wei_series_mttf decreases with more components", {
    shapes <- c(2)
    scales <- c(3)
    mttf1 <- wei_series_mttf(shapes, scales)

    shapes2 <- c(2, 2)
    scales2 <- c(3, 3)
    mttf2 <- wei_series_mttf(shapes2, scales2)

    expect_true(mttf2 < mttf1)
})

# -----------------------------------------------------------------------------
# Exponential as special case of Weibull
# -----------------------------------------------------------------------------

test_that("Weibull with shape=1 is exponential", {
    rate <- 2
    shapes <- c(1)
    scales <- c(1 / rate)
    t <- c(0.5, 1, 2, 3)

    # PDF
    expect_equal(sapply(t, function(ti) dwei_series(ti, shapes, scales)),
                 dexp(t, rate = rate),
                 tolerance = 1e-10)

    # CDF
    expect_equal(pwei_series(t, shapes, scales),
                 pexp(t, rate = rate),
                 tolerance = 1e-10)

    # Survival
    expect_equal(surv_wei_series(t, shapes, scales),
                 pexp(t, rate = rate, lower.tail = FALSE),
                 tolerance = 1e-10)
})

test_that("Weibull series system with shape=1 is exponential series", {
    rates <- c(1, 2, 3)
    shapes <- rep(1, 3)
    scales <- 1 / rates

    t <- c(0.5, 1, 2)

    # System survival = exp(-sum(rates) * t) for exponential
    surv_exp <- exp(-sum(rates) * t)
    surv_wei <- surv_wei_series(t, shapes, scales)
    expect_equal(surv_wei, surv_exp, tolerance = 1e-10)
})
