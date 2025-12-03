## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages, message=FALSE---------------------------------------------
library(mdrelax)

## ----simulate-data------------------------------------------------------------
set.seed(42)

# True parameters
lambda_true <- c(1.0, 1.5, 2.0)
n <- 100  # sample size
m <- 3    # number of components

# Generate component failure times
T_mat <- matrix(rexp(n * m, rate = rep(lambda_true, each = n)),
                nrow = n, ncol = m)

# System failure time is minimum of component times
t_sys <- apply(T_mat, 1, min)

# Identify which component failed (latent)
k_failed <- apply(T_mat, 1, which.min)

# Create data frame with component times
md <- data.frame(
  t = t_sys,
  t1 = T_mat[, 1],
  t2 = T_mat[, 2],
  t3 = T_mat[, 3],
  delta = rep(1, n)  # all observed (no censoring): delta=1 means failure observed
)

## ----apply-masking------------------------------------------------------------
# Apply Bernoulli masking with p = 0.3
md <- md_bernoulli_cand_C1_C2_C3(md, p = 0.3)

# Sample candidate sets
md <- md_cand_sampler(md)

# View first few rows
head(md[, c("t", "x1", "x2", "x3")])

## ----mle----------------------------------------------------------------------
fit <- md_mle_exp_series_C1_C2_C3(md)
print(fit)

