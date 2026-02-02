#' Exponential Series Masked Data (Dataset 1)
#'
#' Simulated masked failure data from a 3-component exponential series system
#' with rates (1, 1, 2) and Bernoulli masking (p=0.3). 200 observations.
#'
#' @format A data frame with 200 rows and 8 columns:
#' \describe{
#'   \item{t}{System failure time}
#'   \item{t1,t2,t3}{Component failure times (latent)}
#'   \item{k}{True failed component (latent)}
#'   \item{x1,x2,x3}{Candidate set indicators}
#' }
"exp_series_md_1"

#' Exponential Series Masked Data (Dataset 2)
#'
#' Simulated masked failure data from a 3-component exponential series system
#' with right-censoring. 500 observations.
#'
#' @format A data frame with 500 rows and 11 columns:
#' \describe{
#'   \item{t}{System failure time}
#'   \item{t1,t2,t3}{Component failure times (latent)}
#'   \item{k}{True failed component (latent)}
#'   \item{tau}{Right-censoring time}
#'   \item{s}{Observed time (min of t, tau)}
#'   \item{delta}{Censoring indicator (TRUE=observed)}
#'   \item{x1,x2,x3}{Candidate set indicators}
#' }
"exp_series_md_2"

#' Exponential Series MLE Statistics
#'
#' Monte Carlo simulation results for MLE performance across sample sizes.
#'
#' @format A data frame with 50 rows and 15 columns containing bias,
#'   MSE, and standard error estimates for various sample sizes.
"exp_series_stats_1"

#' Guo Weibull Series Masked Data
#'
#' Masked failure data from a 3-component Weibull series system,
#' based on the dataset from Guo et al. (2013).
#'
#' @format A data frame with 30 rows and 4 columns:
#' \describe{
#'   \item{t}{System failure time}
#'   \item{x1,x2,x3}{Candidate set indicators}
#' }
"guo_weibull_series_md"

#' Guo Weibull Series MLE Results
#'
#' Maximum likelihood estimation results for the Guo et al. (2013) dataset.
#'
#' @format A list with components:
#' \describe{
#'   \item{mle}{MLE parameter estimates}
#'   \item{loglike}{Log-likelihood at MLE}
#'   \item{data}{The input data}
#' }
"guo_weibull_series_mle"
