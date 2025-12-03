#' md_bernoulli_cand_C1_kld
#'
#' For each row (observation) in `md`, a probability `Q = (q1, q2, ..., qm)` is
#' constructed such that `qj` is the probability that the j-th component is in
#' the candidate set, `qk = 1`, where `k` is failed component.
#'
#' `Q` is an *informed* candidate model that uses `informative_masking_by_rank`
#' to assign higher probabilities to components that failed earlier (which is
#' something we typically only know in, say, a simulation study).
#'
#' The probabilities `Q` have two constraints on them. Let
#' `P = (p, ..., p, 1, p, ..., p)` be the bernoulli candidate model that
#' satisfies conditions C1, C2, and C3. Then, the KL-divergence between `P`
#' and `Q` is as close as possible to `d` while satisfying `sum(P) == sum(Q)`.
#'
#' For `d = 0`, `Q == P`. As `d` increases, `Q` becomes more informative about
#' the components. Given the structure of `informative_masking_by_rank`, it may
#' not be possible to satisfy every `d` specified, but we get as close as
#' we can, which should permit useful experiments.
#'
#' @param md component failure times for the series system
#' @param p numeric, defines `P = (p, ..., p, 1, p, ..., p)`.
#' @param d numeric, the KL divergence from P = (p, p, ..., p, 1, p, ..., p)
#'          to try to obtain
#' @param debug Logical, whether to output debugging information while running
#' @param eps numeric, stopping condition.
#' @param alpha0 numeric, initial guess for `alpha` parameter of
#'               `informative_masking_by_rank`.
#' @param beta0 numeric, initial guess for `beta` parameter of
#'              `informative_masking_by_rank`.
#' @param lambda numeric, controls how much the two constraints are weighted.
#'               Lower value specifies more enforcement of the KL-divergence
#'               constraint being closer to `d`. Defaults to 1.
#' @param max_iter Integer, maximum number of iterations before giving up.
#' @param lr numeric, learning rate.
#' @importFrom md.tools md_decode_matrix
#' @export
md_bernoulli_cand_C1_kld <- function(
        md, p, d, eps=1e-4, max_iter=100000L, lr=1, lambda=1, alpha0=5, beta0=.5, debug=F)
{
    t <- md_decode_matrix(md,"t")
    m <- ncol(t)
    n <- nrow(t)
    stopifnot(n > 0, m > 0)

    for (i in 1:nrow(md))
    {
        if (debug)
            cat("trial ", i, "\n")

        row <- md[i,]
        tryCatch(expr = {
            res <- generate_relaxed_cand_C1(d=d,
                                            ts=t[i,],
                                            p=p,
                                            debug=debug,
                                            eps=eps,
                                            lr=lr,
                                            max_iter=max_iter,
                                            lambda=lambda,
                                            alpha0=alpha0,
                                            beta0=beta0)
            if (res$converged && abs(sum(res$Q) - sum(res$P)) < .05)
            {
                kl.d <- kl_divergence_bernoulli(res$P,res$Q)
                md[i, "d"] <- kl.d
                md[i, "p"] <- p
                md[i, paste0("q",1:m)] <- t(res$Q)

                if (abs(kl.d - d) > 1e-2)
                {
                    cat("Violation of `KL - d` soft constraint: ",
                        kl.d - d, "\n")
                }
            }
            else
            {
                cat("Failed to converge: stopped at (alpha,beta) = (",
                    res$alpha, ", ", res$beta, ")\n")
            }
        },
        error = function(e) {
            cat("Caught an error: ", e$message, "\n")
        })
    }
    md
}

#' md_bernoulli_cand_C1_C2_C3
#'
#' Bernoulli candidate set model is a particular type of *uninformed* model.
#' Note that we do not generate candidate sets with this function. See
#' `md_cand_sampler` for that.
#'
#' This model satisfies conditions C1, C2, and C3.
#' The failed component will be in the corresponding candidate set with
#' probability 1, and the remaining components will be in the candidate set
#' with probability `p`.
#'
#' @param md masked data.
#' @param p a vector of probabilities
#' @param compvar column name of the component lifetime variables, defaults to
#'                `t`, e.g., `t1, t2, ..., tm`.
#' @param qvar column prefix for component probabilities, defaults to `q`,
#'             e.g., `q1, q2, ..., qm`.
#' @param deltavar column name of the right-censoring indicator variable, 
#'                 defaults to `delta`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @export
md_bernoulli_cand_C1_C2_C3 <- function(md, p,
    compvar = "t",
    qvar = "q",
    deltavar = "delta")
{
    n <- nrow(md)
    if (n == 0) {
        return(md)
    }
    p <- rep(p, length.out = n)
    Tm <- md_decode_matrix(md, compvar)
    if (is.null(Tm)) {
        stop("No component lifetime variables found")
    }
    m <- ncol(Tm)
    Q <- matrix(rep(p, m), nrow = n)
    Q[cbind(1:n, apply(Tm, 1, which.min))] <- 1
    if (!is.null(deltavar) && deltavar %in% colnames(md)) {
        # Set Q=0 for censored observations (delta=0), since there's no failure
        # to generate candidate sets for
        censored <- which(md[[deltavar]] == 0)
        if (length(censored) > 0) {
            Q[censored, ] <- 0
        }
    }
    # remove in case it already has columns for q1,...,qm
    md[ ,paste0(qvar, 1:m)] <- NULL
    md %>% bind_cols(md_encode_matrix(Q, qvar)) %>%
           md_mark_latent(paste0(qvar, 1:m))
}

#' md_cand_sampler
#'
#' Candidate set generator. Requires columns for component probabilities
#' e.g., `q1,...,qm` where `qj` is the probability that the jth component
#' will be in the corresponding candidate set generated for that observation
#' in the `md` table.
#'
#' @param md masked data.
#' @param qvar column prefix for component probabilities, defaults to `q`,
#'             e.g., `q1, q2, ..., qm`.
#' @param setvar column prefix for candidate sets (as Boolean matrix), defaults
#'               to `x`, e.g., `x1, x2, ..., xm`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @importFrom stats runif
#' @export
md_cand_sampler <- function(md,qvar="q",setvar="x")
{
    Q <- md_decode_matrix(md,qvar)
    m <- ncol(Q)
    n <- nrow(Q)
    stopifnot(n > 0, m > 0)

    X <- matrix(rep(NA, m*n), nrow=n)
    for (i in 1:n)
        X[i,] <- runif(m) <= Q[i,]
    md %>% bind_cols(md_encode_matrix(X,setvar))
}

#' md_block_candidate_m3
#'
#' Block candidate model. It produces an MLE that is non-unique, i.e.,
#' if estimating exponential series system, the MLE
#' for λ₁, λ₂, and λ₃ satisfies
#'     `λ̂₁ + λ̂₂ = λ₁ + λ₂`
#' and
#'     `λ̂₃ = λ₃`.
#'
#' This illustrates one of the conditions for the MLE to be unique:
#'
#'     (1) multiple components must not always show together. if they
#'     do, then the best we can learn is that they satisfy some hyper-plane
#'     constraint.
#'
#' @param md masked data
#' @param compvar Prefix-code for the component lifetimes, defaults to `t`,
#'                e.g., `t1, t2, ..., tm`.
#' @param qvar column prefix for component probabilities, defaults to `q`,
#'             e.g., `q1, q2, ..., qm`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @export
md_block_candidate_m3 <- function(md,qvar="q",compvar="t")
{
    Tm <- md_decode_matrix(md,"t")
    stopifnot(!is.null(Tm))
    m <- ncol(Tm)
    n <- nrow(Tm)
    stopifnot(m == 3, n > 0)

    block <- function(ts)
    {
        k <- which.min(ts)
        if (k == 1)
            return(c(1,1,0))
        if (k == 2)
            return(c(1,1,0))
        if (k == 3)
        {
            if (runif(1) < 0.1)
                return(c(1,1,1))
            else
                return(c(0,0,1))
        }
    }


    Q <- matrix(rep(NA,m*n), nrow=n)
    for (i in 1:n)
        Q[i,] <- block(Tm[i,])

    md %>% bind_cols(md_encode_matrix(Q,qvar))
}


# =============================================================================
# Exponential Series System Likelihood Functions under C1, C2, C3
# =============================================================================

#' Log-likelihood function for exponential series system with masked data (C1,C2,C3)
#'
#' Returns a log-likelihood function for an exponential series system with masked
#' component cause of failure under candidate sets satisfying conditions C1, C2,
#' and C3.
#'
#' For a series system with m components having exponential lifetimes with rate
#' parameters `theta = (lambda_1, ..., lambda_m)`:
#' - System reliability: `R(t) = exp(-sum(theta) * t)`
#' - System hazard: `h(t) = sum(theta)`
#'
#' For observation i with system lifetime `s_i`, right-censoring indicator `delta_i`,
#' and candidate set `C_i`:
#' - If uncensored (delta_i = 1): `L_i(theta) = R(s_i) * sum_{j in C_i} h_j(s_i)`
#' - If censored (delta_i = 0): `L_i(theta) = R(s_i)`
#'
#' @param md Masked data frame containing:
#'   - System lifetime column (default: "t")
#'   - Candidate set columns (default: "x1", "x2", ..., "xm")
#'   - Optional right-censoring indicator (default: "delta")
#' @param sysvar Column name for system lifetime. Default is "t".
#' @param setvar Column prefix for candidate set (Boolean matrix). Default is "x".
#' @param deltavar Column name for right-censoring indicator. Default is "delta".
#'   If NULL or column doesn't exist, assumes no censoring.
#'
#' @return A function `f(theta)` that computes the log-likelihood given
#'   rate parameters `theta = (lambda_1, ..., lambda_m)`.
#'
#' @details
#' The log-likelihood is:
#' \deqn{\ell(\theta) = \sum_i \log L_i(\theta)}
#'
#' For uncensored observations:
#' \deqn{\log L_i(\theta) = -s_i \sum_{j=1}^m \lambda_j + \log\left(\sum_{j \in C_i} \lambda_j\right)}
#'
#' For censored observations:
#' \deqn{\log L_i(\theta) = -s_i \sum_{j=1}^m \lambda_j}
#'
#' @examples
#' \dontrun{
#' # Generate some masked data
#' md <- tibble::tibble(
#'   t = c(0.5, 1.2, 0.8),
#'   x1 = c(TRUE, TRUE, FALSE),
#'   x2 = c(TRUE, FALSE, TRUE),
#'   x3 = c(FALSE, TRUE, TRUE),
#'   delta = c(1, 1, 0)  # 1 = observed failure, 0 = censored
#' )
#' ll <- md_loglike_exp_series_C1_C2_C3(md)
#' ll(c(1, 1.5, 2))  # Evaluate at theta = (1, 1.5, 2)
#' }
#'
#' @importFrom md.tools md_decode_matrix
#' @export
md_loglike_exp_series_C1_C2_C3 <- function(md,
                                            sysvar = "t",
                                            setvar = "x",
                                            deltavar = "delta") {
    # Input validation
    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }
    if (!(sysvar %in% colnames(md))) {
        stop("System lifetime column '", sysvar, "' not found in md")
    }

    # Extract data
    s <- md[[sysvar]]  # System lifetimes
    C <- md_decode_matrix(md, setvar)  # Candidate sets as Boolean matrix
    if (is.null(C)) {
        stop("No candidate set columns with prefix '", setvar, "' found")
    }
    m <- ncol(C)

    # Handle censoring indicator
    # Standard survival convention: delta=1 (TRUE) means observed failure (uncensored)
    #                               delta=0 (FALSE) means censored
    has_censoring <- !is.null(deltavar) && deltavar %in% colnames(md)
    if (has_censoring) {
        delta <- as.numeric(md[[deltavar]])  # Convert to numeric for consistent handling
    } else {
        delta <- rep(1, n)  # No censoring = all observed failures
    }

    # Return log-likelihood function
    function(theta) {
        if (length(theta) != m) {
            stop("length(theta) = ", length(theta), " but expected m = ", m)
        }
        if (any(theta <= 0)) {
            return(-Inf)
        }

        # Log-likelihood computation
        sum_theta <- sum(theta)
        ll <- -sum(s) * sum_theta  # Survival contribution for all obs

        # Add hazard contribution for uncensored observations
        for (i in seq_len(n)) {
            if (delta[i] == 1) {  # Uncensored (observed failure)
                hazard_sum <- sum(theta[C[i, ]])
                if (hazard_sum <= 0) {
                    return(-Inf)  # Invalid: no component in candidate set
                }
                ll <- ll + log(hazard_sum)
            }
        }
        ll
    }
}


#' Score function (gradient of log-likelihood) for exponential series system (C1,C2,C3)
#'
#' Returns the score function (gradient of log-likelihood with respect to theta)
#' for an exponential series system with masked data under C1, C2, C3 conditions.
#'
#' @inheritParams md_loglike_exp_series_C1_C2_C3
#'
#' @return A function `g(theta)` that computes the score (gradient) vector
#'   at parameter values `theta = (lambda_1, ..., lambda_m)`.
#'
#' @details
#' The score for component j is:
#' \deqn{\frac{\partial \ell}{\partial \lambda_j} = -\sum_i s_i + \sum_{i: \delta_i=0} \frac{I(j \in C_i)}{\sum_{k \in C_i} \lambda_k}}
#'
#' where `I(j in C_i)` is 1 if component j is in candidate set C_i, 0 otherwise.
#'
#' @examples
#' \dontrun{
#' md <- tibble::tibble(
#'   t = c(0.5, 1.2, 0.8),
#'   x1 = c(TRUE, TRUE, FALSE),
#'   x2 = c(TRUE, FALSE, TRUE),
#'   x3 = c(FALSE, TRUE, TRUE)
#' )
#' score <- md_score_exp_series_C1_C2_C3(md)
#' score(c(1, 1.5, 2))
#' }
#'
#' @importFrom md.tools md_decode_matrix
#' @export
md_score_exp_series_C1_C2_C3 <- function(md,
                                          sysvar = "t",
                                          setvar = "x",
                                          deltavar = "delta") {
    # Input validation
    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }
    if (!(sysvar %in% colnames(md))) {
        stop("System lifetime column '", sysvar, "' not found in md")
    }

    # Extract data
    s <- md[[sysvar]]
    C <- md_decode_matrix(md, setvar)
    if (is.null(C)) {
        stop("No candidate set columns with prefix '", setvar, "' found")
    }
    m <- ncol(C)

    # Handle censoring
    # Standard survival convention: delta=1 (TRUE) means observed failure (uncensored)
    has_censoring <- !is.null(deltavar) && deltavar %in% colnames(md)
    if (has_censoring) {
        delta <- as.numeric(md[[deltavar]])
    } else {
        delta <- rep(1, n)  # No censoring = all observed failures
    }

    # Return score function
    function(theta) {
        if (length(theta) != m) {
            stop("length(theta) = ", length(theta), " but expected m = ", m)
        }
        if (any(theta <= 0)) {
            return(rep(NA_real_, m))
        }

        # Score computation
        # Survival contribution: -sum(s) for each component
        v <- rep(-sum(s), m)

        # Hazard contribution for uncensored observations
        for (i in seq_len(n)) {
            if (delta[i] == 1) {  # Uncensored
                C_i <- C[i, ]
                hazard_sum <- sum(theta[C_i])
                if (hazard_sum <= 0) {
                    return(rep(NA_real_, m))
                }
                # Add 1/sum(theta_C) for each j in C_i
                for (j in seq_len(m)) {
                    if (C_i[j]) {
                        v[j] <- v[j] + 1 / hazard_sum
                    }
                }
            }
        }
        v
    }
}


#' Observed Fisher information matrix for exponential series system (C1,C2,C3)
#'
#' Returns a function that computes the observed Fisher information matrix (FIM)
#' for an exponential series system with masked data under C1, C2, C3 conditions.
#'
#' The observed FIM is the negative Hessian of the log-likelihood, which under
#' regularity conditions converges to the expected FIM and can be inverted to
#' obtain the asymptotic variance-covariance matrix of the MLE.
#'
#' @inheritParams md_loglike_exp_series_C1_C2_C3
#'
#' @return A function `I(theta)` that computes the m x m observed Fisher
#'   information matrix at parameter values `theta`.
#'
#' @details
#' The (j,k) element of the observed FIM is:
#' \deqn{I_{jk}(\theta) = \sum_{i: \delta_i=0} \frac{I(j \in C_i) I(k \in C_i)}{(\sum_{l \in C_i} \lambda_l)^2}}
#'
#' Note: The survival contribution to the Hessian is zero (second derivative
#' of linear term is zero), so only uncensored observations contribute.
#'
#' @examples
#' \dontrun{
#' md <- tibble::tibble(
#'   t = c(0.5, 1.2, 0.8),
#'   x1 = c(TRUE, TRUE, FALSE),
#'   x2 = c(TRUE, FALSE, TRUE),
#'   x3 = c(FALSE, TRUE, TRUE)
#' )
#' fim <- md_fim_exp_series_C1_C2_C3(md)
#' fim(c(1, 1.5, 2))
#' }
#'
#' @importFrom md.tools md_decode_matrix
#' @export
md_fim_exp_series_C1_C2_C3 <- function(md,
                                        sysvar = "t",
                                        setvar = "x",
                                        deltavar = "delta") {
    # Input validation
    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }

    # Extract data
    C <- md_decode_matrix(md, setvar)
    if (is.null(C)) {
        stop("No candidate set columns with prefix '", setvar, "' found")
    }
    m <- ncol(C)

    # Handle censoring
    # Standard survival convention: delta=1 (TRUE) means observed failure (uncensored)
    has_censoring <- !is.null(deltavar) && deltavar %in% colnames(md)
    if (has_censoring) {
        delta <- as.numeric(md[[deltavar]])
    } else {
        delta <- rep(1, n)  # No censoring = all observed failures
    }

    # Return FIM function
    function(theta) {
        if (length(theta) != m) {
            stop("length(theta) = ", length(theta), " but expected m = ", m)
        }
        if (any(theta <= 0)) {
            return(matrix(NA_real_, m, m))
        }

        # FIM computation
        I <- matrix(0, nrow = m, ncol = m)

        for (i in seq_len(n)) {
            if (delta[i] == 1) {  # Only uncensored observations contribute
                C_i <- C[i, ]
                hazard_sum <- sum(theta[C_i])
                if (hazard_sum <= 0) {
                    return(matrix(NA_real_, m, m))
                }
                inv_sq <- 1 / hazard_sum^2

                # Outer product contribution
                for (j in seq_len(m)) {
                    if (C_i[j]) {
                        for (k in seq_len(m)) {
                            if (C_i[k]) {
                                I[j, k] <- I[j, k] + inv_sq
                            }
                        }
                    }
                }
            }
        }
        I
    }
}


#' Maximum likelihood estimator for exponential series system (C1,C2,C3)
#'
#' Computes the MLE of component failure rates for an exponential series
#' system with masked data under C1, C2, C3 candidate set conditions.
#'
#' Uses L-BFGS-B optimization with optional simulated annealing for
#' finding a good starting point.
#'
#' @inheritParams md_loglike_exp_series_C1_C2_C3
#' @param theta0 Initial parameter values for optimization. If NULL, uses
#'   the method of moments estimator based on total system hazard.
#' @param use_annealing Logical. If TRUE, uses simulated annealing to find
#'   a good starting point before local optimization. Default is FALSE.
#' @param control List of control parameters passed to optim(). Default uses
#'   L-BFGS-B with reasonable settings.
#' @param lower Lower bounds for parameters. Default is 1e-10.
#' @param upper Upper bounds for parameters. Default is Inf.
#' @param hessian Logical. If TRUE, compute Hessian at the MLE. Default is TRUE.
#'
#' @return An object of class `mle` (from `algebraic.mle` package) containing:
#'   - `theta.hat`: The MLE of rate parameters
#'   - `sigma`: Estimated variance-covariance matrix
#'   - `loglike`: Log-likelihood at the MLE
#'   - `info`: Observed Fisher information matrix
#'   - `nobs`: Number of observations
#'   - `converged`: Whether optimization converged
#'
#' @details
#' The MLE for exponential series systems with masked data has a closed-form
#' solution only in special cases. In general, numerical optimization is
#' required. This function uses the L-BFGS-B algorithm with analytically
#' computed gradients for efficiency.
#'
#' @examples
#' \dontrun{
#' # Generate masked data
#' library(tibble)
#' set.seed(42)
#' n <- 100
#' theta_true <- c(1, 1.5, 2)
#' # ... (data generation code)
#'
#' mle <- md_mle_exp_series_C1_C2_C3(md, theta0 = theta_true)
#' print(mle)
#' confint(mle)
#' }
#'
#' @importFrom md.tools md_decode_matrix
#' @importFrom algebraic.mle mle_numerical mle
#' @importFrom stats optim
#' @importFrom MASS ginv
#' @export
md_mle_exp_series_C1_C2_C3 <- function(md,
                                        theta0 = NULL,
                                        sysvar = "t",
                                        setvar = "x",
                                        deltavar = "delta",
                                        use_annealing = FALSE,
                                        control = list(),
                                        lower = 1e-10,
                                        upper = Inf,
                                        hessian = TRUE) {
    # Input validation
    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }

    # Extract candidate sets to determine m
    C <- md_decode_matrix(md, setvar)
    if (is.null(C)) {
        stop("No candidate set columns with prefix '", setvar, "' found")
    }
    m <- ncol(C)

    # Get likelihood and score functions
    ll_fn <- md_loglike_exp_series_C1_C2_C3(md, sysvar, setvar, deltavar)
    score_fn <- md_score_exp_series_C1_C2_C3(md, sysvar, setvar, deltavar)
    fim_fn <- md_fim_exp_series_C1_C2_C3(md, sysvar, setvar, deltavar)

    # Initialize theta0 if not provided
    if (is.null(theta0)) {
        # Method of moments: estimate total hazard, divide equally
        s <- md[[sysvar]]
        has_censoring <- !is.null(deltavar) && deltavar %in% colnames(md)
        if (has_censoring) {
            # Standard convention: delta=1 means observed failure (uncensored)
            delta <- as.numeric(md[[deltavar]])
            n_uncensored <- sum(delta == 1)
        } else {
            n_uncensored <- n
        }
        # Total hazard rate estimate: n_uncensored / sum(s)
        total_hazard <- n_uncensored / sum(s)
        theta0 <- rep(total_hazard / m, m)
    }

    # Ensure theta0 is positive
    theta0 <- pmax(theta0, lower)

    # Default control settings
    default_control <- list(
        fnscale = -1,      # Maximize (log-likelihood)
        maxit = 1000,
        factr = 1e7,       # Convergence tolerance
        pgtol = 1e-8
    )
    control <- modifyList(default_control, control)

    # Optional simulated annealing for starting point
    start <- theta0
    if (use_annealing) {
        tryCatch({
            annealing_result <- algebraic.mle:::sim_anneal(
                par = theta0,
                fn = ll_fn,
                control = list(
                    fnscale = -1,
                    maxit = 5000,
                    t_init = 1,
                    t_end = 1e-3,
                    alpha = 0.95,
                    it_per_temp = 10,
                    sup = function(x) all(x > lower)
                )
            )
            start <- annealing_result$par
        }, error = function(e) {
            warning("Simulated annealing failed: ", e$message,
                    ". Using initial theta0.")
        })
    }

    # Optimize using L-BFGS-B
    result <- tryCatch({
        optim(
            par = start,
            fn = ll_fn,
            gr = score_fn,
            method = "L-BFGS-B",
            lower = rep(lower, m),
            upper = rep(upper, m),
            control = control,
            hessian = FALSE  # We compute our own
        )
    }, error = function(e) {
        stop("Optimization failed: ", e$message)
    })

    # Compute observed FIM at MLE
    fim <- fim_fn(result$par)

    # Compute variance-covariance matrix
    sigma <- tryCatch({
        ginv(fim)
    }, error = function(e) {
        warning("Could not invert FIM: ", e$message)
        matrix(NA_real_, m, m)
    })

    # Create MLE object using algebraic.mle interface
    mle(
        theta.hat = result$par,
        loglike = result$value,
        score = score_fn(result$par),
        sigma = sigma,
        info = fim,
        obs = md,
        nobs = n,
        superclasses = c("mle_exp_series_C1_C2_C3", "mle_numerical")
    )
}


#' Generate masked data for exponential series system (C1,C2,C3)
#'
#' Generates synthetic masked failure data for an exponential series system
#' using the Bernoulli candidate set model satisfying C1, C2, C3 conditions.
#'
#' @param n Sample size (number of systems).
#' @param rates Vector of exponential rate parameters for each component.
#' @param p Probability that each non-failed component is in the candidate set.
#' @param tau Right-censoring time. Can be scalar (same for all) or vector.
#'
#' @return A masked data frame with columns:
#'   - t1, t2, ..., tm: Component failure times (latent)
#'   - t: System failure time
#'   - k: Index of failed component (latent)
#'   - delta: Censoring indicator (1 = observed failure, 0 = censored)
#'   - q1, ..., qm: Candidate set probabilities
#'   - x1, ..., xm: Candidate set indicators
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' md <- md_exp_series_system_bernoulli_cand_C1_C2_C3(
#'   n = 100,
#'   rates = c(1, 1.5, 2),
#'   p = 0.3,
#'   tau = 2.0
#' )
#' }
#'
#' @importFrom stats rexp
#' @importFrom md.tools md_encode_matrix md_mark_latent
#' @importFrom tibble tibble
#' @importFrom dplyr %>% bind_cols
#' @export
md_exp_series_system_bernoulli_cand_C1_C2_C3 <- function(n, rates, p, tau) {
    m <- length(rates)
    if (m == 0) {
        stop("rates must have at least one element")
    }
    if (any(rates <= 0)) {
        stop("rates must be positive")
    }

    # Generate component failure times
    Tm <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
        Tm[, j] <- rexp(n, rate = rates[j])
    }

    # System lifetime and failed component
    sys_time <- apply(Tm, 1, min)
    failed_comp <- apply(Tm, 1, which.min)

    # Handle censoring
    # Standard convention: delta=1 means observed failure (uncensored), delta=0 means censored
    tau <- rep(tau, length.out = n)
    delta <- as.numeric(sys_time <= tau)  # 1 if observed failure, 0 if censored
    obs_time <- pmin(sys_time, tau)

    # Create base data frame
    md <- tibble(t = obs_time, k = failed_comp, delta = delta) %>%
        bind_cols(md_encode_matrix(Tm, "t"))

    # Mark component times and k as latent
    md <- md_mark_latent(md, c(paste0("t", 1:m), "k"))

    # Add Bernoulli candidate set model
    md <- md_bernoulli_cand_C1_C2_C3(md, p = p, compvar = "t", deltavar = "delta")

    # Sample candidate sets
    md <- md_cand_sampler(md, qvar = "q", setvar = "x")

    md
}


#' Decorate masked data with series system lifetime and right-censoring
#'
#' Takes component failure times and adds system failure time and
#' right-censoring indicators.
#'
#' @param md Data frame with component failure times (columns t1, t2, ..., tm).
#' @param tau Right-censoring times. Can be scalar or vector of length n.
#' @param compvar Column prefix for component lifetimes. Default is "t".
#' @param sysvar Column name for system lifetime. Default is "t".
#' @param deltavar Column name for censoring indicator. Default is "delta".
#' @param failvar Column name for failed component index. Default is "k".
#'
#' @return Modified data frame with additional columns:
#'   - System lifetime (sysvar)
#'   - Censoring indicator (deltavar): 1 if observed failure, 0 if censored
#'   - Failed component index (failvar)
#'
#' @importFrom md.tools md_decode_matrix md_mark_latent
#' @export
md_series_lifetime_right_censoring <- function(md,
                                                tau,
                                                compvar = "t",
                                                sysvar = "t",
                                                deltavar = "delta",
                                                failvar = "k") {
    n <- nrow(md)
    if (n == 0) {
        return(md)
    }

    # Extract component times
    Tm <- md_decode_matrix(md, compvar)
    if (is.null(Tm)) {
        stop("No component time columns with prefix '", compvar, "' found")
    }
    m <- ncol(Tm)

    # System lifetime is minimum component lifetime
    sys_time <- apply(Tm, 1, min)
    failed_comp <- apply(Tm, 1, which.min)

    # Apply right-censoring
    # Standard convention: delta=1 means observed failure (uncensored), delta=0 means censored
    tau <- rep(tau, length.out = n)
    delta <- as.numeric(sys_time <= tau)  # 1 if observed failure, 0 if censored
    obs_time <- pmin(sys_time, tau)

    # Add columns
    md[[sysvar]] <- obs_time
    md[[deltavar]] <- delta
    md[[failvar]] <- failed_comp

    # Mark component times and failed component as latent
    md <- md_mark_latent(md, c(paste0(compvar, 1:m), failvar))

    md
}


#' Component failure probability under conditions C1 and C2
#'
#' Computes the probability that component j caused system failure
#' given the candidate set C and system failure time t.
#'
#' Under conditions C1 and C2, this probability is:
#' \deqn{P(K=j | C, T=t) = \frac{\lambda_j}{\sum_{k \in C} \lambda_k}}
#'
#' @param theta Rate parameters for components.
#' @param C Candidate set (Boolean vector or matrix).
#' @param j Component index to compute probability for.
#'
#' @return Probability that component j caused the failure.
#'
#' @export
md_series_component_failure_probability_C1_C2 <- function(theta, C, j = NULL) {
    # Handle matrix input (multiple candidate sets)
    if (is.matrix(C)) {
        n <- nrow(C)
        m <- ncol(C)
        if (is.null(j)) {
            # Return all probabilities
            probs <- matrix(0, nrow = n, ncol = m)
            for (i in seq_len(n)) {
                C_i <- C[i, ]
                denom <- sum(theta[C_i])
                if (denom > 0) {
                    for (k in seq_len(m)) {
                        if (C_i[k]) {
                            probs[i, k] <- theta[k] / denom
                        }
                    }
                }
            }
            return(probs)
        } else {
            # Return probability for specific component
            probs <- numeric(n)
            for (i in seq_len(n)) {
                C_i <- C[i, ]
                if (C_i[j]) {
                    probs[i] <- theta[j] / sum(theta[C_i])
                }
            }
            return(probs)
        }
    }

    # Single candidate set (vector)
    if (is.null(j)) {
        m <- length(C)
        probs <- numeric(m)
        denom <- sum(theta[C])
        if (denom > 0) {
            for (k in seq_len(m)) {
                if (C[k]) {
                    probs[k] <- theta[k] / denom
                }
            }
        }
        return(probs)
    } else {
        if (!C[j]) return(0)
        return(theta[j] / sum(theta[C]))
    }
}

