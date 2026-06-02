test_that("derivatives of loglikelihood are correct", {
    data <- generate_test_data_1()$data
    
    nbasis <- 5
    k <- 2
    sp <- 10
    basis <- find_orthogonal_spline_basis(nbasis, data$x)

    alpha_components <- find_alpha_components(nbasis, k)

    set.seed(1)
    par <- c(rnorm(nbasis + length(alpha_components), sd = 0.1), 1)

    lp_grad <- loglikelihood_pen_grad(par, basis$X, data$y, data$c - 1, sp, basis$S, k, 1)
    lp_hess <- loglikelihood_pen_hess(par, basis$X, data$y, data$c - 1, sp, basis$S, k, 1)

    lp_fun <- function(par) {
        loglikelihood_pen(par, basis$X, data$y, data$c - 1, sp, basis$S, k, 1)
    }

    lp_grad_man <- numDeriv::grad(lp_fun, par)
    lp_hess_man <- numDeriv::hessian(lp_fun, par)

    expect_equal(lp_grad, lp_grad_man)
    expect_equal(lp_hess, lp_hess_man)

    #library(microbenchmark)
    #microbenchmark(
    #    loglikelihood_pen(par, basis$X, data$y, data$c - 1, sp, basis$S, k, 1),
    #    loglikelihood_pen_grad(par, basis$X, data$y, data$c - 1, sp, basis$S, k, 1),
    #    loglikelihood_pen_hess(par, basis$X, data$y, data$c - 1, sp, basis$S, k, 1)
    #)
    #' timings similar to passing in beta: doing transform does not add too much cost
    #' can do ~ 900 iterations per second with the gradient



})

test_that("same beta with different psi_index gives same penalised log-likelihood", {
    data <- generate_test_data_1()$data
    
    nbasis <- 5
    k <- 2
    sp <- 10
    basis <- find_orthogonal_spline_basis(nbasis, data$x)

    alpha_components <- find_alpha_components(nbasis, k)

    set.seed(1)
    lsigma <- 1
    par <- c(rnorm(nbasis + length(alpha_components), sd = 0.1), lsigma)
    beta0 <- par[1:nbasis]
    alpha_1 <- par[-c((1:nbasis), length(par))]
    beta <- find_beta(alpha_1, nbasis, k, 1)
    alpha_2 <- find_alpha_from_beta(beta, nbasis, k, 2)
    beta_from_alpha_2 <- find_beta(alpha_2, nbasis, k, 2)

    expect_equal(beta, beta_from_alpha_2)

    theta_1 <- c(beta0, alpha_1, lsigma)
    theta_2 <- c(beta0, alpha_2, lsigma)

    sp <- 1
    l_1 <- loglikelihood_pen(theta_1, basis$X, data$y, data$c - 1, sp, basis$S, k, psi_index = 1)
    l_2 <- loglikelihood_pen(theta_2, basis$X, data$y, data$c - 1, sp, basis$S, k, psi_index = 2)

    expect_equal(l_1, l_2)

})


test_that("checking for discontinuities", {

    nbasis <- 3
    K <- 2

    X <- diag(nbasis)
    y <- c(1, 0, 1)
    cluster <- rep(0L, length(y))

    S <- matrix(0, nbasis, nbasis)
    sp <- 0

    make_alpha <- function(t) {
        c(t, 1, 0,   # alpha_1
          1, 1)      # alpha_2
    }
    
    theta_at_t <- function(t, sigma = 0.2) {
        beta0 <- rep(0, nbasis)
        alpha <- make_alpha(t)
        c(beta0, alpha, log(sigma))
    }

    l <- function(theta) {
        loglikelihood_pen(
            theta,
            X = X,
            y = y,
            c = cluster,
            sp = sp,
            S = S,
            K = K,
            psi_index = 1
        )
    }
    
    ll <- function(t) {
        l(theta_at_t(t))
    }

    t_grid <- c(-1e-3, -1e-4, -1e-6, -1e-8, 1e-8, 1e-6, 1e-4, 1e-3)

    l_grid <- sapply(t_grid, ll)
    plot(t_grid, l_grid, type = "l")

    theta0_pos <- theta_at_t(1e-8)
    theta0_neg <- theta_at_t(-1e-8)

    opt_pos <- optim(theta0_pos, l, method = "BFGS", control = list(fnscale = -1))
    opt_neg <- optim(theta0_neg, l, method = "BFGS", control = list(fnscale = -1))
    #' get convergence to different values

    theta_neg <- opt_neg$par
    
    par_split_neg <- split_par(theta_neg, nbasis)

    diag_neg <- householder_diagnostic(
        alpha = par_split_neg$alpha,
        nbasis = nbasis,
        k = K
    )
    #' problematic!
    expect_lt(diag_neg$min_ratio, 1e-5)
    
    theta_pos <- opt_pos$par
    
    par_split_pos <- split_par(theta_pos, nbasis)

    diag_pos <- householder_diagnostic(
        alpha = par_split_pos$alpha,
        nbasis = nbasis,
        k = K
    )
    #' fine
    expect_gt(diag_pos$min_ratio, 1e-5)
    
})
