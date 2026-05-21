add_hessian_and_log_ml <- function(fit, basis, data) {
    fit$hessian <- loglikelihood_pen_hess(fit$par,
                                          X = basis$X, y = data$y, c = data$c - 1,
                                          sp = fit$sp, S = basis$S, K = fit$k,
                                          alpha_index = fit$alpha_index)
    fit$var_par <- tryCatch(solve(-fit$hessian),
                            error = function(cond) {
                                pracma::pinv(-fit$hessian)
                            })
                            
    fit$log_ml <- approx_log_ml(fit, fit$hessian, basis)
    fit
}



split_par <- function(par, nbasis) {
    components <- rep("alpha", length(par))
    components[1:nbasis] <- "beta0"
    components[length(par)] <- "lsigma"
    split(par, components)
}

find_par_cluster <- function(beta0, beta, u_hat) {
    beta0 + tcrossprod(beta, u_hat)
}


find_fit_info <- function(opt, k, basis, sp, data, alpha_index) {
    par <- opt$par
    l_pen <- opt$value
    par_split <- split_par(par, basis$nbasis)


    f0_x <- basis$X %*% par_split$beta0
    f0 <- find_spline_fun(par_split$beta0, basis)
    
    if(k > 0) {
        beta <- find_beta(par_split$alpha, basis$nbasis, k, alpha_index)
        lambda <- colSums(beta^2)
        f_x <- basis$X %*% beta
        u_hat_full <- find_u_hat(exp(par_split$lsigma), data, f0_x, f_x, var = TRUE)
        u_hat <- u_hat_full$u_hat
        var_u_hat <- u_hat_full$var_u_hat
        f <- find_spline_fun(beta, basis)
        par_cluster <- find_par_cluster(par_split$beta0, beta, u_hat)
    } else {
        f <- NULL
        f_x <- matrix(nrow = length(data$x), ncol = 0)
        u_hat <- NULL
        var_u_hat <- NULL
        beta <- NULL
        lambda <- NULL
        par_cluster <- NULL
    }

    list(k = k,
         sp = sp,
         par = par,
         alpha_index = alpha_index,
         l_pen = l_pen,
         opt = opt,
         beta0 = par_split$beta0,
         alpha = par_split$alpha,
         lsigma = par_split$lsigma,
         beta = beta,
         sigma = exp(par_split$lsigma),
         par_cluster = par_cluster,
         f0 = f0,
         f0_x = f0_x,
         f = f,
         f_x = f_x,
         u_hat = u_hat,
         var_u_hat = var_u_hat,
         lambda = lambda,
         data = data,
         basis = basis)    
}

fit_given_par0 <- function(data, sp, k, par0, basis, alpha_index) {
     opt <- stats::optim(par0, loglikelihood_pen, loglikelihood_pen_grad,
                         X = basis$X, y = data$y, c = data$c - 1,
                         sp = sp, S = basis$S, K = k, alpha_index = alpha_index,
                         method = "BFGS", control = list(fnscale = -1, maxit = 10000))
    if(opt$convergence != 0)
        warning("optim has not converged")
     fit <- find_fit_info(opt, k, basis, sp, data, alpha_index)
    
     fit
}

fit_given_par0_nlm <- function(data, sp, k, par0, basis, alpha_index) {
    ml <- function(par) {
        result <- -loglikelihood_pen(par, X = basis$X, y = data$y,
                                     c = data$c - 1,
                                     sp = sp, S = basis$S,
                                     K = k, alpha_index = alpha_index)
        gr <- -loglikelihood_pen_grad(par,
                                      X = basis$X, y = data$y, c = data$c - 1,
                                      sp = sp, S = basis$S, K = k,
                                      alpha_index = alpha_index)
        attr(result, "gradient") <- gr
        result
    }
    opt_nlm <- stats::nlm(ml, par0, iterlim = 100000, check.analyticals = FALSE)

    # convert to format from optim
    list(par = opt_nlm$estimate,
         value = -opt_nlm$minimum,
         convergence = ifelse(opt_nlm$code < 2, 0, opt_nlm$code))
                
}


# Fit the mean-only model
fit_0 <- function(data, sp, basis, alpha_index = 1) {
    X_0 <- basis$X
    
    Xt_y <- crossprod(X_0, data$y)
    XtX <- crossprod(X_0, X_0)

    
    beta_0 <- as.numeric(solve(XtX + sp * basis$S, Xt_y))
    y_hat_0 <- X_0 %*% beta_0
    resid <- data$y - y_hat_0
    sigma <- stats::sd(resid)

    par0 <- c(beta_0, log(sigma))
    
    fit_given_par0(data, sp, 0, par0, basis, alpha_index)
}



find_par0_given_fit_km1 <- function(fit_km1, k, nbasis, fit_k_other_sp = NULL) {
    if(k == 0) {
        if(is.null(fit_k_other_sp))
            return(c(rep(0.01, nbasis), 0))
        else
            return(fit_k_other_sp$par)
    }
    if(is.null(fit_k_other_sp))
        alpha_k0 <- rep(0.01, nbasis - k + 1)
    else
        alpha_k0 <- fit_k_other_sp$alpha[find_alpha_components(nbasis, k) == k]

    c(fit_km1$beta0, fit_km1$alpha, alpha_k0, fit_km1$lsigma)
}


# fit model with k eigenfunctions, given model fit with k - 1 eigenfunctions, same sp
fit_given_fit_km1 <- function(data, sp, k, fit_km1, basis, fit_k_other_sp = NULL,
                              alpha_index) {
    if(!is.null(fit_km1$alpha_index) && fit_km1$alpha_index != alpha_index) {
        stop("fit_km1 was fitted with a different alpha_index")
    }

    if(!is.null(fit_k_other_sp) &&
       !is.null(fit_k_other_sp$alpha_index) &&
       fit_k_other_sp$alpha_index != alpha_index) {
        stop("fit_k_other_sp was fitted with a different alpha_index")
    }
    
    par0 <- find_par0_given_fit_km1(fit_km1, k, basis$nbasis, fit_k_other_sp)        
    fit_given_par0(data, sp, k, par0, basis, alpha_index)
   
}


find_FVE <- function(mod) {
    cumsum(mod$lambda) / sum(mod$lambda)
}


is_k_larger_than_required <- function(mod, k_tol) {
    FVE <-  find_FVE(mod)
    FVE[length(FVE) - 1] > 1 - k_tol
}


fits_given_sp <- function(sp, kmax, data, basis, k_tol, fits_other_sp = NULL,
                          alpha_index) {
    fits <- list(fit_0(data, sp, basis, alpha_index))
    for(k in 1:kmax) {
        if(length(fits_other_sp) > k)
            fit_k_other_sp <- fits_other_sp[[k+1]]
        else
            fit_k_other_sp <- NULL
        
        fits[[k+1]] <- fit_given_fit_km1(data, sp, k, fits[[k]], basis,
                                         fit_k_other_sp, alpha_index)
        if(k > 1)
            if(is_k_larger_than_required(fits[[k+1]], k_tol))
                break
    }
    fits
}
