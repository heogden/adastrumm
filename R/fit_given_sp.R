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

reparameterise_fit_without_optim <- function(fit, basis, data, sp, alpha_index_new) {
    par_new <- par_from_beta_parameterisation(
        beta0 = fit$beta0,
        beta = fit$beta,
        lsigma = fit$lsigma,
        nbasis = basis$nbasis,
        k = fit$k,
        alpha_index = alpha_index_new
    )

    opt_new <- list(
        par = par_new,
        value = fit$l_pen,
        convergence = fit$opt$convergence,
        counts = fit$opt$counts,
        message = "Reparameterised without reoptimisation"
    )

    fit_new <- find_fit_info(
        opt = opt_new,
        k = fit$k,
        basis = basis,
        sp = sp,
        data = data,
        alpha_index = alpha_index_new
    )

    add_hessian_and_log_ml(fit_new, basis, data)
}

choose_alpha_index_from_beta_ci <- function(fit, basis, data, sp) {
    candidates <- seq_len(basis$nbasis - fit$k + 1)

    mods <- lapply(candidates, function(alpha_index) {
        reparameterise_fit_without_optim(
            fit = fit,
            basis = basis,
            data = data,
            sp = sp,
            alpha_index_new = alpha_index
        )
    })

    scores <- vapply(mods, function(mod) {
        if(!is_neg_def(mod$hessian) ||
           matrixcalc::is.singular.matrix(mod$hessian)) {
            return(-Inf)
        }

        householder_ci_diagnostic(mod)$min_z_to_boundary
    }, numeric(1))


    if(all(!is.finite(scores))) {
        return(list(
            alpha_index = fit$alpha_index,
            mod = fit,
            scores = scores,
            candidates = candidates
        ))
    }
    
    best <- which.max(scores)

    list(
        alpha_index = candidates[best],
        mod = mods[[best]],
        scores = scores,
        candidates = candidates
    )
}

maybe_reparameterise_after_hessian <- function(fit, data, sp, basis,
                                               alpha_tol = 1e-2,
                                               alpha_ci_tol = 2,
                                               auto_alpha = TRUE) {
    if(!auto_alpha || fit$k <= 1) {
        return(fit)
    }

    bad_hessian <- !is_neg_def(fit$hessian) ||
        matrixcalc::is.singular.matrix(fit$hessian)

    point_diag <- householder_diagnostic(
        fit$alpha,
        nbasis = basis$nbasis,
        k = fit$k,
        alpha_index = fit$alpha_index
    )

    bad_point_chart <- point_diag$min_ratio < alpha_tol

    ## If the fitted point/Hessian is bad, reoptimise in a better chart.
    if(bad_hessian || bad_point_chart) {
        alpha_index_new <- choose_alpha_index_from_beta(
            fit$beta,
            nbasis = basis$nbasis,
            k = fit$k
        )

        if(alpha_index_new == fit$alpha_index) {
            return(fit)
        }

        par0_new <- par_from_beta_parameterisation(
            beta0 = fit$beta0,
            beta = fit$beta,
            lsigma = fit$lsigma,
            nbasis = basis$nbasis,
            k = fit$k,
            alpha_index = alpha_index_new
        )

        fit_new <- fit_given_par0(
            data = data,
            sp = sp,
            k = fit$k,
            par0 = par0_new,
            basis = basis,
            alpha_index = alpha_index_new
        )

        return(add_hessian_and_log_ml(fit_new, basis, data))
    }

    ## Only use CI diagnostic when the Hessian is already usable.
    ci_diag <- householder_ci_diagnostic(fit)
    bad_ci_chart <- ci_diag$min_z_to_boundary < alpha_ci_tol

    if(!bad_ci_chart) {
        fit$alpha_ci_diagnostic <- ci_diag
        return(fit)
    }

    choice <- choose_alpha_index_from_beta_ci(
        fit = fit,
        basis = basis,
        data = data,
        sp = sp
    )

    if(choice$alpha_index == fit$alpha_index) {
        fit$alpha_ci_scores <- choice$scores
        fit$alpha_ci_candidates <- choice$candidates
        fit$alpha_ci_diagnostic <- ci_diag
        return(fit)
    }

    choice$mod
}


maybe_reparameterise_and_refit <- function(fit, data, sp, basis,
                                           alpha_tol = 1e-2,
                                           auto_alpha = TRUE) {
    if(!auto_alpha || fit$k <= 1) {
        return(fit)
    }

    diag <- householder_diagnostic(
        fit$alpha,
        nbasis = basis$nbasis,
        k = fit$k,
        alpha_index = fit$alpha_index
    )

    if(diag$min_ratio >= alpha_tol) {
        return(fit)
    }

    alpha_index_new <- choose_alpha_index_from_beta(
        fit$beta,
        nbasis = basis$nbasis,
        k = fit$k
    )

    if(alpha_index_new == fit$alpha_index) {
        return(fit)
    }

    par0_new <- par_from_beta_parameterisation(
        beta0 = fit$beta0,
        beta = fit$beta,
        lsigma = fit$lsigma,
        nbasis = basis$nbasis,
        k = fit$k,
        alpha_index = alpha_index_new
    )

    fit_given_par0(
        data = data,
        sp = sp,
        k = fit$k,
        par0 = par0_new,
        basis = basis,
        alpha_index = alpha_index_new
    )
}

fit_given_start_beta <- function(data, sp, k,
                                 beta0, beta, lsigma,
                                 basis,
                                 alpha_index = 1,
                                 auto_alpha = TRUE,
                                 alpha_tol = 1e-2) {
    start <- maybe_switch_alpha_index_start(
        beta0 = beta0,
        beta = beta,
        lsigma = lsigma,
        basis = basis,
        k = k,
        alpha_index = alpha_index,
        alpha_tol = alpha_tol,
        auto_alpha = auto_alpha
    )

    fit <- fit_given_par0(
        data = data,
        sp = sp,
        k = k,
        par0 = start$par0,
        basis = basis,
        alpha_index = start$alpha_index
    )

    fit <- maybe_reparameterise_and_refit(
        fit = fit,
        data = data,
        sp = sp,
        basis = basis,
        alpha_tol = alpha_tol,
        auto_alpha = auto_alpha
    )

    fit
}

fit_given_par0 <- function(data, sp, k, par0, basis, alpha_index = 1) {
     opt <- stats::optim(par0, loglikelihood_pen, loglikelihood_pen_grad,
                         X = basis$X, y = data$y, c = data$c - 1,
                         sp = sp, S = basis$S, K = k, alpha_index = alpha_index,
                         method = "BFGS", control = list(fnscale = -1, maxit = 10000))
    if(opt$convergence != 0)
        warning("optim has not converged")
     fit <- find_fit_info(opt, k, basis, sp, data, alpha_index)
    
     fit
}

fit_given_par0_nlm <- function(data, sp, k, par0, basis, alpha_index = 1) {
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

find_start_beta_given_fit_km1 <- function(fit_km1, k, nbasis,
                                          fit_k_other_sp = NULL,
                                          alpha_index = 1) {
    if(!is.null(fit_k_other_sp)) {
        return(list(
            beta0 = fit_k_other_sp$beta0,
            beta = fit_k_other_sp$beta,
            lsigma = fit_k_other_sp$lsigma
        ))
    }

    beta0 <- fit_km1$beta0
    lsigma <- fit_km1$lsigma

    if(k == 1) {
        alpha_start <- rep(0.01, nbasis)
    } else {
        ## Express previous beta in the target alpha chart.
        alpha_km1 <- find_alpha_from_beta(
            fit_km1$beta,
            nbasis = nbasis,
            k = k - 1,
            alpha_index = alpha_index
        )

        alpha_k0 <- rep(0.01, nbasis - k + 1)

        alpha_start <- c(alpha_km1, alpha_k0)
    }

    beta_start <- find_beta(
        alpha_start,
        nbasis = nbasis,
        k = k,
        alpha_index = alpha_index
    )

    list(
        beta0 = beta0,
        beta = beta_start,
        lsigma = lsigma
    )
}

fit_given_fit_km1 <- function(data, sp, k, fit_km1, basis,
                              fit_k_other_sp = NULL,
                              alpha_index = 1,
                              auto_alpha = TRUE,
                              alpha_tol = 1e-2) {
    start <- find_start_beta_given_fit_km1(
        fit_km1 = fit_km1,
        k = k,
        nbasis = basis$nbasis,
        fit_k_other_sp = fit_k_other_sp,
        alpha_index = alpha_index
    )

    fit_given_start_beta(
        data = data,
        sp = sp,
        k = k,
        beta0 = start$beta0,
        beta = start$beta,
        lsigma = start$lsigma,
        basis = basis,
        alpha_index = alpha_index,
        auto_alpha = auto_alpha,
        alpha_tol = alpha_tol
    )
}

find_FVE <- function(mod) {
    cumsum(mod$lambda) / sum(mod$lambda)
}


is_k_larger_than_required <- function(mod, k_tol) {
    FVE <-  find_FVE(mod)
    FVE[length(FVE) - 1] > 1 - k_tol
}


fits_given_sp <- function(sp, kmax, data, basis, k_tol,
                          fits_other_sp = NULL,
                          alpha_index = 1,
                          auto_alpha = TRUE,
                          alpha_tol = 1e-2) {
    fits <- list(fit_0(data, sp, basis, alpha_index))
    for(k in 1:kmax) {
        if(length(fits_other_sp) > k)
            fit_k_other_sp <- fits_other_sp[[k+1]]
        else
            fit_k_other_sp <- NULL
        
        fits[[k+1]] <- fit_given_fit_km1(
            data = data,
            sp = sp,
            k = k,
            fit_km1 = fits[[k]],
            basis = basis,
            fit_k_other_sp = fit_k_other_sp,
            alpha_index = alpha_index,
            auto_alpha = auto_alpha,
            alpha_tol = alpha_tol
        )
        
        if(k > 1)
            if(is_k_larger_than_required(fits[[k+1]], k_tol))
                break
    }
    fits
}
