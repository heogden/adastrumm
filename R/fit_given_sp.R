add_hessian_and_log_ml <- function(fit, basis, data) {
    fit$hessian <- loglikelihood_pen_hess(fit$psi,
                                          X = basis$X, y = data$y, c = data$c - 1,
                                          sp = fit$sp, S = basis$S, K = fit$k,
                                          psi_index = fit$psi_index)
    var_par <- tryCatch(solve(-fit$hessian),
                        error = function(cond) {
                            pracma::pinv(-fit$hessian)
                        })

    fit$var_par <- 0.5 * (var_par + t(var_par))
    
    fit$log_ml <- approx_log_ml(fit, fit$hessian, basis)
    fit
}



split_psi <- function(psi, nbasis) {
    if(length(psi) < nbasis + 1) {
        stop("psi is too short for the specified nbasis")
    }
    
    components <- rep("alpha", length(psi))
    components[1:nbasis] <- "beta0"
    components[length(psi)] <- "lsigma"
    split(psi, components)
}

find_par_cluster <- function(beta0, beta, u_hat) {
    beta0 + tcrossprod(beta, u_hat)
}


find_fit_info <- function(opt, k, basis, sp, data, psi_index) {
    psi <- opt$par
    l_pen <- opt$value
    psi_split <- split_psi(psi, basis$nbasis)


    f0_x <- basis$X %*% psi_split$beta0
    f0 <- find_spline_fun(psi_split$beta0, basis)
    
    if(k > 0) {
        beta <- find_beta(psi_split$alpha, basis$nbasis, k, psi_index)
        lambda <- colSums(beta^2)
        f_x <- basis$X %*% beta
        u_hat_full <- find_u_hat(exp(psi_split$lsigma), data, f0_x, f_x, var = TRUE)
        u_hat <- u_hat_full$u_hat
        var_u_hat <- u_hat_full$var_u_hat
        f <- find_spline_fun(beta, basis)
        par_cluster <- find_par_cluster(psi_split$beta0, beta, u_hat)
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
         psi = psi,
         psi_index = psi_index,
         l_pen = l_pen,
         opt = opt,
         beta0 = psi_split$beta0,
         alpha = psi_split$alpha,
         lsigma = psi_split$lsigma,
         beta = beta,
         sigma = exp(psi_split$lsigma),
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


reparameterise_fit_without_optim <- function(fit, basis, data, sp, psi_index) {
    psi_new <- psi_from_theta_parameterisation(
        beta0 = fit$beta0,
        beta = fit$beta,
        lsigma = fit$lsigma,
        nbasis = basis$nbasis,
        k = fit$k,
        psi_index = psi_index
    )

    opt_new <- list(
        par = psi_new,
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
        psi_index = psi_index
    )

    add_hessian_and_log_ml(fit_new, basis, data)
}


maybe_switch_psi_for_ci_approx <- function(fit, data, sp, basis,
                                           psi_ci_tol = 2) {
    ## This assumes the Hessian/var_par for the current fit is usable.

    current_diag <- psi_ci_diagnostic(fit)
    fit$psi_ci_diagnostic <- current_diag

    ## If the current index is OK, do not search/switch.
    if(current_diag$min_z_to_boundary >= psi_ci_tol) {
        return(fit)
    }

    ## Current index is bad: cheaply approximate candidate scores.
    approx <- approx_psi_ci_scores(fit)

    fit$psi_ci_scores_approx <- approx$scores
    fit$psi_index_candidates_approx <- approx$candidates

    valid <- is.finite(approx$scores)

    if(!any(valid)) {
        return(fit)
    }

    best_candidate <- approx$candidates[valid][which.max(approx$scores[valid])]

    ## If the approximate best is the current index, keep current fit.
    if(best_candidate == fit$psi_index) {
        return(fit)
    }

    ## Now do the expensive exact check for the best approximate candidate only.
    candidate_fit <- reparameterise_fit_without_optim(
        fit = fit,
        basis = basis,
        data = data,
        sp = sp,
        psi_index = best_candidate
    )

    candidate_diag <- psi_ci_diagnostic(candidate_fit)
    candidate_fit$psi_ci_diagnostic <- candidate_diag
    candidate_fit$psi_ci_scores_approx <- approx$scores
    candidate_fit$psi_index_candidates_approx <- approx$candidates

    ## Switch only if the exact diagnostic improves.
    if(candidate_diag$min_z_to_boundary >
       current_diag$min_z_to_boundary) {
        return(candidate_fit)
    }

    fit
}

maybe_reparameterise_after_hessian <- function(fit, data, sp, basis,
                                               psi_tol = 1e-5,
                                               psi_ci_tol = 2,
                                               auto_psi = TRUE) {
    if(!auto_psi || fit$k <= 1) {
        return(fit)
    }

    bad_hessian <- !is_neg_def(fit$hessian) ||
        matrixcalc::is.singular.matrix(fit$hessian)

    point_diag <- psi_diagnostic(
        fit$alpha,
        nbasis = basis$nbasis,
        k = fit$k,
        psi_index = fit$psi_index
    )

    fit$psi_diagnostic <- point_diag

    bad_point <- point_diag$min_ratio < psi_tol

    if(bad_hessian || bad_point) {
        psi_index_new <- choose_psi_index_from_beta(
            fit$beta,
            nbasis = basis$nbasis,
            k = fit$k
        )

        if(psi_index_new == fit$psi_index) {
            if(bad_hessian) {
                return(fit)
            }

            return(
                maybe_switch_psi_for_ci_approx(
                    fit = fit,
                    data = data,
                    sp = sp,
                    basis = basis,
                    psi_ci_tol = psi_ci_tol
                )
            )
        }

        psi0_new <- psi_from_theta_parameterisation(
            beta0 = fit$beta0,
            beta = fit$beta,
            lsigma = fit$lsigma,
            nbasis = basis$nbasis,
            k = fit$k,
            psi_index = psi_index_new
        )

        fit_new <- fit_given_psi0(
            data = data,
            sp = sp,
            k = fit$k,
            psi0 = psi0_new,
            basis = basis,
            psi_index = psi_index_new
        )

        fit_new <- add_hessian_and_log_ml(fit_new, basis, data)

        if(!is_neg_def(fit_new$hessian) ||
           matrixcalc::is.singular.matrix(fit_new$hessian)) {
            return(fit_new)
        }
        
        return(
            maybe_switch_psi_for_ci_approx(
                fit = fit_new,
                data = data,
                sp = sp,
                basis = basis,
                psi_ci_tol = psi_ci_tol
            )
        )
    }

    maybe_switch_psi_for_ci_approx(
        fit = fit,
        data = data,
        sp = sp,
        basis = basis,
        psi_ci_tol = psi_ci_tol
    )
}

maybe_reparameterise_and_refit <- function(fit, data, sp, basis,
                                           psi_tol = 1e-5,
                                           auto_psi = TRUE) {
    if(!auto_psi || fit$k <= 1) {
        return(fit)
    }

    diag <- psi_diagnostic(
        fit$alpha,
        nbasis = basis$nbasis,
        k = fit$k,
        psi_index = fit$psi_index
    )

    if(diag$min_ratio >= psi_tol) {
        return(fit)
    }

    psi_index_new <- choose_psi_index_from_beta(
        fit$beta,
        nbasis = basis$nbasis,
        k = fit$k
    )

    if(psi_index_new == fit$psi_index) {
        return(fit)
    }

    psi0_new <- psi_from_theta_parameterisation(
        beta0 = fit$beta0,
        beta = fit$beta,
        lsigma = fit$lsigma,
        nbasis = basis$nbasis,
        k = fit$k,
        psi_index = psi_index_new
    )

    fit_given_psi0(
        data = data,
        sp = sp,
        k = fit$k,
        psi0 = psi0_new,
        basis = basis,
        psi_index = psi_index_new
    )
}

fit_given_start_beta <- function(data, sp, k,
                                 beta0, beta, lsigma,
                                 basis,
                                 psi_index = 1,
                                 auto_psi = TRUE,
                                 psi_tol = 1e-5) {
    start <- maybe_switch_psi_index_start(
        beta0 = beta0,
        beta = beta,
        lsigma = lsigma,
        basis = basis,
        k = k,
        psi_index = psi_index,
        psi_tol = psi_tol,
        auto_psi = auto_psi
    )

    fit <- fit_given_psi0(
        data = data,
        sp = sp,
        k = k,
        psi0 = start$psi0,
        basis = basis,
        psi_index = start$psi_index
    )

    fit <- maybe_reparameterise_and_refit(
        fit = fit,
        data = data,
        sp = sp,
        basis = basis,
        psi_tol = psi_tol,
        auto_psi = auto_psi
    )

    fit
}

fit_given_psi0 <- function(data, sp, k, psi0, basis, psi_index = 1) {
     opt <- stats::optim(psi0, loglikelihood_pen, loglikelihood_pen_grad,
                         X = basis$X, y = data$y, c = data$c - 1,
                         sp = sp, S = basis$S, K = k, psi_index = psi_index,
                         method = "BFGS", control = list(fnscale = -1, maxit = 10000))
    if(opt$convergence != 0)
        warning("optim has not converged")
     fit <- find_fit_info(opt, k, basis, sp, data, psi_index)

     fit <- order_fit_components_by_lambda(fit, basis, data, sp)
     
     fit
}

fit_given_psi0_nlm <- function(data, sp, k, psi0, basis, psi_index = 1) {
    ml <- function(psi) {
        result <- -loglikelihood_pen(psi, X = basis$X, y = data$y,
                                     c = data$c - 1,
                                     sp = sp, S = basis$S,
                                     K = k, psi_index = psi_index)
        gr <- -loglikelihood_pen_grad(psi,
                                      X = basis$X, y = data$y, c = data$c - 1,
                                      sp = sp, S = basis$S, K = k,
                                      psi_index = psi_index)
        attr(result, "gradient") <- gr
        result
    }
    opt_nlm <- stats::nlm(ml, psi0, iterlim = 100000, check.analyticals = FALSE)

    # convert to format from optim
    list(par = opt_nlm$estimate,
         value = -opt_nlm$minimum,
         convergence = ifelse(opt_nlm$code < 2, 0, opt_nlm$code))
                
}


# Fit the mean-only model
fit_0 <- function(data, sp, basis, psi_index = 1) {
    X_0 <- basis$X
    
    Xt_y <- crossprod(X_0, data$y)
    XtX <- crossprod(X_0, X_0)

    
    beta_0 <- as.numeric(solve(XtX + sp * basis$S, Xt_y))
    y_hat_0 <- X_0 %*% beta_0
    resid <- data$y - y_hat_0
    sigma <- stats::sd(resid)

    psi0 <- c(beta_0, log(sigma))
    
    fit_given_psi0(data, sp, 0, psi0, basis, psi_index)
}

find_start_beta_given_fit_km1 <- function(fit_km1, k, nbasis,
                                          fit_k_other_sp = NULL,
                                          psi_index = 1) {
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
        alpha_km1 <- find_alpha_from_beta(
            fit_km1$beta,
            nbasis = nbasis,
            k = k - 1,
            psi_index = psi_index
        )

        alpha_k0 <- rep(0.01, nbasis - k + 1)

        alpha_start <- c(alpha_km1, alpha_k0)
    }

    beta_start <- find_beta(
        alpha_start,
        nbasis = nbasis,
        k = k,
        psi_index = psi_index
    )

    list(
        beta0 = beta0,
        beta = beta_start,
        lsigma = lsigma
    )
}

fit_given_fit_km1 <- function(data, sp, k, fit_km1, basis,
                              fit_k_other_sp = NULL,
                              psi_index = 1,
                              auto_psi = TRUE,
                              psi_tol = 1e-5) {
    start <- find_start_beta_given_fit_km1(
        fit_km1 = fit_km1,
        k = k,
        nbasis = basis$nbasis,
        fit_k_other_sp = fit_k_other_sp,
        psi_index = psi_index
    )

    fit_given_start_beta(
        data = data,
        sp = sp,
        k = k,
        beta0 = start$beta0,
        beta = start$beta,
        lsigma = start$lsigma,
        basis = basis,
        psi_index = psi_index,
        auto_psi = auto_psi,
        psi_tol = psi_tol
    )
}

order_fit_components_by_lambda <- function(fit, basis, data, sp) {
    if(fit$k <= 1 || is.null(fit$beta)) {
        return(fit)
    }

    ordered <- order_beta_by_lambda(fit$beta)

    if(identical(ordered$order, seq_len(fit$k))) {
        return(fit)
    }

    psi_new <- psi_from_theta_parameterisation(
        beta0 = fit$beta0,
        beta = ordered$beta,
        lsigma = fit$lsigma,
        nbasis = basis$nbasis,
        k = fit$k,
        psi_index = fit$psi_index
    )

    opt_new <- fit$opt
    opt_new$par <- psi_new
    opt_new$value <- fit$l_pen

    if(!is.null(opt_new$message)) {
        opt_new$message <- paste(
            opt_new$message,
            "Components reordered by decreasing lambda",
            sep = "; "
        )
    } else {
        opt_new$message <- "Components reordered by decreasing lambda"
    }

    find_fit_info(
        opt = opt_new,
        k = fit$k,
        basis = basis,
        sp = sp,
        data = data,
        psi_index = fit$psi_index
    )
}

is_ordered_lambda <- function(lambda, tol = 1e-10) {
    if(length(lambda) <= 1) {
        return(TRUE)
    }

    all(diff(lambda) <= tol)
}

find_FVE <- function(mod) {
    if(!is_ordered_lambda(mod$lambda)) {
        stop("lambda is not ordered decreasingly; reorder components before computing FVE")
    }

    cumsum(mod$lambda) / sum(mod$lambda)
}


is_k_larger_than_required <- function(mod, k_tol) {
    FVE <-  find_FVE(mod)
    FVE[length(FVE) - 1] > 1 - k_tol
}


fits_given_sp <- function(sp, kmax, data, basis, k_tol,
                          fits_other_sp = NULL,
                          psi_index = 1,
                          auto_psi = TRUE,
                          psi_tol = 1e-5) {
    fits <- list(fit_0(data, sp, basis, psi_index))
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
            psi_index = psi_index,
            auto_psi = auto_psi,
            psi_tol = psi_tol
        )
        
        if(k > 1)
            if(is_k_larger_than_required(fits[[k+1]], k_tol))
                break
    }
    fits
}
