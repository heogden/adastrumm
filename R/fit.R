normalise_data <- function(data, norm = NULL) {
    if(length(norm) == 0) {
        m_y <- mean(data$y)
        s_y <- stats::sd(data$y)
        
        m_x <- mean(data$x)
        s_x <- stats::sd(data$x)

        norm <- list(m_y = m_y, s_y = s_y, m_x = m_x, s_x = s_x)
    }
    

    data_norm <- data.frame(y =  (data$y - norm$m_y) / norm$s_y,
                            x =  (data$x - norm$m_x) / norm$s_x,
                            c = data$c)
    
    list(norm = norm, data = data_norm)
}



#' Fit an Adaptively-Structure Mixed Model (AdaStruMM)
#'
#' @param data A data frame, with columns c (identifying the
#'     individual subjects), x (the time) and y (the response).
#' @param nbasis The number of spline basis functions.
#' @param kmax The maximum number of functional principal components
#'     to allow.
#' @param k_tol. The tolerance to use in selecting k, to explain at
#'     least 1-k_tol of the variation in trajectories.
#' @param lsp_poss The grid of possible values to consider for
#'     log(gamma), the log of the smoothing parameter.
#' @param trace If TRUE, print out extra information.
#' @param psi_index Initial index for the choice of psi
#'     parameterisation. Defaults to 1.
#' @param auto_psi Logical; if TRUE, automatically switches psi-
#'     parameterisation when diagnostics indicate switching is needed.
#' @param psi_tol Threshold for the diagnostic used to trigger
#'     automatic psi-parameterisation switching.
#' @param psi_ci_tol Threshold for the confidence interval diagnostic
#'     used to trigger automatic psi-parameterisation switching.
#' @param normalise Logical; if TRUE, automatically normalise the data
#'     (response and covariate) to mean 0, variance 1.
#' @return The fitted model.
#' @examples
#' data_full <- simulate_1dv(1, -0.5, 0.1, 0.5, 0.1, 20, 10)
#' data <- data_full$data
#' mod <- fit_adastrumm(data)
#' @export
fit_adastrumm <- function(data, nbasis = 10, kmax = 10, k_tol = 1e-4,
                          lsp_poss = -5:15, trace = FALSE,
                          psi_index = 1,
                          auto_psi = TRUE,
                          psi_tol = 1e-5,
                          psi_ci_tol = 2,
                          normalise = TRUE) {
    if(any(is.na(data)))
        stop("There are missing values in the data, which adastrumm cannot handle")

    if(psi_index != as.integer(psi_index) || psi_index < 1) {
        stop("psi_index must be a positive integer")
}

    if(psi_index > nbasis) {
        stop("psi_index must be no larger than nbasis")
    }

    
    if(normalise) {
        norm <- NULL
    } else {
        norm <- list(m_y = 0, s_y = 1, m_x = 0, s_x = 1)
    }
    
    data_norm_full <- normalise_data(data, norm)
    data <- data_norm_full$data
    norm <- data_norm_full$norm
    
    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    
    sp_poss <- exp(lsp_poss)

    fits_list <- list()
    fit_sp_poss <- list()
    
    for(i in seq_along(sp_poss)) {
        sp <- sp_poss[i]
        if(trace)
            cat("sp = ", sp, "\n")
        if(i == 1)
            fits_other_sp <- NULL
        else
            fits_other_sp <- fits_list[[i-1]]
        
        fits <- fits_given_sp(
            sp = sp,
            kmax = kmax,
            data = data,
            basis = basis,
            k_tol = k_tol,
            fits_other_sp = fits_other_sp,
            psi_index = psi_index,
            auto_psi = auto_psi,
            psi_tol = psi_tol
        )

        
        if(is_k_larger_than_required(fits[[length(fits)]], k_tol))
            fit <- fits[[length(fits) - 1]]
        else
            fit <- fits[[length(fits)]]
                                        
        if(trace) {
            cat("k = ", fit$k, "\n")
            cat("lambda = ", fit$lambda, "\n")
            cat("FVE = ", find_FVE(fit), "\n")
        }
        
        fit <- add_hessian_and_log_ml(fit, basis, data)
        
        fit <- maybe_reparameterise_after_hessian(
            fit = fit,
            data = data,
            sp = sp,
            basis = basis,
            psi_tol = psi_tol,
            psi_ci_tol = psi_ci_tol,
            auto_psi = auto_psi
        )

        fits[[fit$k + 1]] <- fit
        fits_list[[i]] <- fits
        
        if(!is_neg_def(fit$hessian) | matrixcalc::is.singular.matrix(fit$hessian)) {
            message("fit from optim gave non-negative definite or singular Hessian. Trying with nlm. ")
            opt <- NULL
            try(opt <- fit_given_psi0_nlm(data, sp, fit$k, fit$psi, basis, fit$psi_index))
            if(length(opt) > 0) {    
                fit <- find_fit_info(opt, fit$k, basis, sp, data, fit$psi_index)
                fit <- order_fit_components_by_lambda(fit, basis, data, sp)
                fit <- add_hessian_and_log_ml(fit, basis, data)
                fit <- maybe_reparameterise_after_hessian(
                    fit = fit,
                    data = data,
                    sp = sp,
                    basis = basis,
                    psi_tol = psi_tol,
                    psi_ci_tol = psi_ci_tol,
                    auto_psi = auto_psi
                )
                fits[[fit$k + 1]] <- fit
                fits_list[[i]] <- fits
            }
            
            if(!is_neg_def(fit$hessian) | matrixcalc::is.singular.matrix(fit$hessian) | length(opt) == 0) {
                if(fit$k > 1) {
                    message("fit from nlm failed or gave a non-negative definite or singular Hessian. Reducing k. ")
                    fit <- fits[[fit$k]]
                    fit <- add_hessian_and_log_ml(fit, basis, data)
                    fit <- maybe_reparameterise_after_hessian(
                        fit = fit,
                        data = data,
                        sp = sp,
                        basis = basis,
                        psi_tol = psi_tol,
                        psi_ci_tol = psi_ci_tol,
                        auto_psi = auto_psi
                    )
                    fits[[fit$k + 1]] <- fit
                    fits_list[[i]] <- fits
                }

            }
            
        }

 
        
        fit$norm <- norm
        
        if(trace)
            cat("log_ml = ", fit$log_ml, "\n")
        fit_sp_poss[[i]] <- fit

        if(i > 2)
            if(fit_sp_poss[[i]]$log_ml < fit_sp_poss[[i-1]]$log_ml)
                break         
    }

    log_ml_sp_poss <- sapply(fit_sp_poss, "[[", "log_ml")
    i_opt <- which.max(log_ml_sp_poss)
    fit_sp_poss[[i_opt]]
}

