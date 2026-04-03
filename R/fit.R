#' Fit an Adaptively-Structure Mixed Model (AdaStruMM)
#'
#' @param data A data frame, with columns c (identifying the
#'     individual subjects), x (the time) and y (the response).
#' @param nbasis The number of spline basis functions.
#' @param kmax The maximum number of functional principal components
#'     to allow.
#' @param k_tol. The tolerance to use in selecting k, don't use
#'     eigenvalues smaller than k_tol
#' @param lsp_poss The grid of possible values to consider for
#'     log(gamma), the log of the smoothing parameter.
#' @param trace If TRUE, print out extra information.
#' @return The fitted model.
#' @examples
#' data_full <- simulate_1dv(1, -0.5, 0.1, 0.5, 0.1, 20, 10)
#' data <- data_full$data
#' mod <- fit_adastrumm(data)
#' @export
fit_adastrumm <- function(data, nbasis = 10, kmax = 10, k_tol = 1e-6,
                          lsp_poss = -5:15, trace = FALSE) {
    if(any(is.na(data)))
        stop("There are missing values in the data, which adastrumm cannot handle")

    m_y <- mean(data$y)
    s_y <- stats::sd(data$y)

    m_x <- mean(data$x)
    s_x <- stats::sd(data$x)

    norm <- list(m_y = m_y, s_y = s_y, m_x = m_x, s_x = s_x)

    data$y <- (data$y - m_y) / s_y
    data$x <- (data$x - m_x) / s_x
    
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
        
        fits <- fits_given_sp(sp, kmax, data, basis, k_tol, fits_other_sp)
        if(is_k_larger_than_required(fits[[length(fits)]], k_tol))
            fit <- fits[[length(fits) - 1]]
        else
            fit <- fits[[length(fits)]]
               

        fits_list[[i]] <- fits
                                        
        if(trace)
            cat("k = ", fit$k, "\n")
        fit <- add_hessian_and_log_ml(fit, basis, data)

        if(trace) {
            cat("Hessian is: \n")
            print(fit$hessian)
        }
        
        
        if(!is_neg_def(fit$hessian) | matrixcalc::is.singular.matrix(fit$hessian)) {
            message("fit from optim gave non-negative definite or singular Hessian. Trying with nlm. ")
            opt <- tryCatch(
                fit_given_par0_nlm(data, sp, fit$k, fit$par, basis),
                error = function(e) {
                    message("nlm fitting failed. ")
                    NA
                }
            )
            if(!is.na(opt)) {    
                fit <- find_fit_info(opt, fit$k, basis, sp, data)
                fit <- add_hessian_and_log_ml(fit, basis, data)
            }
            
            if(!is_neg_def(fit$hessian) | matrixcalc::is.singular.matrix(fit$hessian) | is.na(opt)) {
                if(fit$k > 1) {
                    message("fit from nlm failed or gave a non-negative definite or singular Hessian. Reducing k. ")
                    fit <- fits[[fit$k]]
                    fit <- add_hessian_and_log_ml(fit, basis, data)
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

