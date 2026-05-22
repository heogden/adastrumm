#' Update an Adaptively-Structure Mixed Model (AdaStruMM)
#'
#' @param data A data frame, with columns c (identifying the
#'     individual subjects), x (the time) and y (the response).
#' @param prev_fit a previous model fit, fit from a subset of the data.
#' @param fix_tuning logical: should the tuning parameters K and gamma be fixed?
#' @param fix_pop_par logical: should the population-level parameters be fixed?
#' @return The fitted model.
#' @export
update_adastrumm <- function(data, prev_fit, fix_tuning = TRUE, fix_pop_par = FALSE) {
    if(!fix_tuning) {
        stop("update_adastrumm with tuning parameters re-estimated not yet available. Use fit_adastrumm instead.")
    }
    
    data_norm_full <- normalise_data(data, norm = prev_fit$norm)
    data <- data_norm_full$data
    basis <- update_basis(prev_fit$basis, data$x)
    
    if(fix_pop_par) {
        fit <- find_fit_info(prev_fit$opt, prev_fit$k, basis, prev_fit$sp, data, alpha_index = prev_fit$alpha_index)
    } else {
        fit <- fit_given_par0(data, prev_fit$sp, prev_fit$k, prev_fit$par, basis, alpha_index = prev_fit$alpha_index)
        
    }
    fit$norm <- prev_fit$norm
    fit
}
