correct_lprior_alpha <- function(fit, basis) {
    k <- fit$k
    if(k > 1) {
        T_list <- find_T_list(fit$alpha, basis$nbasis, fit$k, fit$psi_index)
    }
    
    S <- basis$S
    nbasis <- basis$nbasis
    contrib <- c()

    for(j in 0:k) {
        if(j > 1) {
            T_j <- T_list[[j]]
            S_j <- t(T_j) %*% S %*% T_j
        } else {
            S_j <- S
        }
        r_j <- min(nbasis - 2, nbasis - j + 1)
        ldg_S_j_spr <- log_det_gen(S_j, r_j) + r_j * (log(fit$sp) - 2 * fit$lsigma)
        contrib[j+1] <- -r_j * log(2 * pi) / 2 + ldg_S_j_spr / 2
    }

    sum(contrib)
}

approx_log_ml <- function(fit, hessian, basis) {
    p <- nrow(hessian)
    fit$l_pen + correct_lprior_alpha(fit, basis) + p/2 * log(2*pi) - 1/2 * log_det(hessian)
}

