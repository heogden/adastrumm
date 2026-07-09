make_J_zeta <- function(fit, basis) {
    nbasis <- basis$nbasis
    k <- fit$k

    if(k == 0) {
        return(diag(nbasis))
    }

    J_alpha <- jac_beta_alpha(
        alpha = fit$alpha,
        K = k,
        n_B = nbasis,
        psi_index = fit$psi_index
    )

    as.matrix(Matrix::bdiag(diag(nbasis), J_alpha))
}

make_P_beta <- function(fit, basis) {
    nbasis <- basis$nbasis
    k <- fit$k

    ## Penalty on beta0, beta1, ..., betaK
    S_blocks <- replicate(k + 1, basis$S, simplify = FALSE)

    P0 <- as.matrix(Matrix::bdiag(S_blocks))

    (fit$sp / fit$sigma^2) * P0
}

approx_log_ml <- function(fit, hessian, basis) {
    k <- fit$k
    
    J_zeta <- make_J_zeta(fit, basis)
    G_zeta <- t(J_zeta) %*% J_zeta
    G_zeta <- (G_zeta + t(G_zeta)) / 2

    P_beta <- make_P_beta(fit, basis)
    P_zeta <- t(J_zeta) %*% P_beta %*% J_zeta
    P_zeta <- (P_zeta + t(P_zeta)) / 2

    H_zeta <- -hessian[seq_len(nrow(hessian) - 1),
                       seq_len(ncol(hessian) - 1),
                       drop = FALSE]
    H_zeta <- (H_zeta + t(H_zeta)) / 2

    P_rel <- log_det_relative(P_zeta, G_zeta, pseudo = TRUE)
    H_rel <- log_det_relative(H_zeta, G_zeta, pseudo = FALSE)

    fit$l_pen +
        0.5 * P_rel$log_det -
        0.5 * H_rel$log_det +
        0.5 * P_rel$null_dim * log(2 * pi) -
        lfactorial(k) -
        k * log(2)
}
