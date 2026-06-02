drop_attributes <- function(x) {
    attributes(x) <- NULL
    x
}


log_det <- function(x) {
    drop_attributes(
        determinant(x, logarithm = TRUE)$modulus
    )
}

# the log of the generalized determinant
log_det_gen <- function(x, rank) {
    evs <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
    sum(log(evs[1:rank]))
}

is_neg_def <- function(A) {
    if(any(is.na(A)))
        return(FALSE)
   
    v <- eigen(A, only.values = TRUE)$values
    if(!is.numeric(v))
        return(FALSE)
    all(v < 0)
}

log_det_relative <- function(A, G, pseudo = FALSE, tol = NULL) {
    A <- (A + t(A)) / 2
    G <- (G + t(G)) / 2

    ## chol returns R such that G = t(R) %*% R
    R <- chol(G)
    R_inv <- backsolve(R, diag(ncol(R)))

    ## Eigenvalues of A v = lambda G v
    C <- t(R_inv) %*% A %*% R_inv
    C <- (C + t(C)) / 2

    ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values

    if(is.null(tol)) {
        tol <- max(dim(A)) * max(abs(ev), 1) * .Machine$double.eps
    }

    if(pseudo) {
        use <- ev > tol
    } else {
        if(any(ev <= tol)) {
            stop("Relative determinant requested, but matrix is not positive definite.")
        }
        use <- rep(TRUE, length(ev))
    }

    list(
        log_det = sum(log(ev[use])),
        eigenvalues = ev,
        rank = sum(use),
        null_dim = length(ev) - sum(use)
    )
}
