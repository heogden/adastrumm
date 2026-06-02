find_alpha_components <- function(nbasis, k) {
    n_each <- nbasis - 0:(k-1)
    rep(1:k, times = n_each)
}


split_alpha <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    split(alpha, component)
}

hh_sign <- function(alpha, psi_index) {
    if (alpha[psi_index] >= 0) 1 else -1
}

find_Hstar_mat <- function(alpha, psi_index) {
    alpha_norm <- sqrt(sum(alpha^2))

    u <- alpha
    u[psi_index] <- u[psi_index] + hh_sign(alpha, psi_index) * alpha_norm

    t <- sum(u^2)
    
    gamma <- 2 / t

    H <- diag(nrow = length(u)) - gamma * outer(u, u)
    H[ , -psi_index, drop = FALSE]
}


find_T_list <- function(alpha, nbasis, k, psi_index = 1) {
    alpha_list <- split_alpha(alpha, nbasis, k)

    T_list <- list()
    # no transformation for k = 1: identity matrix
    T_list[[1]] <- diag(nrow = length(alpha_list[[1]]))
    
    if(k > 1) {
        for(j in 2:k) {
            T_list[[j]] <- T_list[[j-1]] %*% find_Hstar_mat(alpha_list[[j-1]], psi_index)
        }
    }
    
    T_list
}


## TODO:: could get transform from C++ code instead
find_Hstar <- function(alpha, psi_index) {
    alpha_norm <- sqrt(sum(alpha^2))

    u <- alpha
    u[psi_index] <- u[psi_index] + hh_sign(alpha, psi_index) * alpha_norm
    t <- sum(u^2)
    gamma <- 2 / t
    
    list(u = u, gamma = gamma, psi_index = psi_index)
}


find_Hstar_x <- function(Hstar, x) {
    psi_index <- Hstar$psi_index

    x_ext <- numeric(length(x) + 1)
    x_ext[-psi_index] <- x

    a <- sum(Hstar$u * x_ext)

    x_ext - a * Hstar$gamma * Hstar$u
}


find_beta_i <- function(alpha_i, Hstar_list) {
    i <- length(Hstar_list) + 1
    s_j <- alpha_i
    
    for(j in (i-1):1) {
        s_j <- find_Hstar_x(Hstar_list[[j]], s_j)
    }
    
    s_j
}

find_beta <- function(alpha, nbasis, k, psi_index = 1) {
    component <- find_alpha_components(nbasis, k)
    Hstar_list <- list()
    beta <- matrix(nrow = nbasis, ncol = k)
        
    for(i in 1:k) {
        alpha_i <- alpha[component == i]

        if(psi_index < 1 || psi_index > length(alpha_i)) {
            stop("psi_index must be between 1 and the length of each alpha block")
        }
        
        if(i == 1)
            beta_i <- alpha_i
        else
            beta_i <- find_beta_i(alpha_i, Hstar_list)

        beta[,i] <- beta_i

        if(i < k) {
            Hstar_list[[i]] <- find_Hstar(alpha_i, psi_index)
        }
    }
    beta
}


find_alpha_from_beta <- function(beta, nbasis, k, psi_index = 1) {
    if(is.vector(beta)) {
        beta <- matrix(beta, nrow = nbasis, ncol = k)
    }

    stopifnot(nrow(beta) == nbasis, ncol(beta) == k)

    alpha_list <- vector("list", k)
    Hstar_list <- vector("list", max(k - 1, 0))

    for(i in seq_len(k)) {
        x <- beta[, i]

        if(i > 1) {
            for(j in seq_len(i - 1)) {
                x <- as.vector(t(Hstar_list[[j]]) %*% x)
            }
        }

        alpha_list[[i]] <- x

        if(i < k) {
            Hstar_list[[i]] <- find_Hstar_mat(x, psi_index)
        }
    }

    unlist(alpha_list, use.names = FALSE)
}
 
choose_psi_index_from_beta <- function(beta, nbasis, k) {
    if(k <= 1) {
        return(1)
    }

    candidates <- seq_len(nbasis - k + 1)

    scores <- sapply(candidates, function(psi_index) {
        alpha <- find_alpha_from_beta(beta, nbasis = nbasis,
                                      k = k, psi_index = psi_index)

        psi_diagnostic(alpha, nbasis = nbasis, k = k,
                       psi_index = psi_index)$min_ratio
    })

    candidates[which.max(scores)]
}

psi_from_theta_parameterisation <- function(beta0, beta, lsigma,
                                           nbasis, k, psi_index) {
    alpha <- find_alpha_from_beta(beta,
                                  nbasis = nbasis,
                                  k = k,
                                  psi_index = psi_index)

    c(beta0, alpha, lsigma)
}

start_from_fit_beta <- function(fit) {
    list(beta0 = fit$beta0,
         beta = fit$beta,
         lsigma = fit$lsigma)
}

maybe_switch_psi_index_start <- function(beta0, beta, lsigma,
                                         basis, k,
                                         psi_index = 1,
                                         psi_tol = 1e-5,
                                         auto_psi = TRUE) {
    nbasis <- basis$nbasis
    
    if(psi_index > basis$nbasis - k + 1) {
        psi_index <- 1
    }
    
    alpha <- find_alpha_from_beta(beta,
                                  nbasis = nbasis,
                                  k = k,
                                  psi_index = psi_index)

    if(auto_psi && k > 1) {
        diag <- psi_diagnostic(alpha,
                               nbasis = nbasis,
                               k = k,
                               psi_index = psi_index)

        if(diag$min_ratio < psi_tol) {
            psi_index <- choose_psi_index_from_beta(beta, nbasis, k)

            alpha <- find_alpha_from_beta(beta,
                                          nbasis = nbasis,
                                          k = k,
                                          psi_index = psi_index)
        }
    }

    list(psi0 = c(beta0, alpha, lsigma),
         psi_index = psi_index)
}

order_beta_by_lambda <- function(beta) {
    if(is.null(beta)) {
        return(list(
            beta = NULL,
            lambda = NULL,
            order = integer(0)
        ))
    }

    lambda <- colSums(beta^2)

    if(ncol(beta) <= 1) {
        return(list(
            beta = beta,
            lambda = lambda,
            order = seq_len(ncol(beta))
        ))
    }

    ord <- order(lambda, decreasing = TRUE)

    list(
        beta = beta[, ord, drop = FALSE],
        lambda = lambda[ord],
        order = ord
    )
}



psi_diagnostic <- function(alpha, nbasis, k, psi_index = 1,
                           eps = .Machine$double.eps) {
    if(k <= 1) {
        return(list(
            min_ratio = Inf,
            ratios = numeric(0),
            alpha_norms = numeric(0),
            alpha_selected = numeric(0),
            psi_index = psi_index
        ))
    }
    alpha_list <- split_alpha(alpha, nbasis, k)
    
    

    alpha_used <- alpha_list[seq_len(k - 1)]

    alpha_norms <- vapply(alpha_used, function(a) sqrt(sum(a^2)), numeric(1))
    alpha_selected <- vapply(alpha_used, function(a) a[psi_index], numeric(1))

    ratios <- abs(alpha_selected) / pmax(alpha_norms, eps)

    list(
        min_ratio = min(ratios),
        ratios = ratios,
        alpha_norms = alpha_norms,
        alpha_selected = alpha_selected,
        psi_index = psi_index
    )
}

psi_ci_diagnostic_from_cov <- function(psi, V, nbasis, k, psi_index,
                                       eps = .Machine$double.eps) {
    if(k <= 1) {
        return(list(
            min_z_to_boundary = Inf,
            z_to_boundary = numeric(0),
            selected_alpha = numeric(0),
            selected_se = numeric(0),
            psi_index = psi_index
        ))
    }

    psi_split <- split_psi(psi, nbasis)
    alpha <- psi_split$alpha

    V <- 0.5 * (V + t(V))
    se <- sqrt(pmax(diag(V), eps))

    alpha_start <- nbasis + 1

    alpha_components <- split(
        seq_along(alpha),
        find_alpha_components(nbasis, k)
    )

    selected_alpha_positions <- vapply(seq_len(k - 1), function(j) {
        alpha_start + alpha_components[[j]][psi_index] - 1
    }, numeric(1))

    selected_alpha <- psi[selected_alpha_positions]
    selected_se <- se[selected_alpha_positions]

    z_to_boundary <- abs(selected_alpha) / pmax(selected_se, eps)

    list(
        min_z_to_boundary = min(z_to_boundary),
        z_to_boundary = z_to_boundary,
        selected_alpha = selected_alpha,
        selected_se = selected_se,
        psi_index = psi_index
    )
}

psi_ci_diagnostic <- function(mod, eps = .Machine$double.eps) {
    psi_ci_diagnostic_from_cov(
        psi = mod$psi,
        V = mod$var_par,
        nbasis = mod$basis$nbasis,
        k = mod$k,
        psi_index = mod$psi_index,
        eps = eps
    )
}

approx_cov_for_psi_index <- function(mod, psi_index_new) {
    nbasis <- mod$basis$nbasis
    k <- mod$k

    if(k <= 1) {
        return(list(
            psi = mod$psi,
            V = mod$var_par
        ))
    }

    V_current <- 0.5 * (mod$var_par + t(mod$var_par))

    G_current <- jac_beta_alpha(
        alpha = mod$alpha,
        K = k,
        n_B = nbasis,
        psi_index = mod$psi_index
    )

    alpha_new <- find_alpha_from_beta(
        beta = mod$beta,
        nbasis = nbasis,
        k = k,
        psi_index = psi_index_new
    )

    G_new <- jac_beta_alpha(
        alpha = alpha_new,
        K = k,
        n_B = nbasis,
        psi_index = psi_index_new
    )

    ## Approximate map d alpha_current -> d alpha_new.
    A_alpha <- qr.solve(G_new, G_current)

    A_full <- diag(length(mod$psi))

    alpha_ind <- (nbasis + 1):(length(mod$psi) - 1)
    A_full[alpha_ind, alpha_ind] <- A_alpha

    V_new <- A_full %*% V_current %*% t(A_full)
    V_new <- 0.5 * (V_new + t(V_new))

    psi_new <- c(mod$beta0, alpha_new, mod$lsigma)

    list(
        psi = psi_new,
        V = V_new,
        A_alpha = A_alpha
    )
}

approx_psi_ci_scores <- function(mod) {
    nbasis <- mod$basis$nbasis
    k <- mod$k

    candidates <- seq_len(nbasis - k + 1)

    scores <- vapply(candidates, function(psi_index_new) {
        tryCatch({
            approx <- approx_cov_for_psi_index(mod, psi_index_new)

            diag <- psi_ci_diagnostic_from_cov(
                psi = approx$psi,
                V = approx$V,
                nbasis = nbasis,
                k = k,
                psi_index = psi_index_new
            )

            diag$min_z_to_boundary
        }, error = function(e) {
            -Inf
        })
    }, numeric(1))

    list(
        candidates = candidates,
        scores = scores
    )
}
