find_alpha_components <- function(nbasis, k) {
    n_each <- nbasis - 0:(k-1)
    rep(1:k, times = n_each)
}


split_alpha <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    split(alpha, component)
}

hh_sign <- function(alpha, alpha_index) {
    if (alpha[alpha_index] >= 0) 1 else -1
}

find_Hstar_mat <- function(alpha, alpha_index) {
    alpha_norm <- sqrt(sum(alpha^2))

    u <- alpha
    u[alpha_index] <- u[alpha_index] + hh_sign(alpha, alpha_index) * alpha_norm

    t <- sum(u^2)
    
    gamma <- 2 / t

    H <- diag(nrow = length(u)) - gamma * outer(u, u)
    H[ , -alpha_index, drop = FALSE]
}


find_T_list <- function(alpha, nbasis, k, alpha_index = 1) {
    alpha_list <- split_alpha(alpha, nbasis, k)

    T_list <- list()
    # no transformation for k = 1: identity matrix
    T_list[[1]] <- diag(nrow = length(alpha_list[[1]]))
    
    if(k > 1) {
        for(j in 2:k) {
            T_list[[j]] <- T_list[[j-1]] %*% find_Hstar_mat(alpha_list[[j-1]], alpha_index)
        }
    }
    
    T_list
}


## TODO:: could get transform from C++ code instead
find_Hstar <- function(alpha, alpha_index) {
    alpha_norm <- sqrt(sum(alpha^2))

    u <- alpha
    u[alpha_index] <- u[alpha_index] + hh_sign(alpha, alpha_index) * alpha_norm
    t <- sum(u^2)
    gamma <- 2 / t
    
    list(u = u, gamma = gamma, alpha_index = alpha_index)
}


find_Hstar_x <- function(Hstar, x) {
    alpha_index <- Hstar$alpha_index

    x_ext <- numeric(length(x) + 1)
    x_ext[-alpha_index] <- x

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

find_beta <- function(alpha, nbasis, k, alpha_index = 1) {
    component <- find_alpha_components(nbasis, k)
    Hstar_list <- list()
    beta <- matrix(nrow = nbasis, ncol = k)
        
    for(i in 1:k) {
        alpha_i <- alpha[component == i]

        if(alpha_index < 1 || alpha_index > length(alpha_i)) {
            stop("alpha_index must be between 1 and the length of each alpha block")
        }
        
        if(i == 1)
            beta_i <- alpha_i
        else
            beta_i <- find_beta_i(alpha_i, Hstar_list)

        beta[,i] <- beta_i

        if(i < k) {
            Hstar_list[[i]] <- find_Hstar(alpha_i, alpha_index)
        }
    }
    beta
}

householder_diagnostic <- function(alpha, nbasis, k,
                                   alpha_index = 1,
                                   eps = .Machine$double.eps) {
    if(k <= 1) {
        return(list(
            min_ratio = Inf,
            ratios = numeric(0),
            alpha_norms = numeric(0),
            alpha_selected = numeric(0),
            alpha_index = alpha_index
        ))
    }
    alpha_list <- split_alpha(alpha, nbasis, k)
    
    

    alpha_used <- alpha_list[seq_len(k - 1)]

    alpha_norms <- vapply(alpha_used, function(a) sqrt(sum(a^2)), numeric(1))
    alpha_selected <- vapply(alpha_used, function(a) a[alpha_index], numeric(1))

    ratios <- abs(alpha_selected) / pmax(alpha_norms, eps)

    list(
        min_ratio = min(ratios),
        ratios = ratios,
        alpha_norms = alpha_norms,
        alpha_selected = alpha_selected,
        alpha_index = alpha_index
    )
}

householder_ci_diagnostic <- function(mod, eps = .Machine$double.eps) {
    if(mod$k <= 1) {
        return(list(
            min_z_to_boundary = Inf,
            z_to_boundary = numeric(0),
            selected_alpha = numeric(0),
            selected_se = numeric(0),
            alpha_index = mod$alpha_index
        ))
    }
    
    V <- if(!is.null(mod$var_par)) {
             mod$var_par
         } else {
             solve(-mod$hessian)
         }

    se <- sqrt(pmax(diag(V), eps))

    alpha_start <- length(mod$beta0) + 1

    alpha_components <- split(
        seq_along(mod$alpha),
        find_alpha_components(mod$basis$nbasis, mod$k)
    )

    selected_alpha_positions <- vapply(seq_len(mod$k - 1), function(j) {
        alpha_start + alpha_components[[j]][mod$alpha_index] - 1
    }, numeric(1))

    selected_alpha <- mod$par[selected_alpha_positions]
    selected_se <- se[selected_alpha_positions]

    z_to_boundary <- abs(selected_alpha) / pmax(selected_se, eps)

    list(
        min_z_to_boundary = min(z_to_boundary),
        z_to_boundary = z_to_boundary,
        selected_alpha = selected_alpha,
        selected_se = selected_se,
        alpha_index = mod$alpha_index
    )
}

find_alpha_from_beta <- function(beta, nbasis, k, alpha_index = 1) {
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
            Hstar_list[[i]] <- find_Hstar_mat(x, alpha_index)
        }
    }

    unlist(alpha_list, use.names = FALSE)
}

choose_alpha_index_from_beta <- function(beta, nbasis, k) {
    if(k <= 1) {
        return(1)
    }

    candidates <- seq_len(nbasis - k + 1)

    scores <- sapply(candidates, function(alpha_index) {
        alpha <- find_alpha_from_beta(beta, nbasis = nbasis,
                                      k = k, alpha_index = alpha_index)

        householder_diagnostic(alpha, nbasis = nbasis, k = k,
                               alpha_index = alpha_index)$min_ratio
    })

    candidates[which.max(scores)]
}

par_from_beta_parameterisation <- function(beta0, beta, lsigma,
                                           nbasis, k, alpha_index) {
    alpha <- find_alpha_from_beta(beta,
                                  nbasis = nbasis,
                                  k = k,
                                  alpha_index = alpha_index)

    c(beta0, alpha, lsigma)
}

start_from_fit_beta <- function(fit) {
    list(beta0 = fit$beta0,
         beta = fit$beta,
         lsigma = fit$lsigma)
}

maybe_switch_alpha_index_start <- function(beta0, beta, lsigma,
                                           basis, k,
                                           alpha_index = 1,
                                           alpha_tol = 1e-5,
                                           auto_alpha = TRUE) {
    nbasis <- basis$nbasis
    
    if(alpha_index > basis$nbasis - k + 1) {
        alpha_index <- 1
    }
    
    alpha <- find_alpha_from_beta(beta,
                                  nbasis = nbasis,
                                  k = k,
                                  alpha_index = alpha_index)

    if(auto_alpha && k > 1) {
        diag <- householder_diagnostic(alpha,
                                       nbasis = nbasis,
                                       k = k,
                                       alpha_index = alpha_index)

        if(diag$min_ratio < alpha_tol) {
            alpha_index <- choose_alpha_index_from_beta(beta, nbasis, k)

            alpha <- find_alpha_from_beta(beta,
                                          nbasis = nbasis,
                                          k = k,
                                          alpha_index = alpha_index)
        }
    }

    list(par0 = c(beta0, alpha, lsigma),
         alpha_index = alpha_index)
}
