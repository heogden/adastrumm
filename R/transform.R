find_alpha_components <- function(nbasis, k) {
    n_each <- nbasis - 0:(k-1)
    rep(1:k, times = n_each)
}


split_alpha <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    split(alpha, component)
}

hh_sign <- function(alpha) {
    if (alpha[1] >= 0) 1 else -1
}

find_Hstar_mat <- function(alpha) {
    alpha_norm <- sqrt(sum(alpha^2))

    u <- alpha
    u[1] <- u[1] + hh_sign(alpha) * alpha_norm

    t <- sum(u^2)
    
    gamma <- 2 / t

    H <- diag(nrow = length(u)) - gamma * outer(u, u)
    H[ , -1, drop = FALSE]
}


find_T_list <- function(alpha, nbasis, k) {
    alpha_list <- split_alpha(alpha, nbasis, k)

    T_list <- list()
    # no transformation for k = 1: identity matrix
    T_list[[1]] <- diag(nrow = length(alpha_list[[1]]))
    
    if(k > 1) {
        for(j in 2:k) {
            T_list[[j]] <- T_list[[j-1]] %*% find_Hstar_mat(alpha_list[[j-1]])
        }
    }
    
    T_list
}


## TODO:: could get transform from C++ code instead
find_Hstar <- function(alpha) {
    alpha_norm <- sqrt(sum(alpha^2))

    u <- alpha
    u[1] <- u[1] + hh_sign(alpha) * alpha_norm
    t <- sum(u^2)
    gamma <- 2 / t
    
    list(u = u, gamma = gamma)
}


find_Hstar_x <- function(Hstar, x) {
    a <- sum(Hstar$u[-1] * x)
    c(0, x) - a * Hstar$gamma * Hstar$u
}


find_beta_i <- function(alpha_i, Hstar_list) {
    i <- length(Hstar_list) + 1
    s_j <- alpha_i
    
    for(j in (i-1):1) {
        s_j <- find_Hstar_x(Hstar_list[[j]], s_j)
    }
    
    s_j
}

find_beta <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    Hstar_list <- list()
    beta <- matrix(nrow = nbasis, ncol = k)
        
    for(i in 1:k) {
        alpha_i <- alpha[component == i]
        if(i == 1)
            beta_i <- alpha_i
        else
            beta_i <- find_beta_i(alpha_i, Hstar_list)

        beta[,i] <- beta_i

        if(i < k) {
            Hstar_list[[i]] <- find_Hstar(alpha_i)
        }
    }
    beta
}

