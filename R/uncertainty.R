find_u_sample_cluster <- function(cluster, sigma, data, f0_x, f_x) {
    inc <- which(data$c == cluster)
    f0_x_c <- f0_x[inc]
    f_x_c <- f_x[inc, , drop = FALSE]
    y_c <- data$y[inc]

    k <- ncol(f_x)
    I <- diag(1, nrow = k, ncol = k)
    A <- crossprod(f_x_c) / sigma^2 + I
    b <- crossprod(f_x_c, y_c - f0_x_c) / sigma^2
    
    u_hat <- as.numeric(solve(A, b))
    mvnfast::rmvn(1, u_hat, solve(A))
}



find_sample <- function(id, mod) {
    set.seed(id)
    V <- 0.5 * (mod$var_par + t(mod$var_par))
    par <- as.numeric(mvnfast::rmvn(1, mod$par, V))

    par_split <- split_par(par, mod$basis$nbasis)
    
    f0_x <- mod$basis$X %*% par_split$beta0

    if(mod$k > 0) {
        beta <- find_beta(par_split$alpha, mod$basis$nbasis, mod$k, mod$alpha_index)
        f_x <- mod$basis$X %*% beta
    } else {
        beta <- matrix(nrow = mod$basis$nbasis, ncol = 0)
        f_x <- matrix(nrow = nrow(mod$data), ncol = 0)
    }

    clusters <- sort(unique(mod$data$c))
    comps <- lapply(clusters, find_u_sample_cluster,
                    sigma = mod$sigma, data = mod$data, f0_x = f0_x, f_x = f_x)
    u <- Reduce(rbind, comps)
    rownames(u) <- clusters
    
    find_par_cluster(par_split$beta0, beta, u)
}

#' Find parametric bootstrap samples from a fitted model
#'
#' @param mod The fitted model, output from fit_adastrumm.
#' @param n_samples The number of samples to generate.
#'
#' @export
find_samples <- function(mod, n_samples) {
    lapply(1:n_samples, find_sample, mod = mod)
}
