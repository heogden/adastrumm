find_u_hat_cluster <- function(cluster, sigma, data, f0_x, f_x, var = FALSE) {
    inc <- which(data$c == cluster)
    f0_x_c <- f0_x[inc]
    f_x_c <- f_x[inc, , drop = FALSE]
    y_c <- data$y[inc]

    k <- ncol(f_x)
    I <- diag(1, nrow = k, ncol = k)
    A <- crossprod(f_x_c) / sigma^2 + I
    b <- crossprod(f_x_c, y_c - f0_x_c) / sigma^2
    u_hat <- as.numeric(solve(A, b))
    if(var) {
        var_u_hat <- solve(A)
        return(list(u_hat = u_hat, var_u_hat = var_u_hat))
    } else {
        return(u_hat)
    }
}


find_u_hat <- function(sigma, data, f0_x, f_x, var = FALSE) {
    clusters <- sort(unique(data$c))
    comps_full <- lapply(clusters, find_u_hat_cluster,
                         sigma = sigma, data = data, f0_x = f0_x, f_x = f_x,
                         var = var)
    if(var) {
        comps <- lapply(comps_full, "[[", "u_hat")
        var_u_hat <- lapply(comps_full, "[[", "var_u_hat")
        
    } else {
        comps <- comps_full
    }
    
    u_hat <- Reduce(rbind, comps)
    rownames(u_hat) <- clusters
    if(var) {
        return(list(u_hat = u_hat, var_u_hat = var_u_hat))
    } else {
        return(u_hat)
    }
}
