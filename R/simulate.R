simulate_adastrumm <- function(mod) {
    # simulate normalised data, using original x and c
    x_norm <- mod$data$x
    c <- mod$data$c
    n_c <- length(unique(c))
    u <- matrix(stats::rnorm(n_c * mod$k), ncol = mod$k)
    u_ext <- u[c, , drop = FALSE]
    mu_norm <- unname(as.numeric(mod$f0_x) + rowSums(u_ext * mod$f_x))
    epsilon <- stats::rnorm(length(x_norm), sd = mod$sigma)
    y_norm <- mu_norm + epsilon
    
    x <- mod$norm$m_x + mod$norm$s_x * x_norm
    # output mu as "truth" for simulation studies
    mu <- mod$norm$m_y + mod$norm$s_y * mu_norm
    y <- mod$norm$m_y + mod$norm$s_y * y_norm
    data.frame(c = c,
               x = x,
               mu = mu,
               y = y)
}
