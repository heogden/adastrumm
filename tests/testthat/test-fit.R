test_that("sensible fit for test data 1 (straight lines)", {
    data_full <- generate_test_data_1()
    data <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta

    mod <- fit_adastrumm(data)

    expect_equal(mod$k, 2)
    expect_gt(mod$sp, 1000)

    library(tidyverse)
    
    pred_data <- bind_cols(x = data$x, c = data$c, eta = eta) %>%
        mutate(eta_hat = predict_adastrumm(mod, newdata = list(x = x, c = c)))

    rmse <- sqrt(mean(pred_data$eta_hat - pred_data$eta)^2)
    expect_lt(rmse, 0.1)

    pred_data %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = eta_hat)) +
        geom_line(aes(y = eta), linetype = "dashed") +
        facet_wrap(vars(c))

    #' check prediction
    y_hat_data <- predict_adastrumm(mod, newdata = data)
    fitted_y <- fitted_adastrumm(mod)
    expect_equal(y_hat_data, fitted_y)

    
    mu_2_fun <- function(x) {
        predict_adastrumm(mod, newdata = data.frame(x = x, c = 2))
    }

    newdata <- data.frame(x = seq(min(data$x), max(data$x), length = 10),
                          c = 2)
    d_mu_hat_data_man <- numDeriv::grad(mu_2_fun, newdata$x)
    d_mu_hat_data <- predict_adastrumm(mod, newdata = newdata, deriv = TRUE)
    expect_equal(d_mu_hat_data, d_mu_hat_data_man)
    
    #' look at uncertainty
    n_samples <- 1000
    samples <- find_samples(mod, n_samples)

    x_pred_data <- crossing(x = seq(from = min(data$x),
                                    to = max(data$x),
                                    length.out = 100),
                            c = unique(data$c))

    pred_data <- x_pred_data  %>%
        mutate(mu_hat = predict_adastrumm(mod, newdata = list(x = x, c = c), interval = TRUE,
                                      samples = samples),
               d_mu_hat = predict_adastrumm(mod, newdata = list(x = x, c = c), deriv = TRUE,
                                        interval = TRUE, samples = samples)) %>%
        group_by(c) %>%
        mutate(mu = data_full$eta_fun(x, c[1]),
               d_mu = numDeriv::grad(data_full$eta_fun, x = x, c = c[1]))

    pred_data %>%
        filter(c <= 12) %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = mu_hat$estimate)) +
        geom_line(aes(y = mu), colour = "red", linetype = "dashed") + 
        geom_ribbon(aes(ymin = mu_hat$lower, ymax = mu_hat$upper), alpha = 0.2) + 
        geom_point(aes(x = x, y = y), data = data %>% filter(c <= 12)) +
        facet_wrap(vars(c))


    coverage <- as.numeric(pred_data %>%
        mutate(covers = ((mu_hat$lower < mu) & (mu_hat$upper > mu)),) %>%
        ungroup() %>%
        summarise(coverage = mean(covers)))

    expect_gt(coverage, 0.9)
    expect_lt(coverage, 1)

    
    pred_data %>%
        filter(c <= 12) %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = d_mu_hat$estimate)) +
        geom_line(aes(y = d_mu), colour = "red", linetype = "dashed") + 
        geom_ribbon(aes(ymin = d_mu_hat$lower, ymax = d_mu_hat$upper), alpha = 0.2) + 
        facet_wrap(vars(c))

    coverage_d <- as.numeric(pred_data %>%
        mutate(d_covers = ((d_mu_hat$lower < d_mu) & (d_mu_hat$upper > d_mu)),) %>%
        ungroup() %>%
        summarise(d_coverage = mean(d_covers)))

    expect_gt(coverage_d, 0.9)
    expect_lt(coverage_d, 1)
    
})


test_that("sensible fit for test data 2 (not straight lines)", {
    data_full <- generate_test_data_2()
    data <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta


    data %>%
        filter(c <= 10) %>%
        ggplot(aes(x = x, y = y)) +
        geom_point() +
        facet_wrap(vars(c))

    mod <- fit_adastrumm(data)

    expect_lt(mod$sp, 1000)

    library(tidyverse)
    
    pred_data <- bind_cols(x = data$x, c = data$c, eta = eta) %>%
        mutate(eta_hat = predict_adastrumm(mod, newdata = list(x = x, c = c)))

    rmse <- sqrt(mean(pred_data$eta_hat - pred_data$eta)^2)
    expect_lt(rmse, 0.1)

    pred_data %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = eta_hat)) +
        geom_line(aes(y = eta), linetype = "dashed") +
        facet_wrap(vars(c))
})



test_that("gives reasonable fit with tricky blip function", {
    
    g1 <- function(x) {
        dnorm(x, mean = 0.5, sd = 0.05)
    }


    const <- integrate(function(x){g1(x)^2}, lower = 0, upper = 1)$value

    f1 <- function(x) {
        g1(x) / sqrt(const)
    }
    #' so f1 is a normalised version of g


    #' Next, generate data
    
    set.seed(1)

    n <- 100
    u <- rnorm(n)
    n_i <- 10
    sigma <- 0.01
    
    subject <- rep(1:n, each = n_i)
    t <- runif(n * n_i, 0, 1)

    mu <- u[subject] * f1(t)

    epsilon <- rnorm(n * n_i, sd = sigma)
    y <- mu + epsilon

    data <- data.frame(x = t,
                       y = y,
                       c = subject)

    mod <- fit_adastrumm(data, nbasis = 30)


    library(tidyverse)

    pred_data<- crossing(x = seq(0, 1, length.out = 100),
                         c = 1:n) %>%
        mutate(mu = u[c] * f1(x)) %>%
        mutate(mu_hat = predict_adastrumm(mod, newdata = list(x = x, c = c)))

    
    rmse <- sqrt(mean(pred_data$mu_hat - pred_data$mu)^2)
    expect_lt(rmse, 0.1)

    pred_data %>%
        filter(c <= 20) %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = mu_hat)) +
        geom_line(aes(y = mu), linetype = "dashed") +
        facet_wrap(vars(c))
    
    
})

test_that("fits the sleepstudy data", {
    library(tidyverse)
    
    data <- lme4::sleepstudy %>%
        as_tibble %>%
        mutate(c = as.integer(as.factor(Subject)),
               y = Reaction,
               x = Days)

    nbasis <- 15

    expect_no_error(mod <- fit_adastrumm(data, nbasis = nbasis))


    x_pred_data <- crossing(x = seq(from = min(data$x),
                                    to = max(data$x),
                                    length.out = 100),
                            c = unique(data$c))

    expect_no_error(pred_data <- x_pred_data  %>%
        mutate(mu_hat = predict_adastrumm(mod, newdata = list(x = x, c = c))))
    
    pred_data %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = mu_hat)) +
        geom_point(aes(x = x, y = y), data = data) +
        facet_wrap(vars(c))
    
    
})



test_that("fitted values and predictions don't depend on cluster ordering", {
    library(tidyverse)
    
     #' modified from refund::ccb.fpc
    #' obtain a subsample of the data with 25 subjects
    set.seed(1236)
    sample = sample(1:dim(refund::cd4)[1], 25)
    Y.sub = refund::cd4[sample,]

    times <- as.numeric(colnames(Y.sub))
    data <- tibble(c = row(Y.sub)[!is.na(Y.sub)],
                   y = Y.sub[!is.na(Y.sub)],
                   x = times[col(Y.sub)[!is.na(Y.sub)]])

    expect_silent(mod <- fit_adastrumm(data))

    y_hat <- fitted_adastrumm(mod)
    y_hat_pred <- predict_adastrumm(mod, newdata = data)
    expect_equal(y_hat, y_hat_pred)

    pred_with_interval <- predict_adastrumm(mod, newdata = data, interval = TRUE)
    expect_true(all(pred_with_interval$upper > pred_with_interval$estimate))
    
    #' (previously had problems predicting for cluster 8)
    
    y_hat_pred_8 <- y_hat_pred[data$c == 8]
    y_hat_8 <- y_hat[data$c == 8]
    y_8 <- data$y[data$c == 8]

    expect_false(all(y_hat_8 > y_8))
    expect_false(all(y_hat_pred_8 > y_8))


    data_ordered <- data %>%
        arrange(c)

    mod_ordered <- fit_adastrumm(data_ordered)
    y_hat_ordered <- fitted_adastrumm(mod_ordered)
    y_hat_8_ordered <- y_hat_ordered[data_ordered$c == 8]

    expect_equal(y_hat_8, y_hat_8_ordered)
    
    y_hat_pred_ordered <- predict_adastrumm(mod_ordered, newdata = data_ordered)
    y_hat_pred_ordered_8 <- y_hat_pred_ordered[data_ordered$c == 8]

    expect_equal(y_hat_pred_8, y_hat_pred_ordered_8)
})

test_that("fits for first problem data generated from rs model", {
    data_full <- simulate_rs(1, -1, 2, 1, 0.5, 0, 0.1, 20, 5)
    data <- data_full$data

    mod <- fit_adastrumm(data)

    pred_data <- data_full$pred_data
    pred_data$mu_c_hat <- predict_adastrumm(mod, newdata = pred_data)

    expect_lt(mean(abs(pred_data$mu_c_hat - pred_data$mu_c)), 1)
})

test_that("fits for problem data from 1dv model", {
    data_full <- simulate_1dv(22, -0.5, 0.1, 0.5, 0.1, 20, 10)
    data <- data_full$data

    mod <- fit_adastrumm(data)
    
    pred_data <- data_full$pred_data
    pred_data$mu_c_hat <- predict_adastrumm(mod, newdata = pred_data, interval = TRUE)

    expect_lt(mean(abs(pred_data$mu_c_hat$estimate - pred_data$mu_c)), 0.1)

})

test_that("can update model with additional data", {
    data_full <- simulate_1dv(1, -0.5, 0.5, 0.5, 0.1, 50, 5)

    data <- data_full$data

    library(tidyverse)
    data_fit1 <- data %>% filter(c <= 49)

    fit1 <- fit_adastrumm(data_fit1)

    t_full <- system.time(fit_full <- fit_adastrumm(data))
    t_given_fit1 <- system.time(fit_given_fit1 <- update_adastrumm(data,
                                                                   prev_fit = fit1))
    t_given_fit1_fixpar <- system.time(fit_given_fit1_fixpar <-
                                           update_adastrumm(data, prev_fit = fit1,
                                                            fix_pop_par = TRUE))
    expect_lt(t_given_fit1[3], t_full[3])

    expect_lt(t_given_fit1_fixpar[3], t_given_fit1[3])

    pred_data_50 <- data_full$pred_data %>% filter(c == 50)

    preds <- pred_data_50 %>%
        mutate(mu_hat_full = predict_adastrumm(fit_full, newdata = list(x = x, c = c)),
               mu_hat_given_fit1 = predict_adastrumm(fit_given_fit1,
                                                     newdata = list(x = x, c = c)),
               mu_hat_given_fit1_fixpar = predict_adastrumm(fit_given_fit1_fixpar,
                                                            newdata = list(x = x, c = c)))

    expect_lt(mean(abs(preds$mu_hat_full - preds$mu_hat_given_fit1)), 0.1)
    expect_lt(mean(abs(preds$mu_hat_full - preds$mu_hat_given_fit1_fixpar)), 0.1)

})


simulate_2re <- function(seed, n_clusters, n_obs_per_cluster, sigma) {
    set.seed(seed)
    c <- rep(1:n_clusters, each = n_obs_per_cluster)
    x <- runif(length(c), 0, 3*pi)
    
    u1 <- rnorm(n_clusters)
    u2 <- rnorm(n_clusters)
    
    mu_c <- function(x, c) {
        (1 + u1[c]) * (x/2 + sin(x)) + u2[c]
    }
    
    pred_data <- tidyr::crossing(x = seq(min(x), max(x), length.out = 100),
                                 c = 1:n_clusters) %>%
        dplyr::mutate(mu_c = mu_c(x, c))
    
    mu <- (1 + u1[c]) * (x/2 + sin(x)) + u2[c]
    epsilon <- rnorm(length(mu), sd = sigma)
    
    y <- mu + epsilon
    
    data <- tibble(c = c,
                   x = x,
                   y = y,
                   mu = mu)

    list(data = data, pred_data = pred_data)
}

test_that("can find CI in problem case from simulations", {
    data_full <- simulate_2re(69, n_clusters = 50, n_obs_per_cluster = 3, sigma = 0.1)

    data <- data_full$data
    pred_data <- data_full$pred_data

    mod <- fit_adastrumm(data)

    pred_data$mu_c_hat <- predict_adastrumm(mod, newdata = pred_data, interval = TRUE)
    expect_true(all(!is.na(pred_data$mu_c_hat)))
})

test_that("avoid k drop problem in cases from simulations", {
    data_full <- simulate_2re(2, n_clusters = 50, n_obs_per_cluster = 3, sigma = 0.1)
    expect_silent(mod <- fit_adastrumm(data_full$data))

    data_full <- simulate_2re(1, n_clusters = 500, n_obs_per_cluster = 2, sigma = 0.1)
    expect_silent(mod <- fit_adastrumm(data_full$data))
})



make_beta_bad_subspace <- function(nbasis,
                                   k,
                                   bad_row = 1,
                                   scales = seq(k, 1)) {
    Z <- matrix(rnorm((nbasis - 1) * k), nrow = nbasis - 1, ncol = k)
    Q <- qr.Q(qr(Z))[, seq_len(k), drop = FALSE]

    beta <- matrix(0, nrow = nbasis, ncol = k)
    beta[-bad_row, ] <- sweep(Q, 2, scales, `*`)
    beta[bad_row, ] <- 0

    beta
}

simulate_bad <- function(seed,
                         nbasis = 8,
                         k = 3,
                         psi_index_bad = 1,
                         psi_index_good = 2,
                         target_signal_sd = 1,
                         d = 500,
                         n_i = 20,
                         lsigma = -2) {
    set.seed(seed)

    n <- d * n_i
    x_each <- seq(-1, 1, length.out = n_i)
    x <- rep(x_each, times = d)
    c <- rep(seq_len(d), each = n_i)

    basis <- find_orthogonal_spline_basis(nbasis, x)

    ## Construct beta so that the fitted subspace is nearly orthogonal to e_1.
    beta_unscaled <- make_beta_bad_subspace(
        nbasis = nbasis,
        k = k,
        bad_row = psi_index_bad,
        scale = 1
    )

    ## Rescale the signal, without changing the subspace geometry.
    f_x_unscaled <- basis$X %*% beta_unscaled
    scale_fac <- target_signal_sd / sqrt(mean(rowSums(f_x_unscaled^2)))
    beta <- scale_fac * beta_unscaled

    alpha_bad <- find_alpha_from_beta(
        beta,
        nbasis = nbasis,
        k = k,
        psi_index = psi_index_bad
    )

    alpha_good <- find_alpha_from_beta(
        beta,
        nbasis = nbasis,
        k = k,
        psi_index = psi_index_good
    )

    beta0 <- rep(0, nbasis)

    u <- matrix(rnorm(d * k), ncol = k)
    u_ext <- u[c, , drop = FALSE]

    f_x <- basis$X %*% beta
    mu <- rowSums(u_ext * f_x)

    y <- as.vector(mu + rnorm(n, sd = exp(lsigma)))

    list(
        data = data.frame(c = c, x = x, y = y, mu = mu),
        beta = beta,
        beta0 = beta0,
        alpha_bad = alpha_bad,
        alpha_good = alpha_good,
        basis = basis,
        scale_fac = scale_fac
    )
}




test_that("fit_given_start_beta switches psi_index for starting point designed to be close to singular for psi_index = 1", {
    data_full <- simulate_bad(1, d = 200, n_i = 10)
    data <- data_full$data
    basis <- data_full$basis
    sp <- exp(-5)

    fit <- fit_given_start_beta(
        data = data,
        sp = sp,
        k = 3,
        beta0 = data_full$beta0,
        beta = data_full$beta,
        lsigma = -2,
        basis = basis,
        psi_index = 1,
        auto_psi = TRUE,
        psi_tol = 1e-5
    )

    expect_true(fit$psi_index != 1)
    
})



test_that("switching psi_index changes confidence intervals", {
    data_full <- simulate_bad(1, d = 200, n_i = 10)

    data <- data_full$data
    basis <- data_full$basis

    get_mod <- function(psi_index, data, basis) {
        fits <- fits_given_sp(
            sp = exp(-5),
            kmax = 3,
            data = data,
            basis = basis,
            k_tol = 1e-4,
            lambda_tol = 1e-7,
            fits_other_sp = NULL,
            psi_index = psi_index,
            auto_psi = FALSE
        )
        fit <- fits[[4]]
        add_hessian_and_log_ml(fit, basis, data)
    }

    mod_1 <- get_mod(1, data, basis)
    mod_2 <- get_mod(2, data, basis)

    mu_hat_1 <- predict_adastrumm(mod_1,
                                  newdata = data_full$data,
                                  interval = TRUE,
                                  n_samples = 100)
    mu_hat_2 <- predict_adastrumm(mod_2,
                                  newdata = data_full$data,
                                  interval = TRUE,
                                  n_samples = 100)

    mean_width_CI_1 <- mean(mu_hat_1$upper - mu_hat_1$lower)
    mean_width_CI_2 <- mean(mu_hat_2$upper - mu_hat_2$lower)

    diag_1 <- psi_ci_diagnostic(mod_1)
    diag_2 <- psi_ci_diagnostic(mod_2)
    better_on_diag <- which.min(c(diag_1$min_z_to_boundary, diag_2$min_z_to_boundary))

    mean_width_CIs <- c(mean_width_CI_1, mean_width_CI_2)
    expect_gt(max(mean_width_CIs), 0.4)
    expect_lt(min(mean_width_CIs), 0.4)
    expect_gt(mean_width_CIs[better_on_diag], mean_width_CIs[-better_on_diag])
    
 })


test_that("automatic index gives reasonable CI", {
    data_full <- simulate_bad(1, d = 200, n_i = 10)

    data <- data_full$data
    basis <- data_full$basis

    #' do choice of psi_index automatically:
    mod <- fit_adastrumm(data, lsp_poss = -5)

    mu_hat <- predict_adastrumm(mod,
                                newdata = data_full$data,
                                interval = TRUE,
                                n_samples = 100 )
    
    mean_width_CI <- mean(mu_hat$upper - mu_hat$lower)
    expect_lt(mean_width_CI, 0.4)
    
    coverage <- mean(data$mu > mu_hat$lower & data$mu < mu_hat$upper)
    expect_gt(coverage, 0.93)
})

test_that("working log ML has limited parameterisation dependence in test example", {
    data_full <- simulate_bad(1, d = 100, n_i = 5)

    data <- data_full$data
    basis <- data_full$basis
    nbasis <- basis$nbasis

    for(sp in c(exp(-5), 1, exp(3))) {
        fits <- fits_given_sp(
            sp = sp,
            kmax = 3,
            data = data,
            basis = basis,
            k_tol = 1e-4,
            lambda_tol = 1e-7,
            fits_other_sp = NULL,
            psi_index = 1,
            auto_psi = FALSE
        )

        mod_1 <- add_hessian_and_log_ml(fits[[4]], basis, data)

        log_ml_values <- sapply(seq_len(nbasis - mod_1$k + 1), function(psi_index) {
            mod_index <- reparameterise_fit_without_optim(
                mod_1, basis, data, sp, psi_index
            )

            mod_index$log_ml
        })

        expect_lt(diff(range(log_ml_values)), 2)
    }
})


test_that("estimated components are in order of size", {
    library(tidyverse)
    
    #' modified from refund::ccb.fpc
    #' obtain a subsample of the data with 25 subjects
    set.seed(1236)
    sample = sample(1:dim(refund::cd4)[1], 25)
    Y.sub = refund::cd4[sample,]

    times <- as.numeric(colnames(Y.sub))
    data_unnorm <- tibble(c = row(Y.sub)[!is.na(Y.sub)],
                          y = Y.sub[!is.na(Y.sub)],
                          x = times[col(Y.sub)[!is.na(Y.sub)]])

    mod <- fit_adastrumm(data_unnorm, trace = TRUE, lsp_poss = -5)

    expect_equal(mod$lambda, sort(mod$lambda, decreasing = TRUE))
})

test_that("approximate psi CI scores are close enough to exact scores", {
    data_full <- simulate_bad(1, d = 100, n_i = 5)
    data <- data_full$data
    basis <- data_full$basis
    sp <- 1

    fits <- fits_given_sp(
        sp = sp,
        kmax = 3,
        data = data,
        basis = basis,
        k_tol = 1e-4,
        lambda_tol = 1e-7,
        fits_other_sp = NULL,
        psi_index = 1,
        auto_psi = FALSE
    )

    mod <- add_hessian_and_log_ml(fits[[4]], basis, data)

    approx <- approx_psi_ci_scores(mod)

    exact_scores <- vapply(approx$candidates, function(psi_index_new) {
        mod_new <- reparameterise_fit_without_optim(
            fit = mod,
            basis = basis,
            data = data,
            sp = sp,
            psi_index = psi_index_new
        )

        psi_ci_diagnostic(mod_new)$min_z_to_boundary
    }, numeric(1))


    best_approx <- approx$candidates[which.max(approx$scores)]
    current_score <- psi_ci_diagnostic(mod)$min_z_to_boundary
    
    best_exact_score <- max(exact_scores)
    best_approx_exact_score <- exact_scores[approx$candidates == best_approx]

    if(best_approx != mod$psi_index) {
        expect_gt(best_approx_exact_score, current_score)
    }
    expect_gt(best_approx_exact_score, best_exact_score - 0.5)

})

test_that("approx_cov_for_psi_index preserves covariance for current psi_index", {
    data_full <- simulate_bad(1, d = 100, n_i = 5)
    data <- data_full$data
    basis <- data_full$basis
    sp <- 1

    fits <- fits_given_sp(
        sp = sp,
        kmax = 3,
        data = data,
        basis = basis,
        k_tol = 1e-4,
        lambda_tol = 1e-7,
        fits_other_sp = NULL,
        psi_index = 1,
        auto_psi = FALSE
    )

    mod <- add_hessian_and_log_ml(fits[[4]], basis, data)

    approx <- approx_cov_for_psi_index(mod, mod$psi_index)

    expect_equal(approx$psi, mod$psi, tolerance = 1e-8)
    expect_equal(approx$V, mod$var_par, tolerance = 1e-6)
})

test_that("maybe_switch_psi_for_ci_approx improves bad CI diagnostic", {
    data_full <- simulate_bad(1, d = 200, n_i = 10)
    data <- data_full$data
    basis <- data_full$basis
    sp <- 1

    fits <- fits_given_sp(
        sp = sp,
        kmax = 3,
        data = data,
        basis = basis,
        k_tol = 1e-4,
        lambda_tol = 1e-7,
        fits_other_sp = NULL,
        psi_index = 2,
        auto_psi = FALSE
    )

    mod <- add_hessian_and_log_ml(fits[[4]], basis, data)

    diag_before <- psi_ci_diagnostic(mod)

    skip_if(diag_before$min_z_to_boundary >= 2)

    mod2 <- maybe_switch_psi_for_ci_approx(
        fit = mod,
        data = data,
        sp = sp,
        basis = basis,
        psi_ci_tol = 2
    )

    diag_after <- psi_ci_diagnostic(mod2)

    expect_gte(diag_after$min_z_to_boundary,
               diag_before$min_z_to_boundary)

    expect_equal(mod2$beta, mod$beta, tolerance = 1e-8)
    expect_equal(mod2$beta0, mod$beta0, tolerance = 1e-8)
    expect_equal(mod2$lsigma, mod$lsigma, tolerance = 1e-8)
})

test_that("point diagnostic branch can reparameterise without error", {
    data_full <- simulate_bad(1, d = 100, n_i = 5)
    data <- data_full$data
    basis <- data_full$basis
    sp <- 1

    fits <- fits_given_sp(
        sp = sp,
        kmax = 3,
        data = data,
        basis = basis,
        k_tol = 1e-4,
        lambda_tol = 1e-7,
        fits_other_sp = NULL,
        psi_index = 1,
        auto_psi = FALSE
    )

    mod <- add_hessian_and_log_ml(fits[[4]], basis, data)

    mod2 <- maybe_reparameterise_after_hessian(
        fit = mod,
        data = data,
        sp = sp,
        basis = basis,
        psi_tol = Inf,       # force bad_point_chart
        psi_ci_tol = 2,
        auto_psi = TRUE
    )

    expect_equal(mod2$k, mod$k)
    expect_true(is.finite(mod2$l_pen))
    expect_true(is.finite(mod2$log_ml))
})

test_that("Using lambda_tol drops very small components for large sp", {
    set.seed(2)
    data <- data.frame(x = rnorm(100), y = rnorm(100), c = rep(1:10, each = 10))

    lsp_poss <- 7
    mod <- fit_adastrumm(data, lsp_poss = lsp_poss)
    
    suppressMessages(mod_no_lambda_tol <- fit_adastrumm(data, lsp_poss = lsp_poss, lambda_tol = 0))
    
    expect_true(mod$k < mod_no_lambda_tol$k)
})
