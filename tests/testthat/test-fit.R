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

    sd(mod$u[,1])
    sd(mod$u[,2])
    cor(mod$u[,1], mod$u[,2])

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

test_that("chooses a reasonably large sp if have straight lines", {
    sigma_u <- 0.5
    sigma <- 0.1
    n_clusters <- 50
    n_obs_per_cluster <- 10
    beta_0 <- -1
    beta_1 <- 2

    set.seed(1)

    u <- rnorm(n_clusters, sd = sigma_u)

    x_full <- seq(0, 1, length.out = 1000)
    x_index <- sample(seq_along(x_full), n_clusters * n_obs_per_cluster, replace = TRUE)

    library(tidyverse)

    data <- bind_cols(x_index = x_index,
                      c = rep(1:n_clusters, each = n_obs_per_cluster)) %>%
        mutate(x = x_full[x_index],
               eta = (beta_0 + beta_1 * x + u[c]),
               epsilon = rnorm(length(x_index), 0, sigma),
               y = eta + epsilon)

    eta_data <- crossing(x = x_full,
                         c = 1:n_clusters) %>%
        mutate(eta = (beta_0 + beta_1 * x + u[c]))

    mod <- fit_adastrumm(data)

    expect_gt(mod$sp, 100)
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
    mod <- fit_adastrumm(data, nbasis = nbasis)


    x_pred_data <- crossing(x = seq(from = min(data$x),
                                    to = max(data$x),
                                    length.out = 100),
                            c = unique(data$c))

    pred_data <- x_pred_data  %>%
        mutate(mu_hat = predict_adastrumm(mod, newdata = list(x = x, c = c)))
    
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

    mod <- fit_adastrumm(data)

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
    expect_equal(mod$k, 2)
    
    pred_data <- data_full$pred_data
    pred_data$mu_c_hat <- predict_adastrumm(mod, newdata = pred_data, interval = TRUE)

    expect_lt(mean(abs(pred_data$mu_c_hat$estimate - pred_data$mu_c)), 1)

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
