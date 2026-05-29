test_that("prediction works with unordered data with repeated x values", {
    library(tidyverse)
    
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

    data <- bind_cols(x_index = x_index,
                      c = rep(1:n_clusters, each = n_obs_per_cluster)) %>%
        mutate(x = x_full[x_index],
               eta = (beta_0 + beta_1 * x + u[c]),
               epsilon = rnorm(length(x_index), 0, sigma),
               y = eta + epsilon)

    mod <- fit_adastrumm(data)
    y_hat <- fitted_adastrumm(mod)
    y_hat_pred <- predict_adastrumm(mod, newdata = data)

    expect_equal(y_hat, y_hat_pred)
})

test_that("prediction reasonable outside of range of data", {
    data_full <- generate_test_data_1()
    data <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta

    mod <- fit_adastrumm(data)

    x <- seq(min(data$x) - 0.1, max(data$x) + 0.1, length.out = 100)
    y_hat <- predict_adastrumm(mod, newdata = data.frame(x = x, c = 1))

    diffs <- y_hat[-1] - y_hat[-length(y_hat)]
    expect_true(all(diffs >= 0))
    
})


test_that("prediction with confidence interval works in cd4 example", {
    library(tidyverse)
    
     #' modified from refund::ccb.fpc
    #' obtain a subsample of the data with 25 subjects
    set.seed(1236)
    sample = sample(1:dim(refund::cd4)[1], 25)
    Y.sub = refund::cd4[sample,]

    times <- as.numeric(colnames(Y.sub))
    data_cd4 <- tibble(c = row(Y.sub)[!is.na(Y.sub)],
                       y = Y.sub[!is.na(Y.sub)],
                       x = times[col(Y.sub)[!is.na(Y.sub)]])
    
    mod <- fit_adastrumm(data_cd4)
    y_hat_pred <- predict_adastrumm(mod, newdata = data_cd4, interval = TRUE)
    
})

test_that("prediction gives error if pred_data is not of correct form", {
    data <- data.frame(x = rnorm(100), y = rnorm(100), c = rep(1:10, each = 10))
    mod <- fit_adastrumm(data)

    expect_error(predict_adastrumm(mod, newdata = NULL),
                 "must contain a column called 'x'")

    expect_error(
        predict_adastrumm(mod, newdata = list(c = 1:20, x = rnorm(20))),
        "contains clusters not present in the fitted model"
    )
})
