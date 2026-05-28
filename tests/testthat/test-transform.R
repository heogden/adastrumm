test_that("Transform gives orthogonal columns", {
    nbasis <- 5
    k <- 3

    components <- find_alpha_components(nbasis, k)

    
    set.seed(1)
    alpha <- rnorm(length(components))

    beta <- find_beta(alpha, nbasis, k)

    expect_equal(sum(beta[,2] * beta[,1]), 0)
    expect_equal(sum(beta[,3] * beta[,1]), 0)
    expect_equal(sum(beta[,3] * beta[,2]), 0)
})

test_that("alpha to beta to alpha round trip works", {
    set.seed(1)

    for(nbasis in c(5, 8)) {
        for(k in 1:3) {
            for(alpha_index in seq_len(nbasis - k + 1)) {
                n_alpha <- sum(nbasis - 0:(k - 1))
                alpha <- rnorm(n_alpha)

                beta <- find_beta(alpha, nbasis, k, alpha_index)
                alpha2 <- find_alpha_from_beta(beta, nbasis, k, alpha_index)
                beta2 <- find_beta(alpha2, nbasis, k, alpha_index)

                expect_equal(beta2, beta, tolerance = 1e-8)
                expect_equal(alpha2, alpha, tolerance = 1e-8)
            }
        }
    }
})
