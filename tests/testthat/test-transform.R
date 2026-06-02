test_that("Transform gives orthogonal columns", {
    nbasis <- 5
    k <- 3

    components <- find_alpha_components(nbasis, k)

    
    set.seed(1)
    alpha <- rnorm(length(components))

    for(psi_index in seq_len(nbasis - k + 1)) {
        beta <- find_beta(alpha, nbasis, k, psi_index = psi_index)
        
        expect_equal(sum(beta[,2] * beta[,1]), 0)
        expect_equal(sum(beta[,3] * beta[,1]), 0)
        expect_equal(sum(beta[,3] * beta[,2]), 0)
    }
    
})

test_that("alpha to beta to alpha round trip works", {
    set.seed(1)

    for(nbasis in c(5, 8)) {
        for(k in 1:3) {
            for(psi_index in seq_len(nbasis - k + 1)) {
                n_alpha <- sum(nbasis - 0:(k - 1))
                alpha <- rnorm(n_alpha)

                beta <- find_beta(alpha, nbasis, k, psi_index)
                alpha2 <- find_alpha_from_beta(beta, nbasis, k, psi_index)
                beta2 <- find_beta(alpha2, nbasis, k, psi_index)

                expect_equal(beta2, beta, tolerance = 1e-8)
                expect_equal(alpha2, alpha, tolerance = 1e-8)
            }
        }
    }
})


test_that("jac_beta_alpha matches numerical Jacobian", {
    set.seed(1)

    for(nbasis in c(5, 8)) {
        for(k in 1:3) {
            nbasis <- case$nbasis
            k <- case$k

            n_alpha <- sum(nbasis - 0:(k - 1))
            alpha <- rnorm(n_alpha, sd = 0.2)

            for(psi_index in seq_len(nbasis - k + 1)) {
                J_ad <- jac_beta_alpha(
                    alpha = alpha,
                    K = k,
                    n_B = nbasis,
                    psi_index = psi_index
                )

                J_num <- numDeriv::jacobian(
                                       func = function(alpha_in) {
                                           as.vector(find_beta(
                                               alpha = alpha_in,
                                               nbasis = nbasis,
                                               k = k,
                                               psi_index = psi_index
                                           ))
                                       },
                                       x = alpha
                                   )

                expect_equal(J_ad, J_num, tolerance = 1e-5)
            }
        }
    }
   
})
