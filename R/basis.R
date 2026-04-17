find_orthogonal_spline_basis <- function(nbasis, x) {
    interior_knots <- seq(min(x) - 0.01 * stats::sd(x),
                          max(x) + 0.01 * stats::sd(x),
                          length = nbasis - 2)
    knots <- orthogonalsplinebasis::expand.knots(interior_knots)
    basis <- orthogonalsplinebasis::OrthogonalSplineBasis(knots)
    list(X = orthogonalsplinebasis::evaluate(basis, x = x),
         S = orthogonalsplinebasis::OuterProdSecondDerivative(basis),
         basis = basis,
         nbasis = nbasis)
}

update_basis <- function(basis, x) {
    basis$X <- orthogonalsplinebasis::evaluate(basis$basis, x = x)
    basis
}


# if x outside knots, extrapolate constant at nearest knot point, rather than NA
# (which is the default in orthogonalsplinebasis:::evaluate.SplineBasis
evaluate_with_extrap <- function(object, x, deriv) {
    knots <- object@knots
    order <- object@order
    lower <- knots[order]
    upper <- knots[length(knots)-order+1]

    extrap_lower <- (x < lower)
    extrap_upper <- (x > upper)
    
    x[extrap_lower] <- lower
    x[extrap_upper] <- upper
   
    if(deriv)
        object <- orthogonalsplinebasis::deriv(object)
    
    result <- orthogonalsplinebasis::evaluate(object, x)

    if(deriv) {
        result[extrap_lower, ] <- 0
        result[extrap_upper, ] <- 0
    }

    result
}


find_spline_fun <- function(beta, basis) {
    spline_fun <- function(x, deriv = FALSE) {
        B <-  evaluate_with_extrap(basis$basis, x = x, deriv = deriv)
        result <- tcrossprod(B, t(beta))
        if(!is.matrix(beta))
            result <- as.numeric(result)
        result
    }
    spline_fun
}

