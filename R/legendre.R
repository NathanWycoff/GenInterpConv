#!/usr/bin/Rscript
#  R/legendre.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.19.2019

#' Calculate a degree P Legendre polynomial evaluated at x
leg_poly <- function(P, x) {
    if (x > 1 | x < -1) {
        stop("Legendre Polynomial defined on [-1,1]")
    }
    n <- P
    ret <- 0
    for (m in 0:floor(P/2)) {
        if (m %% 2 == 0) {
            sgn <- 1
        } else {
            sgn <- -1
        }
        num <- factorial(2*(n-m))
        denom <-  2^n * factorial(m) * factorial(n-m) * factorial(n-2*m) 
        ret <- ret + sgn*num/denom*x^(n-2*m)
    }
    return(ret)
}
