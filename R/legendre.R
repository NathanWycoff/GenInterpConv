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

#' Polynomial Basis Conversion
#'
#' Convert from a Legendre polynomial representation to a canonical/monomial one.
#' 
#' Translation of Mathematica code form J.M here: https://math.stackexchange.com/questions/86298/the-relationship-between-legendre-polynomials-and-monomial-basis-polynomials
#' Implements Clenshaw's algorithm.
#'
#' @param lc The Legendre coefficients.
#' @return A vector of length length(lc) giving the monomial basis of the input polynomial.
leg_to_mon <- function(lc) {
    n <- length(lc)-1#Poly degree

    pc <- c(lc[length(lc)], rep(0, n))

    z <- rep(0, n+1)
    v <- w <- 0

    for (j in n:1) {
        w <- pc[1]
        pc[1] <- lc[j] - v * z[1]
        z[1] = w
        for (i in 2:(n-j+2)) {
            w <- pc[i]
            pc[i] <- (2 * j - 1) * z[i-1] / j - v * z[i]
            z[i] <- w
        }
        v <- (j-1) / j
    }

    return(pc)
}

#' Polynomial Basis Conversion
#'
#' Convert from a polynomial in monomial / canonical form to one in Legendre form.
#' 
#' Translation of Mathematica code form J.M here: https://math.stackexchange.com/questions/86298/the-relationship-between-legendre-polynomials-and-monomial-basis-polynomials
#' Implements Clenshaw's algorithm.
#'
#' @param pc The monomial coefficients.
#' @return A vector of length length(lc) giving the Legendre basis of the input polynomial.
mon_to_leg <- function(pc) {
    n <- length(pc) - 1
    q <- lc <- rep(0, n+1)

    lc[1] <- pc[n]
    lc[2] <- pc[length(pc)]
    for (k in 2:n) {
        q[1] <- lc[1]
        lc[1] <- pc[n-k+1] + lc[2] / 3
        if (k-1 >= 2) {
            for (j in 2:(k-1)) {
                q[j] <- lc[j]
                lc[j] <- j * lc[j+1] / (2*j+1) + (j-1) * q[j-1]/(2*j-3)
            }
        }
        q[k] <- lc[k]
        lc[k] <- (k-1) * q[k-1] / (2*k-3)
        lc[k+1] <- k * q[k] / (2*k-1)
    }

    return(lc)
}
