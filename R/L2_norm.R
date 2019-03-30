#!/usr/bin/Rscript
#  R/L2_norm.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.21.2019

#' Approximate the L2([-1,1]) norm of a function via l2 norm of a uniform discretization.
#' @param f A function of which the norm is desired.
#' @param N The number of points to evaluate it at.
#' @return The scalar nonnegative L2 norm of the function.
norm_L2 <- function(f, N = 1e3) {
    x <- seq(-1, 1, length.out = N)
    return(2*sum(f(x)) / N)
}

#' Approximate the L2([-1,1]) norm of the derivative of a function via l2 norm of a uniform discretization together with finite differencing.
#' @param f A function of which the derivative's norm is desired.
#' @param N The number of points to evaluate it at.
#' @return The scalar nonnegative L2 norm of the function.
dnorm_L2 <- function(f, N = 1e3, h = 1e-6) {
    x <- seq(-1, 1, length.out = N)
    return(2*sum((f(x+h) - f(x))/h) / N)
}
