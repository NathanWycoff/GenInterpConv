#!/usr/bin/Rscript
#  R/tests.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.21.2019

source('R/L2_norm.R')
source('R/diff.R')

## Test L2 norm
f1 <- function(x) {
    1 + x + x^2 + x^3
}
f2 <- function(x) {
    sin(x)
}


norm_L2(f1)
8/3
dnorm_L2(f1, N = 1e6)
f1(1) - f1(-1)

norm_L2(f2)
0
dnorm_L2(f2)
2*sin(1)

## Test differentiation matrices.
P <- 30
D <- dmat(P)
c <- rnorm(P+1)
dc <- D %*% c

f <- function(x) sapply(0:P, function(p) leg_poly(p, x)) %*% c
fp <- function(x) sapply(0:P, function(p) leg_poly(p, x)) %*% dc

h <- 1e-6
x <- 0.8
(f(x + h) - f(x)) / h
fp(x)
