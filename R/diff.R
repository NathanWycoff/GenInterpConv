#!/usr/bin/Rscript
#  R/diff.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.29.2019


## Forms the differentiation matrix in the Legendre Basis.
#' @param d The polynomial degree.
dmat <- function(d) {
    alpha <- function(n) (n+1) / (2*n+1)
    gamma <- function(n) n / (2*n+1)

    Q <- matrix(0, nrow = d, ncol = d)
    for (i in 1:d) {
        for (j in 1:i) {
            if (i == j) {
                Q[i,i] <- i / alpha(i-1)
            } else {
                Q[i,j] <- 1 / alpha(i-1) * 
                    (ifelse(j>1, alpha(j-2) * Q[i-1,j-1],0) + 
                     gamma(j) * Q[i-1,j+1] - 
                     ifelse(i>2, gamma(i-1) * Q[i-2,j], 0))
            }
        }
    }

    D <- t(cbind(rbind(0, Q),0))
    return(D)
}
