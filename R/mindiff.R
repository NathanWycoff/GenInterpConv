#!/usr/bin/Rscript
#  R/mindiff.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.29.2019

source('R/diff.R')

#' Find the polynomial interpolant with the least derivative norm.
#' @param A A matrix giving the value of the legendre interpolants at the observed locations. Each row represents an observation.
#' @param y The response at the locations of A.
mindiff_interp <- function(A, y) {
    N <- nrow(A)
    d <- ncol(A)-1

    if (N > d) {
        stop("The system should be undertermined to make use of mindiff_interp")
    }

    D <- dmat(d)

    P <- diag(d+1) - t(A) %*% ginv(t(A))

    DtD <- crossprod(D)

    #return(P %*% (-ginv(P %*% DtD %*% P) %*% P %*% DtD %*% ginv(A) %*% y))
    return((-ginv(P %*% DtD %*% P) %*% P %*% DtD %*% ginv(A) %*% y))
}

#' Find the polynomial interpolant with the least derivative norm.
#' @param A A matrix giving the value of the legendre interpolants at the observed locations. Each row represents an observation.
#' @param y The response at the locations of A.
mindiff_interp2 <- function(A, y) {
    N <- nrow(A)
    d <- ncol(A)-1

    if (N > d) {
        stop("The system should be undertermined to make use of mindiff_interp")
    }

    D <- dmat(d)

    DtDgi <- ginv(crossprod(D))
    M <- DtDgi %*% t(A)
    return(M %*% (ginv(A %*% M) %*% y))
}
