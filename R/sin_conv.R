#!/usr/bin/Rscript
#  R/sin_conv.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 03.30.2019

## See if we can get the function to converge to the sin function in some sense.
library(MASS)
source('R/legendre.R')
source('R/mindiff.R')
source('R/sin_coefs.R')

set.seed(123)
N <- 10
P <- 100
if (P > length(sin_coefs)) {
    stop("Get more sine coefs pls.")
}
sigma <- 1

# Generate params
beta <- rnorm(P+1)

# Generate x obs and test points
x <- seq(-1, 1, length.out = N)
xx <- seq(-1, 1, length.out = 1000)

# Generate response
mu <- sin(x)
y <- mu + rnorm(N,0,sigma)

# Calculate polynomial matrices
X_legendre <- do.call(cbind, lapply(0:P, function(p) sapply(x, function(xi) leg_poly(p,xi))))

# Obtain solution vectors
beta_mind2 <- mindiff_interp2(X_legendre, y)
pred_expansion <- c(beta_mind2, rep(0, length(sin_coefs)-P-1))

# Plot error at each stage
plot(abs(pred_expansion - sin_coefs))
