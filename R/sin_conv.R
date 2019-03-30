#!/usr/bin/Rscript
#  R/sin_conv.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 03.30.2019

## See if we can get the function to converge to the sin function in some sense.
library(MASS)
source('R/legendre.R')
source('R/mindiff.R')
source('R/sin_coefs.R')

set.seed(123)
N <- 10
P <- 2*N
if (P > length(sin_coefs)) {
    stop("Get more sine coefs pls.")
}
sigma <- 0.1

# Generate x obs and test points
x <- seq(-1, 1, length.out = N)
xx <- seq(-1, 1, length.out = 1000)

# Generate response
mu <- sin(x)
y <- mu + rnorm(N,0,sigma)

# Calculate polynomial matrices
X_legendre <- do.call(cbind, lapply(0:P, function(p) sapply(x, function(xi) leg_poly(p,xi))))
XX_legendre <- do.call(cbind, lapply(0:P, function(p) sapply(xx, function(xxi) leg_poly(p,xxi))))

# Obtain solution vectors
beta_mind2 <- mindiff_interp2(X_legendre, y)
leg_pred_expansion <- c(beta_mind2, rep(0, length(sin_coefs)-P-1))
mon_pred_expansion <- c(leg_to_mon(beta_mind2), rep(0, length(sin_coefs)-P-1))
pred_mind2 <- XX_legendre %*% beta_mind2

## In the monomial basis
# Plot error at each stage
par(mfrow=c(1,2))
plot(abs(mon_pred_expansion - sin_coefs))
abline(v = P, col = 'blue')
abline(v = N, col = 'red')
plot(abs(mon_pred_expansion - sin_coefs)[1:5])
abline(v = P, col = 'blue')
abline(v = N, col = 'red')

# Do we every figure out the linear term?
mon_pred_expansion[2]

#### In the legendre basis
# Plot error at each stage
sin_leg <- mon_to_leg(sin_coefs)
par(mfrow=c(2,2))
plot(abs(leg_pred_expansion - sin_leg))
abline(v = P, col = 'blue')
abline(v = N, col = 'red')
plot(abs(leg_pred_expansion - sin_leg)[1:5])
abline(v = P, col = 'blue')
abline(v = N, col = 'red')
## Plot the points/fit to confirm interpolation
plot(x, y, ylim = 1.1*range(pred_mind2))
points(xx, pred_mind2, type = 'l', col = 'orange')
points(xx, sapply(xx, sin), type = 'l', col = 'red')

fl_norm(leg_pred_expansion - sin_leg)
