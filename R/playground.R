#!/usr/bin/Rscript
#  R/playground.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.19.2019

library(MASS)
source('R/legendre.R')
source('R/mindiff.R')

set.seed(123)
N <- 20
P <- 40
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
X_naive <- do.call(cbind, lapply(0:P, function(p) x^p))
XX_naive <- do.call(cbind, lapply(0:P, function(p) xx^p))
X_legendre <- do.call(cbind, lapply(0:P, function(p) sapply(x, function(xi) leg_poly(p,xi))))
XX_legendre <- do.call(cbind, lapply(0:P, function(p) sapply(xx, function(xxi) leg_poly(p,xxi))))

# Obtain solution vectors
beta_naive <- ginv(crossprod(X_naive)) %*% t(X_naive) %*% y
beta_legendre <- ginv(crossprod(X_legendre)) %*% t(X_legendre) %*% y
beta_mind <- mindiff_interp(X_legendre, y)
beta_mind2 <- mindiff_interp2(X_legendre, y)

# Obtain predictions
pred_naive <- XX_naive %*% beta_naive
pred_legendre <- XX_legendre %*% beta_legendre
pred_mind <- XX_legendre %*% beta_mind
pred_mind2 <- XX_legendre %*% beta_mind2

# Plot results
plot(x, y, ylim = 1.1*range(c(pred_naive, pred_legendre)))
points(xx, pred_naive, type = 'l', col = 'red')
points(xx, pred_legendre, type = 'l', col = 'blue')
points(xx, pred_mind, type = 'l', col = 'green')
points(xx, pred_mind2, type = 'l', col = 'orange')

plot(x, y, ylim = 1.1*range(pred_mind2))
points(xx, pred_mind2, type = 'l', col = 'orange')
points(xx, sapply(xx, sin), type = 'l', col = 'red')
