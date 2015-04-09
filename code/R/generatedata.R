# model update testing
# independence:
#   beta: PASS
#   xi: PASS
#   beta and xi: PASS
# spatial dependence:
#   beta: PASS
#   xi: PASS
#   beta and xi: PASS
#   alpha: PASS
#   rho: PASS
#   alpha and rho: PASS (both at alpha = 0.4 and alpha = 0.8)
#       moves fairly slowly in the chain though
#   A:
#   beta, xi, alpha:
#   beta, xi, rho:
#   beta, xi, alpha, rho

# currently trying to figure out what's wrong with the dependent version
# of the model...issues with beta, xi, alpha, and rho updates when data
# are created using the model with dependence.

rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)

source("auxfunctions.R")
source("updateModel.R")
set.seed(10)
ns   <- 2000
nt   <- 1
s    <- cbind(runif(ns, 0, 10), runif(ns, 0, 10))
x <- array(1, dim=c(ns, nt, 3))
x[, , 2] <- s[, 1]
x[, , 3] <- s[, 2]
knots.1 <- seq(0, 10, length=9)
knots.2 <- seq(0, 10, length=9)
knots <- expand.grid(knots.1, knots.2)
# knots  <- s
nknots <- nrow(knots)

set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
source("mcmc.R")
iters <- 100000
burn  <- 80000
xi.t <- 0.1
beta.t <- c(-3, 0, 0)
alpha.t <- 0.1
rho.t   <- 3
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t)

plot(s[, 1], s[, 2], type="p", main=bquote(alpha==.(alpha.t)), xlab="", ylab="")
observed <- which(data$y == 1)
points(s[observed, 1], s[observed, 2], bg="red", pch=21)
mean(data$y)

data <- rRareBinaryInd(x=x, beta=beta.t, xi=xi.t, thresh=0)
plot(s[, 1], s[, 2], type="p", main=bquote(alpha==1), xlab="", ylab="")
observed <- which(data$y == 1)
points(s[observed, 1], s[observed, 2], bg="red", pch=21)
mean(data$y)

