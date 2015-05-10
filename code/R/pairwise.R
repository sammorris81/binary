rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)

source("auxfunctions.R")
source("updateModel.R")
source("mcmc.R")

set.seed(7483)  # site
ns    <- 3000
s     <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0, 1, 0.2), seq(0, 1, 0.2))
x     <- matrix(1, ns, 1)

set.seed(3282)  # data
data <- rRareBinarySpat(x, s = s, knots = knots, beta = 0, xi = 0.1,
                        alpha = 0.5, rho = 0.2, prob.success = 0.01)

y <- data$y
thresh <- data$thresh

dw2 <- rdist(s, knots)
alpha <- 0.5
rho <- 0.2
xi <- 0.1
beta <- -thresh

par.true <- c(alpha, rho, xi, beta)
pairwise.rarebinary(par.true, y, dw2, x)

pairwise.rarebinary <- function(par, y, dw2, x) {
  # par: parameter vector (alpha, rho, xi, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)
  # beta: regresion coefficients
  # xi: shape parameter
  # rho: bandwidth for kernel
  # alpha: spatial dependence

  npars <- length(par)
  alpha <- par[1]
  rho   <- par[2]
  xi    <- par[3]
  beta  <- par[4:npars]

  ns     <- length(y)
  np     <- ncol(x)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  x.beta <- x %*% beta             # should be ns x np
  z      <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (alpha < 0 | alpha > 1) {
    return(1e99)
  }
  if (rho < 0) {
    return(1e99)
  }
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- 0
  for (i in 1:(ns - 1)) {
    for (j in (i + 1):ns) {
      # in all cases, we need the joint
      # first take the kernel bits and add i and j for each knot
      joint <- -sum(apply(kernel[c(i, j), ], 2, sum)^alpha)
      if (y[i] == 0 & y[j] == 0) {
        ll <- ll + joint
      } else if (y[i] == 1 & y[j] == 0) {  # add in marginal for i
        ll <- ll + log(exp(-1 / z[j]) - exp(joint))
      } else if (y[i] == 0 & y[j] == 1) {
        ll <- ll + log(exp(-1 / z[i]) - exp(joint))
      } else if (y[i] == 1 & y[j] == 1) {  # do 1,1 likelihood
        ll <- ll + log(1 - exp(-1 / z[i]) - exp(-1 / z[j]) + exp(joint))
      }
    }
  }

  return(-ll)
}

optim(c(0.5, 0.5, 0, -4), pairwise.rarebinary, y=y, dw2=dw2, x=x)
system.time(optim(c(0.5, 0.5, 0, -4), pairwise.rarebinary, y=y, dw2=dw2, x=x))