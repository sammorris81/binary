rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)
library(microbenchmark)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp(file = "./pairwise.cpp")

source("auxfunctions.R")
source("updateModel.R")
source("mcmc.R")

pairwise.rarebinaryR <- function(par, y, dw2, cov) {
  # par: parameter vector (alpha, rho, xi, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

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
      # joint <- -sum(apply(kernel[c(i, j), ], 2, sum)^alpha)
      joint <- -getJoint2(kernel = kernel[c(i, j), ], alpha = alpha)
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

pairwise.rarebinaryCPP <- function(par, y, dw2, x, threads = 1) {
  # par: parameter vector (alpha, rho, xi, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  alpha <- par[1]
  rho   <- par[2]
  xi    <- par[3]
  beta  <- par[4:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  ns     <- length(y)
  np     <- ncol(x)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  x.beta <- x %*% beta             # should be ns x np
  z      <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (alpha < 1e-10 | alpha > (1 - 1e-10)) {
    return(1e99)
  }
  if (rho < 1e-10) {
    return(1e99)
  }
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    threads = threads)

  return(-ll)
}

getJoint2 <- function(kernel, alpha) {
  knot.contribute <- colSums(kernel)^alpha
  # print(knot.contribute)
  joint <- sum(knot.contribute)
  return(joint)
}

set.seed(7483)  # site
ns    <- 50
s     <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0.01, 0.99, length=12), seq(0.01, 0.99, length=12))
x     <- matrix(1, ns, 1)

set.seed(3282)  # data
data <- rRareBinarySpat(x, s = s, knots = knots, beta = 0, xi = 0.25,
                        alpha = 0.9, rho = 0.1, prob.success = 0.05)

y <- data$y
plot(s[which(y == 1), ], pch=19, col=2, ylim=c(0, 1), xlim=c(0, 1))
points(knots)

thresh <- data$thresh
dw2 <- rdist(s, knots)

alpha <- 0.3
rho <- 0.1
xi <- 0.25
beta <- -thresh

par.true <- c(alpha, rho, xi, beta)
pairwise.rarebinary1(par.true, y, dw2, x)
pairwise.rarebinary2(par.true, y, dw2, x, threads = 6)
microbenchmark(pairwise.rarebinary1(par.true, y, dw2, x),
               pairwise.rarebinary2(par.true, y, dw2, x, threads = 6), times = 50)


results <- optim(c(0.5, 0.2, 0, -4), pairwise.rarebinaryCPP, y = y, dw2 = dw2, x = x,
                 threads = 4)

# get gradient and hessian matrices - unnecessary calculations currently causing this to be very slow
library(numDeriv)
d <- matrix(0, 1, 4)
h <- matrix(0, 4, 4)
for (i in 1:(ns - 1)) {
  for (j in i:(ns)) {
    d <- d + jacobian(pairwise.rarebinaryR, x=results$par, y=y[c(i, j)], dw2=dw2, cov=x)
    h <- h + hessian(pairwise.rarebinaryR, x=results$par, y=y[c(i, j)], dw2=dw2, cov=x)
  }
  print(paste("i = ", i))
}




system.time(optim(c(0.5, 0.2, 0, -4), pairwise.rarebinary4, y=y, dw2=dw2, x=x, threads = 4))

# timing
# around 5 seconds when ns = 200
# around 35 seconds when ns = 500
# around 30 seconds when ns = 1000 and nthread = 4





# originally testing speed of 4 different methods
# r with function to get kernel combinations
# r with apply to get kernel combinations
# cpp without pointers used for kernel combinations
# cpp with pointers used for kernel combinations
pairwise.rarebinary1(par.true, y, dw2, x)
pairwise.rarebinary2(par.true, y, dw2, x)
pairwise.rarebinary3(par.true, y, dw2, x, threads = 4)
pairwise.rarebinary4(par.true, y, dw2, x, threads = 4)

microbenchmark(pairwise.rarebinary1(par.true, y, dw2, x),
               pairwise.rarebinary2(par.true, y, dw2, x),
               pairwise.rarebinary3(par.true, y, dw2, x),
               pairwise.rarebinary4(par.true, y, dw2, x), times = 50)

# when ns = 100
# Unit: milliseconds
#                                      expr       min        lq      mean    median        uq       max neval
# pairwise.rarebinary1(par.true, y, dw2, x) 135.24563 136.01825 151.13502 136.77547 152.98440 268.17594    50
# pairwise.rarebinary2(par.true, y, dw2, x) 372.12786 374.45486 400.07225 375.38344 389.75383 577.38213    50
# pairwise.rarebinary3(par.true, y, dw2, x)  15.78985  15.90479  17.25929  16.04462  16.45333  24.27989    50
# pairwise.rarebinary4(par.true, y, dw2, x)  10.69880  10.71617  11.54585  10.79183  11.06266  15.96124    50

# when ns = 200
# Unit: milliseconds
#                                      expr       min        lq      mean    median        uq       max neval
# pairwise.rarebinary1(par.true, y, dw2, x)  540.62006  548.09158  557.58586  550.75619  558.26123  597.28543    50
# pairwise.rarebinary2(par.true, y, dw2, x) 1490.19606 1514.30555 1529.62580 1524.22436 1540.11635 1618.19104    50
# pairwise.rarebinary3(par.true, y, dw2, x)   81.09228   81.92453   84.71521   82.88241   88.15506   92.18456    50
# pairwise.rarebinary4(par.true, y, dw2, x)   42.20400   42.28760   42.59847   42.32965   42.62128   48.97245    50

xplot <- c(100, 200)
y1 <- c(151.135, 557.586)
y2 <- c(400.072, 1529.626)
y3 <- c(17.259, 84.715)
y4 <- c(11.546, 42.598)
ylim <- range(y1, y2, y3, y4)
plot(xplot, y1, type = "l", ylim = ylim)
lines(xplot, y2, lty=2)
lines(xplot, y3, lty=3)
lines(xplot, y4, lty=4)
