rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)
library(microbenchmark)
library(mvtnorm)
library(Rcpp)
library(numDeriv)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

source("../../code/R/auxfunctions.R", chdir = TRUE)

# sourceCpp(file = "./getKernel.cpp")
set.seed(2020)
ns <- 1500; nknots <- 144; nt <- 1
wz.star <- array(rnorm(ns * nknots * nt), dim=c(ns, nknots, nt))
a <- matrix(rnorm(nknots*nt), nknots, nt)
IDs <- vector(mode = "list", length = ns)
for (i in 1:ns) {
  IDs[[i]] <- 1:144
}
temp1 <- getKernelCPPwithID(wz_star = wz.star, a = a, nt = nt, ns = ns, nknots = nknots, IDs = IDs)
temp2 <- getKernel(wz.star = wz.star, a = a, IDs = IDs)
temp3 <- getKernelCPP(wz_star = wz.star, a = a)
temp4 <- getKernel(wz.star = wz.star, a = a)
mean(temp1 / temp2)
sd(temp1 / temp2)
mean(temp3 / temp4)
sd(temp3 / temp4)

microbenchmark(getKernelCPP(wz_star = wz.star, a = a),
               getKernel(wz.star = wz.star, a = a, IDs = IDs), 
               getKernelCPPwithID(wz_star = wz.star, a = a, nt = nt, ns = ns, nknots = nknots, IDs = IDs), 
               getKernel(wz.star = wz.star, a = a))

sourceCpp(file = "./getKernel.cpp")
set.seed(2020)
ns <- 1500; nt <- 1
knots <- expand.grid(seq(0, 1, length=41), seq(0, 1, length=41))
nknots <- nrow(knots)
wz.star <- array(rnorm(ns * nknots * nt), dim=c(ns, nknots, nt))
a <- matrix(rnorm(nknots*nt), nknots, nt)
s <- cbind(runif(1500), runif(1500))
d <- rdist(s, knots)^2
d[d < 0.0001] <- 0
IDs1 <- IDs2 <- vector(mode = "list", length = 1500)  # tells us which knots to include
bw <- 0.025
for (i in 1:ns) {
  IDs1[[i]] <- which(d[i, ] < 5 * bw)
  IDs2[[i]] <- 1:nknots
}

microbenchmark(getKernelCPP(wz_star = wz.star, a = a),
               getKernelCPPwithID(wz_star = wz.star, a = a, nt = nt, ns = ns, nknots = nknots, IDs = IDs2),
               getKernel(wz.star = wz.star, a = a),
               getKernel(wz.star = wz.star, a = a, IDs = IDs))

temp1 <- getKernelCPPwithID(wz_star = wz.star, a = a, nt = nt, ns = ns, nknots = nknots, IDs = IDs)
temp2 <- getKernel(wz.star = wz.star, a = a, IDs = IDs)
temp3 <- getKernelCPP(wz_star = wz.star, a = a)
temp4 <- getKernel(wz.star = wz.star, a = a)

mean(temp1 / temp2)
mean(temp3 / temp4)
sd(temp3 / temp4)

sourceCpp(file = "./getKernel.cpp")
z <- abs(matrix(rnorm(ns * nt), ns, nt))
w <- matrix(runif(nknots * ns), ns, nknots)
alpha <- 0.5
tempwzstar1 <- getwzStarCPP(z = z, w = w, alpha = alpha)
tempwzstar2 <- getwzStar(z = z, w = w, alpha = alpha)
mean(tempwzstar1 - tempwzstar2)

microbenchmark(getwzStarCPP(z = z, w = w, alpha = alpha), 
               getwzStar(z = z, w = w, alpha = alpha))

source("../../code/R/auxfunctions.R", chdir = TRUE)
set.seed(2087)
ns <- 1500; nknots <- 144; nt <- 1
alpha <- 0.5
wz.star <- abs(array(rnorm(ns * nknots * nt), dim=c(ns, nknots, nt)))
wz <- wz.star^alpha
a <- abs(matrix(rnorm(nknots*nt), nknots, nt))
a.star <- a^alpha

temp3 <- getKernelCPP(wz = wz, a_star = a.star, alpha = alpha)
temp4 <- getKernel(wz.star = wz.star, a = a)
mean(temp3 / temp4)
sd(temp3 / temp4)

microbenchmark(getKernelCPP(wz = wz, a_star = a.star, alpha = alpha),
               getKernel(wz.star = wz.star, a = a))


set.seed(2087)
z <- abs(matrix(rnorm(100 * 1), 100, 5))
w <- abs(matrix(rnorm(100 * 45), 100, 45))
temp1 <- getwzCPP(z = z, w = w)
temp2 <- getwz(z = z, w = w)
mean(temp1 / temp2)
sd(temp1 / temp2)
microbenchmark(getwzCPP(z = z, w = w), getwz(z = z, w = w))
