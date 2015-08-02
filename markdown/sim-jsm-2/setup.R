#### Exploring really small bandwidth and large number of knots

rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)
# library(microbenchmark)  # comment out for beowulf
library(mvtnorm)
library(Rcpp)
library(numDeriv)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

source("../../code/R/auxfunctions.R", chdir = TRUE)
source("../../code/R/updateModel.R")
source("../../code/R/mcmc.R")
source("../../code/R/probit.R", chdir=T)

set.seed(7483)  # site
ns      <- 1000
nsettings <- 8
s       <- cbind(runif(ns), runif(ns))
knots   <- expand.grid(seq(0, 1, length=21), seq(0, 1, length=21))
knots   <- as.matrix(knots) 
knots.log <- expand.grid(seq(0, 1, length=21), seq(0, 1, length=21))
knots.log <- as.matrix(knots.log)
knots.h <- abs(knots[1, 1] - knots[2, 1])
x       <- matrix(1, ns, 1)

alpha.t <- 0.3
rho.t   <- 0.05
xi.t    <- 0
prob.t  <- c(0.025, 0.050)
int.log <- log(prob.t / (1 - prob.t))

nsets <- 50
set.seed(3282)  # data
y   <- array(data = NA, dim = c(ns, nsets, nsettings))
int <- matrix(data = NA, nrow = nsets, ncol = nsettings)
for (i in 1:nsets) {
  data <- rRareBinarySpat(x, s = s, knots = knots, beta = 0, xi = xi.t,
                          alpha = alpha.t, rho = rho.t, 
                          prob.success = prob.t[1])

  y[, i, 1] <- data$y
  int[i, 1] <- -data$thresh
  
  data <- rRareBinarySpat(x, s = s, knots = knots, beta = 0, xi = xi.t,
                          alpha = alpha.t, rho = rho.t, 
                          prob.success = prob.t[2])
  
  y[, i, 2] <- data$y
  int[i, 2] <- -data$thresh
  
  data <- rLogitSpat(x = x, s = s, knots = knots, beta = int.log[1], 
                     rho = rho.t, sigma.sq = 1, nu = 0.5)
  
  y[, i, 3] <- data$y
  int[i, 3] <- int.log[1]
  
  data <- rLogitSpat(x = x, s = s, knots = knots, beta = int.log[2], 
                     rho = rho.t, sigma.sq = 1, nu = 0.5)
  
  y[, i, 4] <- data$y
  int[i, 4] <- int.log[2]
  
  data <- rProbitSpat(x, s = s, knots = knots, beta = 0, rho = rho.t, 
                      sigma.sq = 1, nu = 0.5, prob.success = prob.t[1])
  
  y[, i, 5] <- data$y
  int[i, 5] <- -data$thresh
  
  data <- rProbitSpat(x, s = s, knots = knots, beta = 0, rho = rho.t, 
                      sigma.sq = 1, nu = 0.5, prob.success = prob.t[2])
  
  y[, i, 6] <- data$y
  int[i, 6] <- -data$thresh
  
  data <- rHotSpotSpat(x, s = s, knots = knots, rho = rho.t, 
                       prob.success = prob.t[1])
  y[, i, 7] <- data$y
  
  data <- rHotSpotSpat(x, s = s, knots = knots, rho = rho.t, 
                       prob.success = prob.t[2])
  y[, i, 8] <- data$y
  
  print(paste("Dataset ", i, " finished", sep = ""))
}

plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "")
points(s[which(y[, 1, 8] != 1), ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[which(y[, 1, 8] == 1), ], pch = 21, col = "firebrick4", bg = "firebrick1")

save(y, x, s, knots, knots.log, knots.h, alpha.t, rho.t, xi.t, prob.t, int, 
     file = "simdata.RData")


