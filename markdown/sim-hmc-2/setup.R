#### Exploring a couple of hotspots
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

source("../../code/R/spatial_gev.R", chdir = TRUE)
source("../../code/R/spatial_logit.R", chdir = TRUE)
source("../../code/R/spatial_probit.R", chdir = TRUE)

set.seed(7483)  # site
ns      <- c(1000, 2000, 2000, 1000)
nsettings <- 4
nhotspots <- c(6, 6, 3, 3)

knots   <- expand.grid(seq(0, 1, length=21), seq(0, 1, length=21))
knots   <- as.matrix(knots) 
knots.h <- abs(knots[1, 1] - knots[2, 1])

simdata <- vector(mode = "list", length = nsettings)

alpha.t <- 0.3
rho.t   <- 0.05
xi.t    <- 0
prob.t  <- c(0.05, 0.05, 0.025, 0.025)

p <- 0.95  # P(Y=1|hot spot)
q <- 0.01  # P(Y=1|background)
r <- 0.05  # Hot spot radius

nsets <- 100
set.seed(3282)  # data

for (setting in 1:nsettings) {
  simdata[[setting]]$y <- matrix(data = NA, nrow = ns[setting], ncol = nsets)
  simdata[[setting]]$s <- cbind(runif(ns[setting]), runif(ns[setting]))
  simdata[[setting]]$x <- matrix(1, ns[setting], 1)
  simdata[[setting]]$hotspots <- vector(mode = "list", length = nsets)
  for (i in 1:nsets) {
    k  <- rpois(1, nhotspots[setting]) + 1
    c1 <- runif(k)
    c2 <- runif(k)
    
    hot <- rep(FALSE, ns[setting])
    for (j in 1:k) {
      d   <- sqrt((simdata[[setting]]$s[, 1] - c1[j])^2 + 
                    (simdata[[setting]]$s[, 2] - c2[j])^2) 
      hot <- ifelse(d < r, TRUE, hot)
    }
    simdata[[setting]]$y[, i] <- rbinom(ns[setting], 1, ifelse(hot, p, q))
    simdata[[setting]]$hotspots[[i]] <- cbind(c1, c2)
    
    if (i %% 20 == 0) {
      print(paste("  Dataset ", i, " finished", sep = ""))
    }
  }
  print(paste("Setting ", setting, " finished", sep = ""))
}

quartz(width = 10, height = 7)
par(mfrow = c(2, 2))
plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
     main = "Setting 1: 5%, ns = 1000")
points(simdata[[1]]$s[which(simdata[[1]]$y[, 8] != 1), ], pch = 21, 
       col = "dodgerblue4", bg = "dodgerblue1")
points(simdata[[1]]$s[which(simdata[[1]]$y[, 8] == 1), ], pch = 21, 
       col = "firebrick4", bg = "firebrick1")
points(simdata[[1]]$hotspots[, , 8], pch = 21, col = "black", bg = "black")

plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
     main = "Setting 2: 5%, ns = 2000")
points(simdata[[2]]$s[which(simdata[[2]]$y[, 8] != 1), ], pch = 21, 
       col = "dodgerblue4", bg = "dodgerblue1")
points(simdata[[2]]$s[which(simdata[[2]]$y[, 8] == 1), ], pch = 21, 
       col = "firebrick4", bg = "firebrick1")
points(simdata[[2]]$hotspots[, , 8], pch = 21, col = "black", bg = "black")

plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
     main = "Setting 3: 2.5%, ns = 2000")
points(simdata[[3]]$s[which(simdata[[3]]$y[, 8] != 1), ], pch = 21, 
       col = "dodgerblue4", bg = "dodgerblue1")
points(simdata[[3]]$s[which(simdata[[3]]$y[, 8] == 1), ], pch = 21, 
       col = "firebrick4", bg = "firebrick1")
points(simdata[[3]]$hotspots[, , 8], pch = 21, col = "black", bg = "black")

plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
     main = "Setting 4: 2.5%, ns = 1000")
points(simdata[[4]]$s[which(simdata[[4]]$y[, 8] != 1), ], pch = 21, 
       col = "dodgerblue4", bg = "dodgerblue1")
points(simdata[[4]]$s[which(simdata[[4]]$y[, 8] == 1), ], pch = 21, 
       col = "firebrick4", bg = "firebrick1")
points(simdata[[4]]$hotspots[, , 8], pch = 21, col = "black", bg = "black")


mean(simdata[[1]]$y)
mean(simdata[[2]]$y)
mean(simdata[[3]]$y)
mean(simdata[[4]]$y)

save(simdata, knots, knots.h, rho.t, prob.t, file = "simdata.RData")