#########################################################################
# A simulation study to determine if a thresholded or skew methods improve
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
#   All: s in [0, 6] x [0, 6]
#   1: GEV link
#      a: alpha = 0.3, 100 knots, 1% rareness, bw = 3
#      b: alpha = 0.7, 100 knots, 1% rareness, bw = 3
#      c: alpha = 0.3, 100 knots, 5% rareness, bw = 3
#      d: alpha = 0.7, 100 knots, 5% rareness, bw = 3
#   2: Logit link
#      a: rho = 1, 100 knots, 1% rareness
#      b: rho = 3, 100 knots, 1% rareness
#      c: rho = 1, 100 knots, 5% rareness
#      d: rho = 3, 100 knots, 5% rareness
#   3: Probit link -- Hold off for now
#      a: gamma = 0.9, 100 knots, 1% rareness
#      b: gamma = 0.1, 100 knots, 1% rareness
#      c: gamma = 0.9, 100 knots, 5% rareness
#      d: gamma = 0.1, 100 knots, 5% rareness
#########################################################################

rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)

source("../../R/auxfunctions.R", chdir=T)
source("../../R/updateModel.R", chdir=T)
source("../../R/mcmc.R", chdir=T)
source("../../R/probit.R", chdir=T)

ns        <- 4000
nt        <- 1
nsets     <- 10
nlinks    <- 2   # using Max-stable and logit
nsettings <- 4   # a, b, c, d from above

knots.1 <- knots.2 <- seq(from = 0, to = 6, length = 10)
knots <- expand.grid(knots.1, knots.2)
nknots <- nrow(knots)

X <- matrix(1, ns, nt)

y         <- array(NA, dim=c(ns, nt, nsets, nsettings * nlinks))
a.t       <- array(NA, dim=c(nknots, nt, nsets, nsettings))
w.tilde.t <- array(NA, dim=c(ns, nt, nsets, nsettings))

set.seed <- 7483 # site
s <- cbind(runif(ns), runif(ns))


# max-stable % rareness is set with prob.ms
prob.ms  <- c(0.01, 0.01, 0.05, 0.05)
alpha.ms <- c(0.3, 0.7, 0.3, 0.7)
xi.ms    <- 0.1
rho.ms   <- 3

# to generate logistic data with around 1% rareness, need intercept = -log(99)
# to generate logistic data with around 5% rareness, need intercept = -log(19)
int.logit   <- c(-log(99), -log(99), -log(19), -log(19)) 
rho.logit   <- c(1, 3, 1, 3)
sigsq.logit <- 1
nu.logit    <- 0.5

set.seed <- 3181 # data
for (i in 1:nsets) {
  for (j in 1:nlinks) {
    for (k in 1:nsettings) {
      setting <- (j - 1) * nsettings + k
      if (j == 1) {  # max-stable
        data <- rRareBinarySpat(x = X, s = s, knots = knots, beta = 0, xi.ms, 
                                alpha = alpha.ms[k], rho.ms, 
                                prob.success = prob.ms[k])
        a.t[, , i, k] <- data$a
      } else {  # logit
        data <- rLogitSpat(x = X, s = s, knots = knots, beta = int.logit[k], 
                           rho = rho.logit[k], sigma.sq = sigsq.logit, 
                           nu = nu.logit)
        w.tilde.t[, , i, k] <- data$w.tilde
      }
      y[, , i, setting] <- data$y
    }   
  }
  cat("Set", i, "finished \n") 
}