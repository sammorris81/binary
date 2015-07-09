#### Making a few adjustments to the MCMC. In particular, tyring to let alpha
#### and rho not be fixed in the MCMC. For more details, see the section
#### in Sim_take_3.pdf on somewhat informative priors.

#### this is the same as pairwise-sim-2.R, but I made some minor changes to
#### the MCMC to use C and do some site excusion based on distance to the knots.

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
source("../../code/R/updateModel.R")
source("../../code/R/mcmc.R")
source("../../code/R/probit.R", chdir=T)

set.seed(7483)  # site
ns    <- 2000
s     <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0.01, 0.99, length=12), seq(0.01, 0.99, length=12))
knots <- as.matrix(knots)
knots.h <- abs(knots[1, 1] - knots[2, 1])
x     <- matrix(1, ns, 1)

alpha.t <- 0.3
rho.t   <- 0.1
xi.t    <- 0

set.seed(3282)  # data
y <- matrix(data = NA, nrow = ns, ncol = 10)
for (i in 1:10) {
  data <- rRareBinarySpat(x, s = s, knots = knots, beta = 0, xi = xi.t,
                          alpha = alpha.t, rho = rho.t, prob.success = 0.05)

  y[, i] <- data$y
}
# plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "")
# points(s[which(y[, 1] != 1), ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
# points(s[which(y[, 1] == 1), ], pch = 21, col = "firebrick4", bg = "firebrick1")

thresh <- data$thresh
dw2    <- rdist(s, knots)
d      <- rdist(s)
diag(d) <- 0

# testing vs training
ntrain <- floor(0.75 * ns)
ntest  <- ns - ntrain
obs   <- c(rep(T, ntrain), rep(F, ntest))
y.o   <- y[obs, , drop = FALSE]
X.o   <- matrix(x[obs], ntrain, 1)
s.o   <- s[obs, ]
dw2.o <- rdist(s.o, knots)
d.o   <- as.matrix(rdist(s.o))
diag(d.o) <- 0
y.validate <- y[!obs, , drop = FALSE]
X.p <- matrix(x[!obs, ], ntest, 1)
s.p <- s[!obs, ]

####################################################################
#### Start MCMC setup: Most of this is used for the spBayes package
####################################################################
iters <- 45000; burn <- 35000; update <- 1000; thin <- 1
# iters <- 100; burn <- 50; update <- 10; thin <- 1
# setup for spGLM
n.report <- update
verbose <- TRUE
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1,
               "beta"=0.1, "w"=0.1)
starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1,
                 "beta"=0, "w"=0)
priors <- list("beta.norm"=list(1, 100),
               "phi.unif"=c(0.1, 1e4), "sigma.sq.ig"=c(1, 1),
               "tau.sq.ig"=c(1, 1))
cov.model <- "exponential"

for (i in 3:4) {
  filename <- paste("sim-results/pairwise-sim-", i, "-3.RData", sep = "")
  y.i.o <- y.o[, i, drop = FALSE]
  y.i.p <- y.validate[, i, drop = FALSE]
  print(paste("Starting: Set ", i, sep = ""))

  print("  start fit.9")
  # fit alpha and rho
  print("    start pcl fit")
  fit.9 <- fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                             alpha.init = 0.5, rho.init = knots.h,
                             xi.fix = TRUE, alpha.fix = FALSE,
                             rho.fix = FALSE, beta.fix = TRUE,
                             y = y.i.o, dw2 = dw2.o, d = d.o,
                             cov = X.o, max.dist = 2.5 * knots.h,
                             alpha.min = 0.1, alpha.max = 0.9,
                             threads = 3)

  # spatial GEV
  print("    start mcmc fit")
  mcmc.seed <- i * 10
  set.seed(mcmc.seed)

  if (fit.9$par[1] < 0.3) {
    alpha.init <- 0.25
  } else {
    alpha.init <- fit.9$par[1]
  }
  fit.gev.9 <- mcmc(y = y.i.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                    beta.init = fit.9$beta, beta.m = 0, beta.s = 100,
                    xi.init = 0, xi.m = 0, xi.s = 0.5,
                    knots = knots, beta.tune = 1, xi.tune = 0.1,
                    alpha.tune = 0.05, alpha.m = fit.9$par[1], alpha.s = 0.05,
                    rho.tune = 0.1, logrho.m = log(fit.9$par[2]), logrho.s = 2,
                    A.tune = 1, A.cutoff = 5 * fit.9$par[2],
                    beta.attempts = 50, xi.attempts = 50,
                    alpha.attempts = 300, rho.attempts = 100,
                    spatial = TRUE, rho.init = fit.9$par[2], rho.upper = 9,
                    alpha.init = alpha.init, a.init = 1000, iterplot = TRUE,
                    alpha.fix = FALSE, rho.fix = FALSE, xibeta.joint = FALSE,
                    xi.fix = TRUE,
                    iters = iters, burn = burn, update = update, thin = 1)

  print("    start mcmc predict")
  post.prob.gev.9 <- pred.spgev(mcmcoutput = fit.gev.9, x.pred = X.p,
                                s.pred = s.p, knots = knots, A.cutoff = 5 * fit.9$par[2],
                                start = 1, end = iters - burn, update = update)

  # fit alpha only with rho set at knot spacing
  print("  start fit.10")
  print("    start pcl fit")
  fit.10 <- fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                                   alpha.init = 0.5, rho.init = knots.h,
                                   xi.fix = TRUE, alpha.fix = FALSE,
                                   rho.fix = TRUE, beta.fix = TRUE,
                                   y = y.i.o, dw2 = dw2.o, d = d.o,
                                   cov = X.o, max.dist = 2.5 * knots.h,
                                   alpha.min = 0.1, alpha.max = 0.9,
                                   threads = 3)

  print("    start mcmc fit 1")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  if (fit.10$par[1] < 0.3) {
    alpha.init <- 0.25
  } else {
    alpha.init <- fit.10$par[1]
  }
  fit.gev.10 <- mcmc(y = y.i.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                     beta.init = fit.10$beta, beta.m = 0, beta.s = 100,
                     xi.init = 0, xi.m = 0, xi.s = 0.5,
                     knots = knots, beta.tune = 1, xi.tune = 0.1,
                     alpha.tune = 0.05, alpha.m = fit.10$par[1], alpha.s = 0.05,
                     rho.tune = 0.1, logrho.m = log(knots.h), logrho.s = 2,
                     A.tune = 1, A.cutoff = 5 * knots.h,
                     beta.attempts = 50, xi.attempts = 50,
                     alpha.attempts = 300, rho.attempts = 100,
                     spatial = TRUE, rho.init = knots.h, rho.upper = 9,
                     alpha.init = alpha.init, a.init = 1000, iterplot = TRUE,
                     alpha.fix = FALSE, rho.fix = FALSE, xibeta.joint = FALSE,
                     xi.fix = TRUE,
                     iters = iters, burn = burn, update = update, thin = 1)

  print("    start mcmc predict")
  post.prob.gev.10 <- pred.spgev(mcmcoutput = fit.gev.10, x.pred = X.p,
                                 s.pred = s.p, knots = knots, A.cutoff = 5 * knots.h,
                                 start = 1, end = iters - burn, update = update)

  print("    start mcmc fit 2")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  if (fit.10$par[1] < 0.3) {
    alpha.init <- 0.25
  } else {
    alpha.init <- fit.10$par[1]
  }
  fit.gev.10a <- mcmc(y = y.i.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                     beta.init = fit.10$beta, beta.m = 0, beta.s = 100,
                     xi.init = 0, xi.m = 0, xi.s = 0.5,
                     knots = knots, beta.tune = 1, xi.tune = 0.1,
                     alpha.tune = 0.05, alpha.m = fit.10$par[1], alpha.s = 0.05,
                     rho.tune = 0.1, logrho.m = log(knots.h), logrho.s = 2,
                     A.tune = 1, A.cutoff = 5 * knots.h,
                     beta.attempts = 50, xi.attempts = 50,
                     alpha.attempts = 300, rho.attempts = 100,
                     spatial = TRUE, rho.init = knots.h, rho.upper = 9,
                     alpha.init = alpha.init, a.init = 1000, iterplot = TRUE,
                     alpha.fix = FALSE, rho.fix = TRUE, xibeta.joint = FALSE,
                     xi.fix = TRUE,
                     iters = iters, burn = burn, update = update, thin = 1)

  print("    start mcmc predict")
  post.prob.gev.10a <- pred.spgev(mcmcoutput = fit.gev.10a, x.pred = X.p,
                                  s.pred = s.p, knots = knots, A.cutoff = 5 * knots.h,
                                  start = 1, end = iters - burn,
                                  update = update)

  print(paste("Finished: Set ", i, sep = ""))
  save(fit.9, fit.gev.9, post.prob.gev.9,
       fit.10, fit.gev.10, post.prob.gev.10,
       fit.gev.10a, post.prob.gev.10a,
       y.i.p, y.i.o, s,
       file = filename)
}
