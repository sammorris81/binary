#### Exploring knot spacing and bandwidth
#### Few of options:
#### 1. Keep true knots on a grid, make it smaller, and decrease rho.
#### 2. Keep true knots on a grid, keep it the same size, but decrease rho
#### 3. Adjust sample size
####    n = 1000 for prop = 0.05 and
####    n = 2000 for prop = 0.01

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
ns    <- 1000
s     <- cbind(runif(ns), runif(ns))
knots.t1 <- expand.grid(seq(0, 1, length=20), seq(0, 1, length=20))
knots.t2 <- expand.grid(seq(0, 1, length=10), seq(0, 1, length=10))
x     <- matrix(1, ns, 1)

alpha.t <- 0.3
rho.t   <- 0.01
xi.t    <- 0
prop.t  <- 0.05

set.seed(3282)  # data
nreps <- 10
nsettings <- 2
y <- matrix(data = NA, nrow = ns, ncol = nreps * nsettings)
thresh <- rep(NA, nreps * nsettings)
for (i in 1:nreps) {
  idx <- (i - 1) * nsettings + 1
  data <- rRareBinarySpat(x, s = s, knots = knots.t1, beta = 0, xi = xi.t,
                          alpha = alpha.t, rho = rho.t, prob.success = prop.t)

  y[, idx] <- data$y
  thresh[idx] <- data$thresh

  idx <- idx + 1
  data <- rRareBinarySpat(x, s = s, knots = knots.t2, beta = 0, xi = xi.t,
                          alpha = alpha.t, rho = rho.t, prob.success = prop.t)

  y[, idx] <- data$y
  thresh[idx] <- data$thresh

  idx <- idx + 1
}

# par(mfrow=c(3, 3))
# for (i in 1:3) {
#   idx <- (i - 1) * 3 + 1
#   plot(knots.t1, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "",
#        main = paste("knots 20 x 20, rho =", rho.t))
#   points(s[which(y[, idx] != 1), ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
#   points(s[which(y[, idx] == 1), ], pch = 21, col = "firebrick4", bg = "firebrick1")
#
#   idx <- idx + 1
#   plot(knots.t2, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "",
#        main = paste("knots 10 x 10, rho =", rho.t))
#   points(s[which(y[, idx] != 1), ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
#   points(s[which(y[, idx] == 1), ], pch = 21, col = "firebrick4", bg = "firebrick1")
#
#   idx <- idx + 1
#   plot(knots.t3, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "",
#        main = "500 random knots")
#   points(s[which(y[, idx] != 1), ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
#   points(s[which(y[, idx] == 1), ], pch = 21, col = "firebrick4", bg = "firebrick1")
# }

# knots for fitting - purposely using
knots <- expand.grid(seq(0.01, 0.99, length=12), seq(0.01, 0.99, length=12))
knots.h <- knots[2, 1] - knots[1, 1]
knots <- as.matrix(knots)
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

for (i in 15:16) {
  filename <- paste("sim-results/sim-knots-", i, "-5-1-r.RData", sep = "")
  y.i.o <- y.o[, i, drop = FALSE]
  y.i.p <- y.validate[, i, drop = FALSE]
  print(paste("Starting: Set ", i, sep = ""))

  # fit alpha only with rho set at knot spacing
  print("  start GEV fit")
  print("    start pcl fit")
  fit.gev.pcl <- fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                                   alpha.init = 0.5, rho.init = 0.02,
                                   xi.fix = TRUE, alpha.fix = FALSE,
                                   rho.fix = FALSE, beta.fix = TRUE,
                                   y = y.i.o, dw2 = dw2.o, d = d.o,
                                   cov = X.o, max.dist = 2.5 * knots.h,
                                   alpha.min = 0.1, alpha.max = 0.9,
                                   threads = 3)

  print("    start mcmc gev")
  mcmc.seed <- i * 10
  set.seed(mcmc.seed)
  if (fit.gev.pcl$par[1] < 0.3) {
    alpha.init <- 0.25
  } else {
    alpha.init <- fit.gev.pcl$par[1]
  }
  fit.gev <- mcmc(y = y.i.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit.gev.pcl$beta, beta.m = 0, beta.s = 100,
                  xi.init = 0, xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.05, alpha.m = fit.gev.pcl$par[1],
                  alpha.s = 0.05, rho.tune = 0.1,
                  logrho.m = log(fit.gev.pcl$par[2]), logrho.s = 2, A.tune = 1,
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 300, rho.attempts = 100,
                  spatial = TRUE, rho.init = fit.gev.pcl$par[2], rho.upper = 9,
                  alpha.init = alpha.init, a.init = 1000, iterplot = TRUE,
                  alpha.fix = FALSE, rho.fix = FALSE, xibeta.joint = FALSE,
                  xi.fix = TRUE,
                  iters = iters, burn = burn, update = update, thin = 1)

  print("    start mcmc predict")
  post.prob.gev <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p, s.pred = s.p,
                              knots = knots, start = 1, end = iters - burn,
                              update = update)

#   # spatial logit
#   print("  start logit")
#
#   print("    start mcmc fit")
#   mcmc.seed <- mcmc.seed + 1
#   set.seed(mcmc.seed)
#   fit.logit <- spGLM(formula = y.i.o ~ 1, family = "binomial", coords = s.o,
#                      knots = knots, starting = starting, tuning = tuning,
#                      priors = priors, cov.model = cov.model,
#                      n.samples = iters, verbose = verbose,
#                      n.report = n.report)
#
#   print("    start mcmc predict")
#   yp.sp.log <- spPredict(sp.obj = fit.logit, pred.coords = s.p,
#                          pred.covars = X.p, start = burn + 1, end = iters,
#                          thin = 1, verbose = TRUE, n.report = 500)
#
#   post.prob.log <- t(yp.sp.log$p.y.predictive.samples)
#
#   # spatial probit
#   print("  start probit")
#
#   print("    start mcmc fit")
#   mcmc.seed <- mcmc.seed + 1
#   set.seed(mcmc.seed)
#   fit.probit <- probit(Y = y.i.o, X = X.o, s = s.o, knots = knots,
#                        iters = iters, burn = burn, update = update)
#
#   print("    start mcmc predict")
#   post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
#                                s.pred = s.p, knots = knots,
#                                start = 1, end = iters - burn, update = update)

  print(paste("Finished: Set ", i, sep = ""))
  save(fit.gev.pcl, fit.gev, post.prob.gev,
#        fit.logit, post.prob.log,
#        fit.probit, post.prob.pro,
       y.i.p, y.i.o, s,
       file = filename)
}
