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
knots   <- expand.grid(seq(0, 1, length=41), seq(0, 1, length=41))
# knots <- expand.grid(seq(0, 1, length=15), seq(0, 1, length=15))
knots   <- as.matrix(knots)
knots.h <- abs(knots[1, 1] - knots[2, 1])
x       <- matrix(1, ns, 1)

alpha.t <- 0.3
rho.t   <- 0.025
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
}
plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "")
points(s[which(data$y[, 1] != 1), ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[which(data$y[, 1] == 1), ], pch = 21, col = "firebrick4", bg = "firebrick1")

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
iters <- 50000; burn <- 40000; update <- 1000; thin <- 1
# iters <- 100; burn <- 50; update <- 10; thin <- 1
# setup for spGLM
n.report     <- update
batch.length <- 100
n.batch      <- floor(iters / batch.length)
verbose <- TRUE
tuning <- list("phi" = 0.1, "sigma.sq" = 0.2, "beta" = 1, "w" = 5)
starting <- list("phi" = 3/0.5, "sigma.sq" = 50, "beta" = 0, "w" = 0)
priors <- list("beta.norm" = list(0, 100),
               "phi.unif" = c(0.1, 1e4), "sigma.sq.ig" = c(1, 1))
cov.model <- "exponential"
amcmc <- list("n.batch" = n.batch, "batch.length" = batch.length,
              "accept.rate" = 0.35)

for (i in 7:10) {
  filename <- paste("sim-results/sim-v4-", i, "-1.RData", sep = "")
  y.i.o <- y.o[, i, drop = FALSE]
  y.i.p <- y.validate[, i, drop = FALSE]
  print(paste("Starting: Set ", i, sep = ""))

  print("  start fit.pcl")
  # fit alpha and rho
  print("    start pcl fit")
  fit.pcl <- fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                               alpha.init = 0.5, rho.init = knots.h,
                               xi.fix = TRUE, alpha.fix = FALSE,
                               rho.fix = FALSE, beta.fix = TRUE,
                               y = y.i.o, dw2 = dw2.o, d = d.o,
                               cov = X.o, max.dist = 3 * knots.h,
                               alpha.min = 0.1, alpha.max = 0.9,
                               threads = 2)

  # spatial GEV
  print("    start mcmc fit")
  mcmc.seed <- i * 10
  set.seed(mcmc.seed)

  if (fit.pcl$par[1] < 0.3) {
    alpha.init <- 0.25
  } else {
    alpha.init <- fit.pcl$par[1]
  }
  fit.gev <- mcmc.gev(y = y.i.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit.pcl$beta, beta.m = 0, beta.s = 100,
                  xi.init = 0, xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.05, alpha.m = fit.pcl$par[1], alpha.s = 0.05,
                  rho.tune = 0.1, logrho.m = log(fit.pcl$par[2]), logrho.s = 2,
                  A.tune = 1, A.cutoff = 5 * fit.pcl$par[2],
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 300, rho.attempts = 100,
                  A.attempts = 100, spatial = TRUE, rho.init = fit.pcl$par[2],
                  rho.upper = 9, alpha.init = alpha.init, a.init = 1000,
                  iterplot = FALSE, alpha.fix = FALSE, rho.fix = FALSE,
                  xibeta.joint = FALSE, xi.fix = TRUE,
                  iters = iters, burn = burn, update = update, thin = 1)

  print("    start mcmc predict")
  post.prob.gev <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                                s.pred = s.p, knots = knots,
                                start = 1, end = iters - burn, update = update)

  # spatial logit
  print("  start logit")

  print("    start mcmc fit")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  # does not converge very well, but adaptive takes 3 days per dataset
  fit.logit <- spGLM(formula = y.i.o ~ 1, family = "binomial",
                     coords = s.o, knots = knots, starting = starting,
                     tuning = tuning, priors = priors,
                     cov.model = cov.model, n.samples = iters,
                     verbose = verbose, n.report = update)

  print("    start mcmc predict")
  yp.sp.log <- spPredict(sp.obj = fit.logit, pred.coords = s.p,
                         pred.covars = X.p, start = burn + 1,
                         end = iters, thin = 1, verbose = TRUE,
                         n.report = 500)

  post.prob.log <- t(yp.sp.log$p.y.predictive.samples)

  # spatial probit
  print("  start probit")

  print("    start mcmc fit")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  fit.probit <- probit(Y = y.i.o, X = X.o, s = s.o, knots = knots,
                       iters = iters, burn = burn, update = update)

  print("    start mcmc predict")
  post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.p, knots = knots,
                               start = 1, end = iters - burn, update = update)

  print(paste("Finished: Set ", i, sep = ""))
  save(fit.pcl, fit.gev, post.prob.gev,
       fit.logit, post.prob.log,
       fit.probit, post.prob.pro,
       y.i.p, y.i.o, s,
       file = filename)
}
