# load packages and source files
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

# get the datasets
load("./simdata.RData")

# which setting should run
setting <- 4  
group   <- 2  # sets the start and end sets (grouped by 10)

# extract info about simulation settings
ns        <- dim(y)[1]
nt        <- 1
nsets     <- dim(y)[2]
nsettings <- dim(y)[3]
nknots    <- nrow(knots)

# some precalculated values for quicker pairwise evaluation
dw2     <- rdist(s, knots)
d       <- rdist(s)
diag(d) <- 0

# testing vs training
ntrain <- floor(0.75 * ns)
ntest  <- ns - ntrain
obs    <- c(rep(T, ntrain), rep(F, ntest))
y.o    <- matrix(y[obs, , setting], ntrain, nsets)
X.o    <- matrix(x[obs], ntrain, 1)
s.o    <- s[obs, ]
y.p    <- matrix(y[!obs, , setting], ntest, nsets)
X.p    <- matrix(x[!obs, ], ntest, 1)
s.p    <- s[!obs, ]
dw2.o  <- rdist(s.o, knots)
d.o    <- as.matrix(rdist(s.o))
diag(d.o) <- 0

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
# with so many knots, adaptive is time prohibitive
# amcmc <- list("n.batch" = n.batch, "batch.length" = batch.length,
#               "accept.rate" = 0.35)

start <- (group - 1) * 10 + 1
end   <- group * 10

for (i in start:end) {
  filename <- paste("sim-results/", setting, "-", i, ".RData", sep = "")
  y.i.o <- matrix(y.o[, i], ntrain, 1)
  y.i.p <- matrix(y.p[, i], ntest, 1)
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