# load packages and source files
rm(list=ls())
options(warn=2)
library(fields)
library(evd)
# library(spBayes)
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

# get the datasets
load("./simdata.RData")

# data setting and sets to include - written by bash script
setMKLthreads(1)
sets <- c(16:20)
setting <- 3

# extract the relevant setting from simdata
y <- simdata[[setting]]$y
s <- simdata[[setting]]$s
x <- simdata[[setting]]$x

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
y.o    <- matrix(y[obs, ], ntrain, nsets)
X.o    <- matrix(x[obs], ntrain, 1)
s.o    <- s[obs, ]
y.p    <- matrix(y[!obs, ], ntest, nsets)
X.p    <- matrix(x[!obs, ], ntest, 1)
s.p    <- s[!obs, ]
dw2.o  <- rdist(s.o, knots)
d.o    <- as.matrix(rdist(s.o))
diag(d.o) <- 0

####################################################################
#### Start MCMC setup: Most of this is used for the spBayes package
####################################################################
iters <- 20000; burn <- 10000; update <- 1000; thin <- 1
# iters <- 100; burn <- 50; update <- 10; thin <- 1
timings <- rep(NA, 3)

for (i in sets) {
  filename <- paste("sim-results/", setting, "-", i, ".RData", sep = "")
  tblname  <- paste("sim-tables/", setting, "-", i, ".txt", sep ="")
  y.i.o <- matrix(y.o[, i], ntrain, 1)
  y.i.p <- matrix(y.p[, i], ntest, 1)
  cat("Starting: Set", i, "\n")
  
  cat("  Start gev \n")
  
  # spatial GEV
  cat("    Start mcmc fit \n")
  mcmc.seed <- i * 10
  set.seed(mcmc.seed)
  
  fit.gev <- spatial_GEV(y = y.i.o, s = s.o, x = X.o, knots = knots, 
                         beta.init = log(-log(1 - mean(y.o))),
                         beta.mn = 0, beta.sd = 10,
                         beta.eps = 0.1, beta.attempts = 50, 
                         xi.init = 0, xi.mn = 0, xi.sd = 0.5, xi.eps = 0.01, 
                         xi.attempts = 50, xi.fix = TRUE, 
                         a.init = 10, a.eps = 0.2, a.attempts = 50, 
                         a.cutoff = 0.1, b.init = 0.5, b.eps = 0.2, 
                         b.attempts = 50, alpha.init = 0.5, alpha.attempts = 50, 
                         a.alpha.joint = TRUE, alpha.eps = 0.0001,
                         rho.init = 0.1, logrho.mn = -2, logrho.sd = 1, 
                         rho.eps = 0.1, rho.attempts = 50, threads = 1, 
                         iters = iters, burn = burn, 
                         update = update, thin = 1, thresh = 0)
  
  cat("    Start mcmc predict \n")
  post.prob.gev <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = update)
  timings[1] <- fit.gev$minutes
  
  bs.gev <- BrierScore(post.prob.gev, y.i.p)
  print(bs.gev * 100)
  
  # copy table to tables folder on beowulf
  bs <- rbind(bs.gev)
  write.table(bs, file = tblname)
  upload.cmd <- paste("scp ", tblname, " samorris@hpc.stat.ncsu.edu:~/rare-binary/markdown/sim-hmc-2/sim-tables", sep = "")
  system(upload.cmd)
  
  # spatial probit
  cat("  Start probit \n")
  
  cat("    Start mcmc fit \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  fit.probit <- probit(Y = y.i.o, X = X.o, s = s.o, knots = knots, 
                       iters = iters, burn = burn, update = update)
  
  cat("    Start mcmc predict \n")
  post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.p, knots = knots,
                               start = 1, end = iters - burn, update = update)
  timings[2] <- fit.probit$minutes
  
  bs.pro <- BrierScore(post.prob.pro, y.i.p)
  print(bs.pro * 100)
  
  # copy table to tables folder on beowulf
  bs <- rbind(bs.gev, bs.pro)
  write.table(bs, file = tblname)
  upload.cmd <- paste("scp ", tblname, " samorris@hpc.stat.ncsu.edu:~/rare-binary/markdown/sim-hmc-2/sim-tables", sep = "")
  system(upload.cmd)
  
  # spatial logit
  cat("  Start logit \n")
  cat("    Start mcmc fit \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  fit.logit <- spatial_logit(Y = y.i.o, s = s.o, eps = 0.1, 
                             a = 1, b = 1, knots = knots, 
                             iters = iters, burn = burn, update = update)
  
  cat("    Start mcmc predict \n")
  post.prob.log <- pred.splogit(mcmcoutput = fit.logit, s.pred = s.p, 
                                knots = knots, start = 1, end = iters - burn, 
                                update = update)
  timings[3] <- fit.logit$minutes
  
  bs.log <- BrierScore(post.prob.log, y.i.p)
  print(bs.log * 100)
  
  # copy table to tables folder on beowulf
  bs <- rbind(bs.gev, bs.pro, bs.log)
  write.table(bs, file = tblname)
  upload.cmd <- paste("scp ", tblname, " samorris@hpc.stat.ncsu.edu:~/rare-binary/markdown/sim-hmc-2/sim-tables", sep = "")
  system(upload.cmd)
  
  cat("Finished: Set", i, "\n")
  save(fit.gev, post.prob.gev, bs.gev,
       fit.probit, post.prob.pro, bs.pro, 
       fit.logit, post.prob.log, bs.log,
       y.i.p, y.i.o, s, timings,
       file = filename)
}
