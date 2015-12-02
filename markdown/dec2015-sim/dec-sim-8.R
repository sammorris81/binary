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
library(pROC)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

source("../../code/R/spatial_gev.R", chdir = TRUE)
source("../../code/R/spatial_logit.R", chdir = TRUE)
source("../../code/R/spatial_probit.R", chdir = TRUE)

# get the datasets
load("./simdata.RData")

# data setting and sets to include - written by bash script
# setMKLthreads(1)
setting <- 2
sets <- 11:20

# extract the relevant setting from simdata
y <- simdata[[setting]]$y
s <- simdata[[setting]]$s
x <- simdata[[setting]]$x

# extract info about simulation settings
ns     <- dim(y)[1]
nt     <- 1
nknots <- 500

# some precalculated values for quicker pairwise evaluation
dw2     <- rdist(s, knots)
d       <- rdist(s)
diag(d) <- 0

# testing vs training
if (setting %in% c(1, 3, 5)) {
  ntrain <- 500
} else {
  ntrain <- 1000
}

ntest   <- ns - ntrain
ntest.1 <- floor(0.05 * ntest)
ntest.0 <- ntest - ntest.1

####################################################################
#### Start MCMC setup: Most of this is used for the spBayes package
####################################################################
iters <- 25000; burn <- 15000; update <- 1000; thin <- 1
# iters <- 100; burn <- 50; update <- 10; thin <- 1
n.report     <- 10
batch.length <- 100
n.batch      <- floor(iters / batch.length)
verbose      <- TRUE
tuning       <- list("phi" = 0.1, "sigma.sq" = 0.2, "beta" = 1, "w" = 5)
starting     <- list("phi" = 3/0.5, "sigma.sq" = 50, "beta" = 0, "w" = 0)
priors       <- list("beta.norm" = list(0, 100),
                     "phi.unif" = c(0.1, 1e4), "sigma.sq.ig" = c(1, 1))
cov.model <- "exponential"
timings   <- rep(NA, 3)
# with so many knots, adaptive is time prohibitive
amcmc     <- list("n.batch" = n.batch, "batch.length" = batch.length,
                  "accept.rate" = 0.35)

upload.pre <- "samorris@hpc.stat.ncsu.edu:~/rare-binary/markdown/dec2015-sim/sim-tables"

timings <- rep(NA, 3)

for (i in sets) {
  # get the sites where we will be predicting
  pred.0 <- sample(which(y[, i] == 0), size = ntest.0)
  pred.1 <- sample(which(y[, i] == 1), size = ntest.1)
  obs    <- rep(T, ns)
  obs[c(pred.0, pred.1)] <- F
  
  y.i.o  <- matrix(y[obs, i], ntrain, 1)
  X.o    <- matrix(x[obs], ntrain, 1)
  s.i.o  <- s[obs, , i]
  y.i.p  <- matrix(y[!obs, i], ntest, 1)
  X.p    <- matrix(x[!obs], ntest, 1)
  s.i.p  <- s[!obs, , i]
  
  ntrain.0 <- sum(y.i.o == 0) - ntest.0
  ntrain.1 <- sum(y.i.o == 1) - ntest.1
  
  filename <- paste("sim-results/", setting, "-", i, ".RData", sep = "")
  tblname  <- paste("sim-tables/", setting, "-", i, ".txt", sep ="")
  
  #### Knot setup
  # for ns.o = 500, using all sites
  # 
  # for ns.o = 1000, using a space filling design with 500 knots
  # There is a bug in cover.design when you have set nn = FALSE. So, to avoid 
  # the warning, setting number of nearest neighbors to 0 and turning off 
  # nearest neighbors. This will take a bit longer, but ultimately comes back 
  # with a similar design
  ####
  if (setting %in% c(1, 3, 5)) {
    knots.i.o <- s.i.o
  } else {
    set.seed(setting * 10 + i)
    knots.i.o.0 <- cover.design(R = s.i.o[y.i.o == 0, ], nd = floor(nknots * 0.95), 
                                nruns = 1, nn = FALSE, num.nn = 0)$design
    knots.i.o.1 <- cover.design(R = s.i.o[y.i.o == 1, ], nd = floor(nknots * 0.05), 
                                nruns = 1, nn = FALSE, num.nn = 0)$design
    knots.i.o <- rbind(knots.i.o.0, knots.i.o.1)
  }
  cat("Starting: Set", i, "\n")
  
  cat("  Start gev \n")
  
  # spatial GEV
  cat("    Start mcmc fit \n")
  mcmc.seed <- i * 10
  set.seed(mcmc.seed)
  
  fit.gev <- spatial_GEV(y = y.i.o, s = s.i.o, x = X.o, knots = knots.i.o, 
                         beta.init = log(-log(1 - mean(y.i.o))),
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
                              s.pred = s.i.p, knots = knots.i.o,
                              start = 1, end = iters - burn, update = update)
  timings[1] <- fit.gev$minutes
  
  bs.gev <- BrierScore(post.prob.gev, y.i.p)
  post.prob.gev.med <- apply(post.prob.gev, 2, median)
  roc.gev <- roc(y.i.p ~ post.prob.gev.med)
  auc.gev <- roc.gev$auc
  
  print(bs.gev * 100)
  
  # copy table to tables folder on beowulf
  bs <- rbind(c(bs.gev, auc.gev))
  rownames(bs) <- "gev"
  colnames(bs) <- c("bs", "auc")
  write.table(bs, file = tblname)
  upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
  system(upload.cmd)
  
  # spatial probit
  cat("  Start probit \n")
  
  cat("    Start mcmc fit \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  fit.probit <- probit(Y = y.i.o, X = X.o, s = s.i.o, knots = knots.i.o, 
                       iters = iters, burn = burn, update = update)
  
  cat("    Start mcmc predict \n")
  post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.i.p, knots = knots.i.o,
                               start = 1, end = iters - burn, update = update)
  timings[2] <- fit.probit$minutes
  
  bs.pro <- BrierScore(post.prob.pro, y.i.p)
  post.prob.pro.med <- apply(post.prob.pro, 2, median)
  roc.pro <- roc(y.i.p ~ post.prob.pro.med)
  auc.pro <- roc.pro$auc
  
  print(bs.pro * 100)
  
  # copy table to tables folder on beowulf
  bs <- rbind(bs, c(bs.pro, auc.pro))
  rownames(bs)[2] <- "pro"
  write.table(bs, file = tblname)
  upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
  system(upload.cmd)
  
  # spatial logit
  print("  start logit")
  print("    start mcmc fit")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  tic       <- proc.time()[3]
  fit.logit <- spGLM(formula = y.i.o ~ 1, family = "binomial",
                     coords = s.i.o, knots = knots.i.o, starting = starting,
                     tuning = tuning, priors = priors,
                     cov.model = cov.model, n.samples = iters,
                     verbose = verbose, n.report = n.report, amcmc = amcmc)
  
  print("    start mcmc predict")
  yp.sp.log <- spPredict(sp.obj = fit.logit, pred.coords = s.i.p,
                         pred.covars = X.p, start = burn + 1,
                         end = iters, thin = 1, verbose = TRUE,
                         n.report = 500)
  
  post.prob.log <- t(yp.sp.log$p.y.predictive.samples)
  toc        <- proc.time()[3]
  timings[3] <- toc - tic
  
  bs.log <- BrierScore(post.prob.log, y.i.p)
  post.prob.log.med <- apply(post.prob.log, 2, median)
  roc.log <- roc(y.i.p ~ post.prob.log.med)
  auc.log <- roc.log$auc
  
  print(bs.log * 100)
  
  # copy table to tables folder on beowulf
  bs <- rbind(bs, c(bs.log, auc.log))
  rownames(bs)[3] <- "log"
  write.table(bs, file = tblname)
  upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
  system(upload.cmd)
  
  cat("Finished: Set", i, "\n")
  save(fit.gev, post.prob.gev, bs.gev, roc.gev, auc.gev,
       fit.probit, post.prob.pro, bs.pro, roc.pro, auc.pro,
       fit.logit, post.prob.log, bs.log, roc.log, auc.log,
       y.i.p, y.i.o, knots.i.o, s.i.o, s.i.p, timings,
       file = filename)
}
