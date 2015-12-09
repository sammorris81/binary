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
setting <- 1
sets    <- 16:20
knot.design <- 3

# extract the relevant setting from simdata
y <- simdata[[setting]]$y
s <- simdata[[setting]]$s
x <- simdata[[setting]]$x
do.upload <- TRUE

# extract info about simulation settings
ns     <- dim(y)[1]
nt     <- 1
nknots <- 441

# testing vs training
ntrain <- 1000

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

upload.pre <- "samorris@hpc.stat.ncsu.edu:~/rare-binary/markdown/dec2015-knots/sim-tables"

for (i in sets) {
  # storage for some of the results
  bs <- matrix(NA, 9, 2)  # place to store brier scores and auc
  rownames(bs) <- c("gev-1", "gev-2", "gev-3", 
                    "probit-1", "probit-2", "probit-3",
                    "logit-1", "logit-2", "logit-3")
  colnames(bs) <- c("bs", "auc")
  bs.gev  <- bs.pro  <- bs.log  <- rep(NA, 3)
  auc.gev <- auc.pro <- auc.log <- rep(NA, 3)
  roc.gev <- roc.pro <- roc.log <- vector(mode = "list", length = 3)
  
  timings <- matrix(NA, 3, 3)
  rownames(timings) <- c("gev", "probit", "logit")
  colnames(timings) <- c("knots 1", "knots 2", "knots 3")
  
  # start the simulation
  set.seed(setting * 10 + i)
  
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
  
  filename <- paste("sim-results/", setting, "-", i, "-", knot.design, ".RData", 
                    sep = "")
  tblname  <- paste("sim-tables/", setting, "-", i, "-", knot.design, ".txt", 
                    sep ="")
  
  #### Knot setup
  # knot design 1: 21 x 21 grid
  # knot design 2: SRS of 441 sites
  # knot design 3: stratified sample of 441 sites where the breakdown matches
  #                % of sites with 1s and 0s from the training set
  #
  # There is a bug in cover.design when you have set nn = FALSE. So, to avoid 
  # the warning, setting number of nearest neighbors to 0 and turning off 
  # nearest neighbors. This will take a bit longer, but ultimately comes back 
  # with a similar design
  ####
  knots.1 <- as.matrix(expand.grid(x = seq(0, 1, length = 21), 
                                   y = seq(0, 1, length = 21)))
  knots.2 <- cover.design(R = s.i.o, nd = nknots, nruns = 1, nn = FALSE, 
                          num.nn = 0)$design
  
  phat.o   <- mean(y.i.o)  # what proportion of the sites in training are 1s
  nknots.1 <- floor(phat.o * nknots)
  nknots.0 <- floor((1 - phat.o) * nknots)
  knots.3.0 <- cover.design(R = s.i.o[y.i.o == 0, ], nd = nknots.0, nruns = 1, 
                            nn = FALSE, num.nn = 0)$design
  knots.3.1 <- cover.design(R = s.i.o[y.i.o == 1, ], nd = nknots.1, nruns = 1, 
                            nn = FALSE, num.nn = 0)$design
  knots.3 <- rbind(knots.3.0, knots.3.1)
  
  if (knot.design == 1) {
    knots <- knots.1
  } else if (knot.design == 2) {
    knots <- knots.2
  } else {
    knots <- knots.3
  }
  
  cat("Starting: Set", i, "\n")
  
  #### spatial GEV
  cat("  Start gev \n")
  
  cat("    Start mcmc fit - Knots", knot.design, " \n")
  mcmc.seed <- i * 10
  set.seed(mcmc.seed)
  
  fit.gev <- spatial_GEV(y = y.i.o, s = s.i.o, x = X.o, knots = knots, 
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
                              s.pred = s.i.p, knots = knots,
                              start = 1, end = iters - burn, update = update)
  timings[1, knot.design] <- fit.gev$minutes
  
  bs.gev[knot.design] <- BrierScore(post.prob.gev, y.i.p)
  post.prob.gev.med <- apply(post.prob.gev, 2, median)
  roc.gev[[knot.design]] <- roc(y.i.p ~ post.prob.gev.med)
  auc.gev[knot.design] <- roc.gev[[knot.design]]$auc
  
  print(bs.gev * 100)
  
  # copy table to tables folder on beowulf
  bs.row <- knot.design
  bs[bs.row, ] <- c(bs.gev[knot.design], auc.gev[knot.design])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ###### spatial probit
  cat("  Start probit \n")
  
  cat("    Start mcmc fit - Knots", knot.design, " \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  fit.probit <- probit(Y = y.i.o, X = X.o, s = s.i.o, knots = knots, 
                       iters = iters, burn = burn, update = update)
  
  cat("    Start mcmc predict \n")
  post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.i.p, knots = knots,
                               start = 1, end = iters - burn, update = update)
  timings[2, knot.design] <- fit.probit$minutes
  
  bs.pro[knot.design] <- BrierScore(post.prob.pro, y.i.p)
  post.prob.pro.med <- apply(post.prob.pro, 2, median)
  roc.pro[[knot.design]] <- roc(y.i.p ~ post.prob.pro.med)
  auc.pro[knot.design] <- roc.pro[[knot.design]]$auc
  
  print(bs.pro * 100)
  
  # copy table to tables folder on beowulf
  bs.row <- bs.row + 3
  bs[bs.row, ] <- c(bs.pro[knot.design], auc.pro[knot.design])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ####### spatial logit
  cat("  start logit \n")
  
  cat("    Start mcmc fit - Knots", knot.design, " \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  tic       <- proc.time()[3]
  fit.logit <- spGLM(formula = y.i.o ~ 1, family = "binomial",
                     coords = s.i.o, knots = knots, starting = starting,
                     tuning = tuning, priors = priors,
                     cov.model = cov.model, n.samples = iters,
                     verbose = verbose, n.report = n.report, amcmc = amcmc)
  toc        <- proc.time()[3]
  
  print("    start mcmc predict")
  yp.sp.log <- spPredict(sp.obj = fit.logit, pred.coords = s.i.p,
                         pred.covars = X.p, start = burn + 1,
                         end = iters, thin = 1, verbose = TRUE,
                         n.report = 500)
  
  post.prob.log <- t(yp.sp.log$p.y.predictive.samples)
  
  timings[3, knot.design] <- toc - tic
  
  bs.log[knot.design] <- BrierScore(post.prob.log, y.i.p)
  post.prob.log.med <- apply(post.prob.log, 2, median)
  roc.log[[knot.design]] <- roc(y.i.p ~ post.prob.log.med)
  auc.log[knot.design] <- roc.log[[knot.design]]$auc
  
  print(bs.log * 100)
  
  # copy table to tables folder on beowulf
  bs.row <- bs.row + 3
  bs[bs.row, ] <- c(bs.log[knot.design], auc.log[knot.design])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  cat("Finished: Set", i, "\n")
  save(fit.gev, bs.gev, roc.gev, auc.gev,
       fit.probit, bs.pro, roc.pro, auc.pro,
       fit.logit, bs.log, roc.log, auc.log,
       y.i.p, y.i.o, knot.design, knots,  
       s.i.o, s.i.p, timings,
       file = filename)
}
