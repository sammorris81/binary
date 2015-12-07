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
setMKLthreads(1)
setting <- 2
sets    <- 6:10 

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

timings <- matrix(NA, 3, 3)
rownames(timings) <- c("gev", "probit", "logit")
colnames(timings) <- c("knots 1", "knots 2", "knots 3")

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
  
  bs <- matrix(NA, 9, 2)  # place to store brier scores and auc
  rownames(bs) <- c("gev-1", "gev-2", "gev-3", 
                    "probit-1", "probit-2", "probit-3",
                    "logit-1", "logit-2", "logit-3")
  colnames(bs) <- c("bs", "auc")
  bs.gev  <- bs.pro  <- bs.log  <- rep(NA, 3)
  auc.gev <- auc.pro <- auc.log <- rep(NA, 3)
  roc.gev <- roc.pro <- roc.log <- vector(mode = "list", length = 3)
  
  #### Knot setup
  # for ns.o = 500, using all sites
  # 
  # for ns.o = 1000, using a space filling design with 500 knots
  # There is a bug in cover.design when you have set nn = FALSE. So, to avoid 
  # the warning, setting number of nearest neighbors to 0 and turning off 
  # nearest neighbors. This will take a bit longer, but ultimately comes back 
  # with a similar design
  ####
  set.seed(setting * 10 + i)
  knots.1 <- as.matrix(expand.grid(x = seq(0, 1, length = 21), 
                                   y = seq(0, 1, length = 21)))
  knots.2 <- cover.design(R = s.i.o, nd = nknots, nruns = 1, nn = FALSE, 
                          num.nn = 0)$design
  # stratified sample
  phat.o   <- mean(y.i.o)  # what proportion of the sites in training are 1s
  nknots.1 <- floor(phat.o * nknots)
  nknots.0 <- floor((1 - phat.o) * nknots)
  knots.3.0 <- cover.design(R = s.i.o[y.i.o == 0, ], nd = nknots.0, nruns = 1, 
                            nn = FALSE, num.nn = 0)$design
  knots.3.1 <- cover.design(R = s.i.o[y.i.o == 1, ], nd = nknots.1, nruns = 1, 
                            nn = FALSE, num.nn = 0)$design
  knots.3 <- rbind(knots.3.0, knots.3.1)
  
  cat("Starting: Set", i, "\n")
  
  #### spatial GEV
  cat("  Start gev \n")
  
  ## Knot setup 1
  cat("    Start mcmc fit - Knots 1 \n")
  mcmc.seed <- i * 10
  set.seed(mcmc.seed)
  
  fit.gev <- spatial_GEV(y = y.i.o, s = s.i.o, x = X.o, knots = knots.1, 
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
                              s.pred = s.i.p, knots = knots.1,
                              start = 1, end = iters - burn, update = update)
  timings[1, 1] <- fit.gev$minutes
  
  bs.gev[1] <- BrierScore(post.prob.gev, y.i.p)
  post.prob.gev.med <- apply(post.prob.gev, 2, median)
  roc.gev[[1]] <- roc(y.i.p ~ post.prob.gev.med)
  auc.gev[1] <- roc.gev[[1]]$auc
  
  print(bs.gev[1] * 100)
  
  # copy table to tables folder on beowulf
  bs[1, ] <- c(bs.gev[1], auc.gev[1])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ## Knot setup 2
  cat("    Start mcmc fit - Knots 2 \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  
  fit.gev <- spatial_GEV(y = y.i.o, s = s.i.o, x = X.o, knots = knots.2, 
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
                              s.pred = s.i.p, knots = knots.2,
                              start = 1, end = iters - burn, update = update)
  timings[1, 2] <- fit.gev$minutes
  
  bs.gev[2] <- BrierScore(post.prob.gev, y.i.p)
  post.prob.gev.med <- apply(post.prob.gev, 2, median)
  roc.gev[[2]] <- roc(y.i.p ~ post.prob.gev.med)
  auc.gev[2] <- roc.gev[[2]]$auc
  
  print(bs.gev * 100)
  
  # copy table to tables folder on beowulf
  bs[2, ] <- c(bs.gev[2], auc.gev[2])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ## Knot setup 3
  cat("    Start mcmc fit - Knots 3 \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  
  fit.gev <- spatial_GEV(y = y.i.o, s = s.i.o, x = X.o, knots = knots.3, 
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
                              s.pred = s.i.p, knots = knots.3,
                              start = 1, end = iters - burn, update = update)
  timings[1, 3] <- fit.gev$minutes
  
  bs.gev[3] <- BrierScore(post.prob.gev, y.i.p)
  post.prob.gev.med <- apply(post.prob.gev, 2, median)
  roc.gev[[3]] <- roc(y.i.p ~ post.prob.gev.med)
  auc.gev[3] <- roc.gev[[3]]$auc
  
  print(bs.gev * 100)
  
  # copy table to tables folder on beowulf
  bs[3, ] <- c(bs.gev[3], auc.gev[3])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ###### spatial probit
  cat("  Start probit \n")
  
  ## Knot setup 1
  cat("    Start mcmc fit - Knots 1 \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  fit.probit <- probit(Y = y.i.o, X = X.o, s = s.i.o, knots = knots.1, 
                       iters = iters, burn = burn, update = update)
  
  cat("    Start mcmc predict \n")
  post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.i.p, knots = knots.1,
                               start = 1, end = iters - burn, update = update)
  timings[2, 1] <- fit.probit$minutes
  
  bs.pro[1] <- BrierScore(post.prob.pro, y.i.p)
  post.prob.pro.med <- apply(post.prob.pro, 2, median)
  roc.pro[[1]] <- roc(y.i.p ~ post.prob.pro.med)
  auc.pro[1] <- roc.pro[[1]]$auc
  
  print(bs.pro * 100)
  
  # copy table to tables folder on beowulf
  bs[4, ] <- c(bs.pro[1], auc.pro[1])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ## Knot setup 2
  cat("    Start mcmc fit - Knots 2 \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  fit.probit <- probit(Y = y.i.o, X = X.o, s = s.i.o, knots = knots.2, 
                       iters = iters, burn = burn, update = update)
  
  cat("    Start mcmc predict \n")
  post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.i.p, knots = knots.2,
                               start = 1, end = iters - burn, update = update)
  timings[2, 2] <- fit.probit$minutes
  
  bs.pro[2] <- BrierScore(post.prob.pro, y.i.p)
  post.prob.pro.med <- apply(post.prob.pro, 2, median)
  roc.pro[[2]] <- roc(y.i.p ~ post.prob.pro.med)
  auc.pro[2] <- roc.pro[[2]]$auc
  
  print(bs.pro * 100)
  
  # copy table to tables folder on beowulf
  bs[5, ] <- c(bs.pro[2], auc.pro[2])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ## Knot setup 3
  cat("    Start mcmc fit - Knots 3 \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  fit.probit <- probit(Y = y.i.o, X = X.o, s = s.i.o, knots = knots.3, 
                       iters = iters, burn = burn, update = update)
  
  cat("    Start mcmc predict \n")
  post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.i.p, knots = knots.3,
                               start = 1, end = iters - burn, update = update)
  timings[2, 3] <- fit.probit$minutes
  
  bs.pro[3] <- BrierScore(post.prob.pro, y.i.p)
  post.prob.pro.med <- apply(post.prob.pro, 2, median)
  roc.pro[[3]] <- roc(y.i.p ~ post.prob.pro.med)
  auc.pro[3] <- roc.pro[[3]]$auc
  
  print(bs.pro * 100)
  
  # copy table to tables folder on beowulf
  bs[6, ] <- c(bs.pro[3], auc.pro[3])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ####### spatial logit
  cat("  start logit \n")
  
  ## Knot setup 1
  cat("    Start mcmc fit - Knots 1 \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  tic       <- proc.time()[3]
  fit.logit <- spGLM(formula = y.i.o ~ 1, family = "binomial",
                     coords = s.i.o, knots = knots.1, starting = starting,
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
  
  timings[3, 1] <- toc - tic
  
  bs.log[1] <- BrierScore(post.prob.log, y.i.p)
  post.prob.log.med <- apply(post.prob.log, 2, median)
  roc.log[[1]] <- roc(y.i.p ~ post.prob.log.med)
  auc.log[1] <- roc.log[[1]]$auc
  
  print(bs.log * 100)
  
  # copy table to tables folder on beowulf
  bs[7, ] <- c(bs.log[1], auc.log[1])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ## Knot setup 2
  cat("    Start mcmc fit - Knots 2 \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  tic       <- proc.time()[3]
  fit.logit <- spGLM(formula = y.i.o ~ 1, family = "binomial",
                     coords = s.i.o, knots = knots.2, starting = starting,
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
  
  timings[3, 2] <- toc - tic
  
  bs.log[2] <- BrierScore(post.prob.log, y.i.p)
  post.prob.log.med <- apply(post.prob.log, 2, median)
  roc.log[[2]] <- roc(y.i.p ~ post.prob.log.med)
  auc.log[2] <- roc.log[[2]]$auc
  
  print(bs.log * 100)
  
  # copy table to tables folder on beowulf
  bs[8, ] <- c(bs.log[2], auc.log[2])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  ## Knot setup 3
  cat("    Start mcmc fit - Knots 3 \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  tic       <- proc.time()[3]
  fit.logit <- spGLM(formula = y.i.o ~ 1, family = "binomial",
                     coords = s.i.o, knots = knots.3, starting = starting,
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
  
  timings[3, 3] <- toc - tic
  
  bs.log[3] <- BrierScore(post.prob.log, y.i.p)
  post.prob.log.med <- apply(post.prob.log, 2, median)
  roc.log[[3]] <- roc(y.i.p ~ post.prob.log.med)
  auc.log[3] <- roc.log[[3]]$auc
  
  print(bs.log * 100)
  
  # copy table to tables folder on beowulf
  bs[9, ] <- c(bs.log[3], auc.log[3])
  write.table(bs, file = tblname)
  if (do.upload) {
    upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
    system(upload.cmd)
  }
  
  cat("Finished: Set", i, "\n")
  save(bs.gev, roc.gev, auc.gev,
       bs.pro, roc.pro, auc.pro,
       bs.log, roc.log, auc.log,
       y.i.p, y.i.o, knots.1, knots.2, knots.3, 
       s.i.o, s.i.p, timings,
       file = filename)
}
