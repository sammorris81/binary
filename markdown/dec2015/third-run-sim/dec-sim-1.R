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

source("../../../code/R/spatial_gev.R", chdir = TRUE)
source("../../../code/R/spatial_logit.R", chdir = TRUE)
source("../../../code/R/spatial_probit.R", chdir = TRUE)

# get the datasets
load("./simdata.RData")

# data setting and sets to include - written by bash script
setting <- 1

if (Sys.info()["nodename"] == "sam-ubuntu") {
  setMKLthreads(1)
  do.upload <- TRUE
} else if (Sys.info()["sysname"] == "Darwin") {
  do.upload <- TRUE
} else {
  do.upload <- FALSE
}

# extract the relevant setting from simdata
y <- simdata[[setting]]$y
s <- simdata[[setting]]$s
x <- simdata[[setting]]$x

# extract info about simulation settings
ns     <- dim(y)[1]
nt     <- 1
nknots <- 441

# testing vs training
if (setting %in% c(1, 3, 5)) {
  ntrain <- 650
} else {
  ntrain <- 1300
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

upload.pre <- "samorris@hpc.stat.ncsu.edu:~/repos-git/rare-binary/markdown/"
upload.pre <- paste(upload.pre, "dec2015/third-run-sim/sim-tables/", sep = "")

# the directory in the sim study that controls what sets are left
if (do.upload) {
  control.dir <- paste("samorris@hpc.stat.ncsu.edu:~/repos-git/rare-binary/",
                       "markdown/dec2015/third-run-sim/sim-control/", sep = "")
} else {
  control.dir <- "./sim-control/"
}

# we need this for locking while other scripts are selecting their set
lock.command <- "touch ~/repos-git/rare-binary/markdown/dec2015/third-run-sim/"
lock.command <- paste(lock.command, "sim-control/lock.txt", sep = "")
if (do.upload) {
  lock.command.ssh <- paste("ssh samorris@hpc.stat.ncsu.edu '", lock.command, "'", 
                            sep = "")
}

unlock.command <- "rm ~/repos-git/rare-binary/markdown/dec2015/third-run-sim/"
unlock.command <- paste(unlock.command, "sim-control/lock.txt", sep = "")
if (do.upload) {
  unlock.command.ssh <- paste("ssh samorris@hpc.stat.ncsu.edu '", 
                              unlock.command, "'", sep = "")
}

# download the remaining sets
if (do.upload) {
  system(paste("scp ", control.dir, "sets-remain.txt ./sim-control", sep = ""))
}
remain <- read.table(file = "./sim-control/sets-remain.txt")
if (any(remain[, setting])) {
  sets.remain <- TRUE
}

while (sets.remain) {
  # first check for the lock file
  if (do.upload) {
    system(paste("scp ", control.dir, "*.* ./sim-control", sep = ""))
  }
  
  if (!file.exists("./sim-control/lock.txt")) {
    unlocked <- TRUE
  } else {
    unlocked <- FALSE
  }
  
  if (unlocked) {  # we can actually run the sim.
    # lock
    if (do.upload) {  # we need to lock both on beowulf and locally
      system(lock.command.ssh)
    }
    system(lock.command)
    
    # load sets
    if (do.upload) {  # copy sets from beowulf
      system(paste("scp ", control.dir, "sets-remain.txt ./sim-control", 
                   sep = ""))
    }
    remain <- read.table(file = "./sim-control/sets-remain.txt")
    
    # update the remaining sets file
    if (any(remain[, setting])) {
      i <- min(which(remain[, setting]))
      remain[i, setting] <- FALSE
      write.table(remain, file = "./sim-control/sets-remain.txt")
      if (do.upload) {  # upload the remaining sets to beowulf
        system(paste("scp ./sim-control/sets-remain.txt ", control.dir, 
                     sep = ""))
      }
      
      # unlock
      system(unlock.command)
      if (do.upload) {  # need to unlock on beowulf too
        system(unlock.command.ssh)
      }
    } else {
      # unlock
      system(unlock.command)
      if (do.upload) {  # need to unlock on beowulf too
        system(unlock.command.ssh)
      }
      
      stop("No more sets to run")
    }
    
    # storage for some of the results
    scores <- matrix(NA, 3, 2)  # place to store brier scores and auc
    rownames(scores) <- c("gev", "probit", "logit")
    colnames(scores) <- c("bs", "auc")
    
    timings <- rep(NA, 3)
    
    # start the simulation
    set.seed(setting * 10 + i)
    
    # get the sites where we will be predicting
#     pred.0 <- sample(which(y[, i] == 0), size = ntest.0)
#     pred.1 <- sample(which(y[, i] == 1), size = ntest.1)
#     obs    <- rep(T, ns)
#     obs[c(pred.0, pred.1)] <- F
    # observations have already been sorted in setup.R
    obs <- rep(F, ns)
    obs[1:ntrain] <- TRUE
    
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
    # We did a simulation study to look at the impact of knot selection on 
    # the results, and it appears that a grid of knots performs similarly 
    # to randomly selecting sites for knots across all the methods.
    ####
    knots <- as.matrix(expand.grid(x = seq(0, 1, length = 21), 
                                   y = seq(0, 1, length = 21)))
    
    rho.init.pcl <- 0.05   
    dw2.o   <- rdist(s.i.o, knots)
    d.o     <- as.matrix(rdist(s.i.o))
    diag(d.o) <- 0
    
    cat("Starting: Set", i, "\n")
    
    #### spatial GEV
    cat("  Start gev \n")
    
    # using pairwise estimates as starting points for rho, alpha, and beta. also 
    # using the the pairwise estimate of alpha as the mean of the prior 
    # distribution along with a standard deviation of 0.05 to allow for some 
    # variability, but hopefully also some better convergence w.r.t. alpha.
    # we set max.dist as 0.15 in order to only consider pairs that are within 
    # a ball of radius 0.15 of one another.
    cat("    Start pairwise fit \n")
    fit.pcl <- tryCatch(
      fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                        alpha.init = 0.5, rho.init = rho.init.pcl,
                        xi.fix = TRUE, alpha.fix = FALSE,
                        rho.fix = FALSE, beta.fix = TRUE,
                        y = y.i.o, dw2 = dw2.o, d = d.o,
                        cov = X.o, method = "BFGS",
                        max.dist = 0.20,  
                        alpha.min = 0.1, alpha.max = 0.9,
                        threads = 2),
      error = function(e) {
        fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                          alpha.init = 0.5, rho.init = rho.init.pcl,
                          xi.fix = TRUE, alpha.fix = FALSE,
                          rho.fix = FALSE, beta.fix = TRUE,
                          y = y.i.o, dw2 = dw2.o, d = d.o,
                          cov = X.o, method = "Nelder-Mead",
                          max.dist = 0.20,
                          alpha.min = 0.1, alpha.max = 0.9,
                          threads = 2)
      }
    )
    cat("    Finish pairwise fit \n")
    
    cat("    Start mcmc fit \n")
    mcmc.seed <- i * 10
    set.seed(mcmc.seed)
    
    alpha.mn <- fit.pcl$par[1]
    
    # for numerical stability with the current set of starting values for the a 
    # terms. if alpha is too small, the algorithm has a very hard time getting 
    # started.
    if (alpha.mn < 0.3) {
      alpha.init <- 0.3  
    } else {
      alpha.init <- alpha.mn
    }
    
    rho.init <- fit.pcl$par[2]
    beta.init <- fit.pcl$beta
    
    fit.gev <- spatial_GEV(y = y.i.o, s = s.i.o, x = X.o, knots = knots, 
                           beta.init = log(-log(1 - mean(y.i.o))),
                           beta.mn = 0, beta.sd = 10,
                           beta.eps = 0.1, beta.attempts = 50, 
                           xi.init = 0, xi.mn = 0, xi.sd = 0.5, xi.eps = 0.01, 
                           xi.attempts = 50, xi.fix = TRUE, 
                           a.init = 10, a.eps = 0.2, a.attempts = 50, 
                           a.cutoff = 0.2, b.init = 0.5, b.eps = 0.2, 
                           b.attempts = 50, 
                           alpha.init = alpha.init, alpha.attempts = 50, 
                           alpha.mn = alpha.mn, alpha.sd = 0.05,
                           a.alpha.joint = TRUE, alpha.eps = 0.001,
                           rho.init = rho.init, logrho.mn = -2, logrho.sd = 1, 
                           rho.eps = 0.1, rho.attempts = 50, threads = 1, 
                           iters = iters, burn = burn, 
                           update = update, 
                           # update = 10, iterplot = TRUE,
                           thin = 1, thresh = 0)
    
    cat("    Start mcmc predict \n")
    post.prob.gev <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                                s.pred = s.i.p, knots = knots,
                                start = 1, end = iters - burn, update = update)
    timings[1] <- fit.gev$minutes
    
    bs.gev <- BrierScore(post.prob.gev, y.i.p)
    post.prob.gev.med <- apply(post.prob.gev, 2, median)
    roc.gev <- roc(y.i.p ~ post.prob.gev.med)
    auc.gev <- roc.gev$auc
    
    print(bs.gev * 100)
    
    # copy table to tables folder on beowulf
    scores[1, ] <- c(bs.gev, auc.gev)
    write.table(scores, file = tblname)
    if (do.upload) {
      upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
      system(upload.cmd)
    }
    
    ###### spatial probit
    cat("  Start probit \n")
    
    cat("    Start mcmc fit \n")
    mcmc.seed <- mcmc.seed + 1
    set.seed(mcmc.seed)
    fit.probit <- probit(Y = y.i.o, X = X.o, s = s.i.o, knots = knots, 
                         iters = iters, burn = burn, update = update)
    
    cat("    Start mcmc predict \n")
    post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                                 s.pred = s.i.p, knots = knots,
                                 start = 1, end = iters - burn, update = update)
    timings[2] <- fit.probit$minutes
    
    bs.pro <- BrierScore(post.prob.pro, y.i.p)
    post.prob.pro.med <- apply(post.prob.pro, 2, median)
    roc.pro <- roc(y.i.p ~ post.prob.pro.med)
    auc.pro <- roc.pro$auc
    
    print(bs.pro * 100)
    
    # copy table to tables folder on beowulf
    scores[2, ] <- c(bs.pro, auc.pro)
    write.table(scores, file = tblname)
    if (do.upload) {
      upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
      system(upload.cmd)
    }
    
    ####### spatial logit
    cat("  start logit \n")
    
    cat("    Start mcmc fit \n")
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
    
    timings[3] <- toc - tic
    
    bs.log <- BrierScore(post.prob.log, y.i.p)
    post.prob.log.med <- apply(post.prob.log, 2, median)
    roc.log <- roc(y.i.p ~ post.prob.log.med)
    auc.log <- roc.log$auc
    
    print(bs.log * 100)
    
    # copy table to tables folder on beowulf
    scores[3, ] <- c(bs.log, auc.log)
    write.table(scores, file = tblname)
    if (do.upload) {
      upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
      system(upload.cmd)
    }
    
    cat("Finished: Set", i, "\n")
    save(fit.gev, bs.gev, roc.gev, auc.gev,
         fit.probit, bs.pro, roc.pro, auc.pro,
         fit.logit, bs.log, roc.log, auc.log,
         y.i.p, y.i.o, knots,  
         s.i.o, s.i.p, timings,
         file = filename)
  } else {
    cat("Waiting for unlock... \n")
    if (do.upload) {
      system("rm ./sim-control/*.*")
    }
    Sys.sleep(10)  # wait a few seconds and try again
  }
  
  sets.remain <- any(remain[, setting])
}
