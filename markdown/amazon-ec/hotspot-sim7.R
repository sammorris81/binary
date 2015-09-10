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
source("../../code/R/probit.R", chdir = TRUE)

# openblas.set.num.threads(1)

# get the datasets
load("./simdata.RData")

# data setting and sets to include - written by bash script
settings <- c(1:4)
sets <- c(81:85)
nthreads <- 1

for (i in sets) {
  for (setting in settings) {
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
    iters <- 50000; burn <- 40000; update <- 1000; thin <- 1
    # iters <- 100; burn <- 50; update <- 10; thin <- 1
    # setup for spGLM
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
    
    
    filename <- paste("sim-results/", setting, "-", i, ".RData", sep = "")
    tblname  <- paste("sim-tables/", setting, "-", i, ".txt", sep ="")
    y.i.o <- matrix(y.o[, i], ntrain, 1)
    y.i.p <- matrix(y.p[, i], ntest, 1)
    print(paste("Starting: Set ", i, sep = ""))
    
    print("  start fit.pcl")
    # fit alpha and rho
    print("    start pcl fit")
    fit.pcl <- tryCatch(
      fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                        alpha.init = 0.5, rho.init = knots.h,
                        xi.fix = TRUE, alpha.fix = FALSE,
                        rho.fix = FALSE, beta.fix = TRUE,
                        y = y.i.o, dw2 = dw2.o, d = d.o,
                        cov = X.o, method = "BFGS",
                        max.dist = 3 * knots.h,
                        alpha.min = 0.1, alpha.max = 0.9,
                        threads = nthreads),
      error = function(e) {
        fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                          alpha.init = 0.5, rho.init = knots.h,
                          xi.fix = TRUE, alpha.fix = FALSE,
                          rho.fix = FALSE, beta.fix = TRUE,
                          y = y.i.o, dw2 = dw2.o, d = d.o,
                          cov = X.o, method = "Nelder-Mead",
                          max.dist = 3 * knots.h,
                          alpha.min = 0.1, alpha.max = 0.9,
                          threads = nthreads)
      }
    )
    
    # spatial GEV
    print("    start mcmc fit")
    mcmc.seed <- i * 10
    set.seed(mcmc.seed)
    
    if (fit.pcl$par[1] < 0.3) {
      alpha.init <- 0.25
    } else {
      alpha.init <- fit.pcl$par[1]
    }
    tic <- proc.time()[3]
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
                        xibeta.joint = FALSE, xi.fix = TRUE, threads = nthreads,
                        iters = iters, burn = burn, update = update, thin = 1)
    
    print("    start mcmc predict")
    post.prob.gev <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                                s.pred = s.p, knots = knots,
                                start = 1, end = iters - burn, update = update)
    toc <- proc.time()[3]
    timings[1] <- toc - tic
    
    bs.gev <- BrierScore(post.prob.gev, y.i.p)
    print(bs.gev * 100)
    
    # copy table to tables folder on beowulf
    bs <- rbind(bs.gev)
    write.table(bs, file = tblname)
    upload.cmd <- paste("scp ", tblname, " samorris@hpc.stat.ncsu.edu:~/rare-binary/markdown/sim-hotspot-1/sim-tables", sep = "")
    system(upload.cmd)
    
    # spatial probit
    print("  start probit")
    
    print("    start mcmc fit")
    mcmc.seed <- mcmc.seed + 1
    set.seed(mcmc.seed)
    tic        <- proc.time()[3]
    fit.probit <- probit(Y = y.i.o, X = X.o, s = s.o, knots = knots,
                         iters = iters, burn = burn, update = update)
    
    print("    start mcmc predict")
    post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                                 s.pred = s.p, knots = knots,
                                 start = 1, end = iters - burn, update = update)
    toc        <- proc.time()[3]
    timings[2] <- toc - tic
    
    bs.pro <- BrierScore(post.prob.pro, y.i.p)
    print(bs.pro * 100)
    
    # copy table to tables folder on beowulf
    bs <- rbind(bs.gev, bs.pro)
    write.table(bs, file = tblname)
    upload.cmd <- paste("scp ", tblname, " samorris@hpc.stat.ncsu.edu:~/rare-binary/markdown/sim-hotspot-1/sim-tables", sep = "")
    system(upload.cmd)
    
    # spatial logit
    print("  start logit")
    print("    start mcmc fit")
    mcmc.seed <- mcmc.seed + 1
    set.seed(mcmc.seed)
    tic       <- proc.time()[3]
    fit.logit <- spGLM(formula = y.i.o ~ 1, family = "binomial",
                       coords = s.o, knots = knots, starting = starting,
                       tuning = tuning, priors = priors,
                       cov.model = cov.model, n.samples = iters,
                       verbose = verbose, n.report = n.report, amcmc = amcmc)
    
    print("    start mcmc predict")
    yp.sp.log <- spPredict(sp.obj = fit.logit, pred.coords = s.p,
                           pred.covars = X.p, start = burn + 1,
                           end = iters, thin = 1, verbose = TRUE,
                           n.report = 500)
    
    post.prob.log <- t(yp.sp.log$p.y.predictive.samples)
    toc        <- proc.time()[3]
    timings[3] <- toc - tic
    
    bs.log <- BrierScore(post.prob.log, y.i.p)
    print(bs.log * 100)
    
    # copy table to tables folder on beowulf
    bs <- rbind(bs.gev, bs.pro, bs.log)
    write.table(bs, file = tblname)
    upload.cmd <- paste("scp ", tblname, " samorris@hpc.stat.ncsu.edu:~/rare-binary/markdown/sim-hotspot-1/sim-tables", sep = "")
    system(upload.cmd)
    
    print(paste("Finished: Set ", i, sep = ""))
    save(fit.pcl, fit.gev, post.prob.gev, bs.gev,
         fit.probit, post.prob.pro, bs.pro, 
         fit.logit, post.prob.log, bs.log,
         y.i.p, y.i.o, s, timings,
         file = filename)
    
  }
}