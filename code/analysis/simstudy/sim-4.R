rm(list=ls())
source("./package_load.R", chdir = TRUE)

# get the datasets
load("./simdata.RData")

# data setting and sets to include - written by bash script
setting <- 4

# extract the relevant setting from simdata
y <- simdata[[setting]]$y
s <- simdata[[setting]]$s
x <- simdata[[setting]]$x

# extract info about simulation settings
ns     <- dim(y)[1]
nt     <- 1

#### Knot setup
# We did a simulation study to look at the impact of knot selection on 
# the results, and it appears that a grid of knots performs similarly 
# to randomly selecting sites for knots across all the methods.
####
knots.grid <- as.matrix(expand.grid(x = seq(1 / 30, 29 / 30, length = 15), 
                                    y = seq(1 / 30, 29 / 30, length = 15)))

rho.lower <- 1 / 30 
rho.upper <- 1
nknots <- nrow(knots.grid)

# testing vs training
if (setting %in% c(1, 3, 5)) {
  ntrain <- 100
} else {
  ntrain <- 250
}

ntest   <- ns - ntrain
ntest.1 <- floor(0.05 * ntest)
ntest.0 <- ntest - ntest.1

####################################################################
#### Start MCMC setup: Most of this is used for the spBayes package
####################################################################
iters <- 25000; burn <- 20000; update <- 500; thin <- 1; iterplot <- FALSE
# iters <- 100; burn <- 20; update <- 50; thin <- 1; iterplot <- TRUE
n.report     <- 10
batch.length <- 100
n.batch      <- floor(iters / batch.length)
verbose      <- TRUE
tuning       <- list("phi" = 0.1, "sigma.sq" = 0.2, "beta" = 1, "w" = 5)
starting     <- list("phi" = 3/0.5, "sigma.sq" = 50, "beta" = 0, "w" = 0)
priors       <- list("beta.norm" = list(0, 100),
                     "phi.unif" = c(1 / rho.upper, 1 / rho.lower), 
                     "sigma.sq.ig" = c(0.1, 0.1))
cov.model <- "exponential"
timings   <- rep(NA, 3)
# with so many knots, adaptive is time prohibitive
amcmc     <- list("n.batch" = n.batch, "batch.length" = batch.length,
                  "accept.rate" = 0.35)

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
    scores <- matrix(NA, 3, 4)  # place to store brier scores and auc
    rownames(scores) <- c("gev", "probit", "logit")
    colnames(scores) <- c("bs", "auc", "bs.1", "bs.0")
    
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
    
    knots <- rbind(knots.grid, s.i.o[y.i.o == 1, ])
    
    ntrain.0 <- sum(y.i.o == 0) - ntest.0
    ntrain.1 <- sum(y.i.o == 1) - ntest.1
    
    table.file   <- paste("./sim-tables/", setting, "-", i, ".txt", sep ="")
    results.file <- paste("./sim-results/", setting, "-", i, ".RData", sep = "")
    fit.file     <- paste("./sim-fit/", setting, "-", i, "-fit.RData", sep = "")
    
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
    # cat("    Start pairwise fit \n")
    # fit.pcl <- tryCatch(
    #   fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
    #                     alpha.init = 0.5, rho.init = rho.init.pcl,
    #                     xi.fix = TRUE, alpha.fix = FALSE,
    #                     rho.fix = FALSE, beta.fix = TRUE,
    #                     y = y.i.o, dw2 = dw2.o, d = d.o,
    #                     cov = X.o, method = "BFGS",
    #                     max.dist = 0.20,  
    #                     alpha.min = 0.1, alpha.max = 0.9,
    #                     threads = 2),
    #   error = function(e) {
    #     fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
    #                       alpha.init = 0.5, rho.init = rho.init.pcl,
    #                       xi.fix = TRUE, alpha.fix = FALSE,
    #                       rho.fix = FALSE, beta.fix = TRUE,
    #                       y = y.i.o, dw2 = dw2.o, d = d.o,
    #                       cov = X.o, method = "Nelder-Mead",
    #                       max.dist = 0.20,
    #                       alpha.min = 0.1, alpha.max = 0.9,
    #                       threads = 2)
    #   }
    # )
    # cat("    Finish pairwise fit \n")
    
    cat("    Start mcmc fit \n")
    mcmc.seed <- i * 10
    set.seed(mcmc.seed)
    
    # alpha.mn <- fit.pcl$par[1]
    
    # for numerical stability with the current set of starting values for the a 
    # terms. if alpha is too small, the algorithm has a very hard time getting 
    # started.
    # if (alpha.mn < 0.3) {
    #   alpha.init <- 0.3  
    # } else {
    #   alpha.init <- alpha.mn
    # }
    # 
    # rho.init <- fit.pcl$par[2]
    # beta.init <- fit.pcl$beta
    
    rho.init <- (rho.upper - rho.lower) / 2
    
    alpha.mn <- 2 / (2 + 5)
    alpha.sd <- sqrt(2 * 5 / (49 * 8))  # beta(2, 5)
    alpha.init <- 0.3
    
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
                           alpha.mn = alpha.mn, alpha.sd = alpha.sd,
                           a.alpha.joint = FALSE, alpha.eps = 0.01,
                           rho.init = rho.init, 
                           rho.lower = rho.lower, rho.upper = rho.upper,
                           rho.eps = 0.1, rho.attempts = 50, threads = 1, 
                           iters = iters, burn = burn, 
                           update = update, iterplot = iterplot,
                           # update = 10, iterplot = TRUE,
                           thin = 1, thresh = 0)
    
    cat("    Start mcmc predict \n")
    y.pred.gev <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                             s.pred = s.i.p, knots = knots, thin = 10,
                             start = 1, end = iters - burn, update = update)
    timings[1] <- fit.gev$minutes
    
    post.prob.gev <- apply(y.pred.gev, 2, mean)
    bs.gev        <- mean((y.i.p - post.prob.gev)^2)
    bs.1.gev      <- mean((y.i.p[y.i.p == 1] - post.prob.gev[y.i.p == 1])^2)
    bs.0.gev      <- mean((y.i.p[y.i.p == 0] - post.prob.gev[y.i.p == 0])^2)
    roc.gev       <- roc(y.i.p ~ post.prob.gev)
    auc.gev       <- roc.gev$auc
    
    print(bs.gev * 100)
    rm(y.pred.gev)  # to help conserve memory
    
    # copy table to tables folder on beowulf
    scores[1, ] <- c(bs.gev, auc.gev, bs.1.gev, bs.0.gev)
    write.table(scores, file = table.file)
    if (do.upload) {
      upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
      system(upload.cmd)
    }
    
    ###### spatial probit
    cat("  Start probit \n")
    
    cat("    Start mcmc fit \n")
    mcmc.seed <- mcmc.seed + 1
    set.seed(mcmc.seed)
    fit.probit <- probit(Y = y.i.o, X = X.o, s = s.i.o, knots = knots, 
                         bw.lower = rho.lower, bw.upper = rho.upper,
                         a = 0.1, b = 0.1,
                         iters = iters, burn = burn, update = update,
                         iterplot = iterplot)
    
    cat("    Start mcmc predict \n")
    y.pred.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                              s.pred = s.i.p, knots = knots, thin = 10,
                              start = 1, end = iters - burn, update = update)
    timings[2] <- fit.probit$minutes
    
    post.prob.pro <- apply(y.pred.pro, 2, mean)
    bs.pro        <- mean((y.i.p - post.prob.pro)^2)
    bs.1.pro      <- mean((y.i.p[y.i.p == 1] - post.prob.pro[y.i.p == 1])^2)
    bs.0.pro      <- mean((y.i.p[y.i.p == 0] - post.prob.pro[y.i.p == 0])^2)
    roc.pro       <- roc(y.i.p ~ post.prob.pro)
    auc.pro       <- roc.pro$auc
    
    print(bs.pro * 100)
    
    # copy table to tables folder on beowulf
    scores[2, ] <- c(bs.pro, auc.pro, bs.1.pro, bs.0.pro)
    write.table(scores, file = table.file)
    if (do.upload) {
      upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
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
    post.prob.log <- spPredict(sp.obj = fit.logit, pred.coords = s.i.p,
                               pred.covars = X.p, start = burn + 1,
                               end = iters, thin = 10, verbose = TRUE,
                               n.report = 500)$p.y.predictive.samples
    
    post.prob.log <- t(post.prob.log)
    y.pred.log <- matrix(
      rbinom(n = length(post.prob.log), size = 1, prob = post.prob.log),
      nrow = nrow(post.prob.log), ncol = ncol(post.prob.log))
    rm(post.prob.log)
    
    timings[3] <- (toc - tic) / 60
    
    post.prob.log <- apply(y.pred.log, 2, mean)
    bs.log        <- mean((y.i.p - post.prob.log)^2)
    bs.1.log      <- mean((y.i.p[y.i.p == 1] - post.prob.log[y.i.p == 1])^2)
    bs.0.log      <- mean((y.i.p[y.i.p == 0] - post.prob.log[y.i.p == 0])^2)
    roc.log       <- roc(y.i.p ~ post.prob.log)
    auc.log       <- roc.log$auc
    
    print(bs.log * 100)
    rm(y.pred.log)
    
    # copy table to tables folder on beowulf
    scores[3, ] <- c(bs.log, auc.log, bs.1.log, bs.0.log)
    write.table(scores, file = table.file)
    if (do.upload) {
      upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
      system(upload.cmd)
    }
    
    cat("Finished: Set", i, "\n")
    save(fit.gev, fit.probit, fit.logit,
         post.prob.gev, post.prob.pro, post.prob.log,
         y.i.p, y.i.o, s.i.o, s.i.p, knots, timings, file = fit.file)
    
    save(post.prob.gev, post.prob.pro, post.prob.log,
         y.i.o, y.i.p, s.i.o, s.i.p, file = results.file)
  } else {
    cat("Waiting for unlock... \n")
    if (do.upload) {
      system("rm ./sim-control/*.*")
    }
    Sys.sleep(10)  # wait a few seconds and try again
  }
  
  sets.remain <- any(remain[, setting])
}