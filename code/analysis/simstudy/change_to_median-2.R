rm(list = ls())
source("./package_load.R", chdir = TRUE)

# get the datasets
load("./simdata.RData")
nsettings <- 6
nsets <- 100
setting <- 2
upload.pre <- paste("samorris@hpc.stat.ncsu.edu:~/repos-git/rare-binary/", 
                    "code/analysis/simstudy/sim-tables-med/", sep = "")

iters <- 25000; burn <- 15000; update <- 500; thin <- 1; iterplot <- FALSE
# iters <- 25000; burn <- 15000; update <- 500; thin <- 1; iterplot <- TRUE
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

# for (setting in 1:nsettings) {
  for (i in 1:nsets) {
    sets.to.change <- read.table("change_to_median.txt")
    if (sets.to.change[i, setting]) {
      sets.to.change[i, setting] <- FALSE
      write.table(x = sets.to.change, file = "change_to_median.txt")
      cat("Starting setting ", setting, ": set ", i, "\n", sep = "")
      filename <- paste("./sim-results/", setting, "-", i, ".RData", sep = "")
      if (file.exists(filename)) {
        load(filename)
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
        
        # storage for some of the results
        scores <- matrix(NA, 3, 2)  # place to store brier scores and auc
        rownames(scores) <- c("gev", "probit", "logit")
        colnames(scores) <- c("bs", "auc")
        
        timings <- rep(NA, 3)
        
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
        
        tblname  <- paste("sim-tables-med/", setting, "-", i, ".txt", sep ="")
        
        cat("    Start gev predict \n")
        post.prob.gev <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                                    s.pred = s.i.p, knots = knots,
                                    start = 1, end = iters - burn, update = update)
        
        bs.gev             <- BrierScore(post.prob.gev, y.i.p, median)
        post.prob.gev.med  <- apply(post.prob.gev, 2, median)
        post.prob.gev.mean <- apply(post.prob.gev, 2, mean)
        roc.gev            <- roc(y.i.p ~ post.prob.gev.med)
        auc.gev            <- roc.gev$auc
        
        print(bs.gev * 100)
        
        # copy table to tables folder on beowulf
        scores[1, ] <- c(bs.gev, auc.gev)
        write.table(scores, file = tblname)
        upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
        system(upload.cmd)

        # probit update
        cat("    Start pro predict \n")
        post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                                     s.pred = s.i.p, knots = knots,
                                     start = 1, end = iters - burn, 
                                     update = update)
        
        bs.pro             <- BrierScore(post.prob.pro, y.i.p, median)
        post.prob.pro.med  <- apply(post.prob.pro, 2, median)
        post.prob.pro.mean <- apply(post.prob.pro, 2, mean)
        roc.pro            <- roc(y.i.p ~ post.prob.pro.med)
        auc.pro            <- roc.pro$auc
        
        print(bs.pro * 100)
        
        # copy table to tables folder on beowulf
        scores[2, ] <- c(bs.pro, auc.pro)
        write.table(scores, file = tblname)
        upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
        system(upload.cmd)
        
        # logit
        cat("    Start log predict \n")
        yp.sp.log <- spPredict(sp.obj = fit.logit, pred.coords = s.i.p,
                               pred.covars = X.p, start = burn + 1,
                               end = iters, thin = 1, verbose = TRUE,
                               n.report = 500)
        post.prob.log <- t(yp.sp.log$p.y.predictive.samples)
        
        bs.log             <- BrierScore(post.prob.log, y.i.p, median)
        post.prob.log.med  <- apply(post.prob.log, 2, median)
        post.prob.log.mean <- apply(post.prob.log, 2, mean)
        roc.log            <- roc(y.i.p ~ post.prob.log.med)
        auc.log            <- roc.log$auc
        
        print(bs.log * 100)
        
        # copy table to tables folder on beowulf
        scores[3, ] <- c(bs.log, auc.log)
        write.table(scores, file = tblname)
        upload.cmd <- paste("scp ", tblname, " ", upload.pre, sep = "")
        system(upload.cmd)
        
        save(fit.gev, post.prob.gev.med, post.prob.gev.mean, 
             bs.gev, roc.gev, auc.gev,
             fit.probit, post.prob.pro.med, post.prob.pro.mean, 
             bs.pro, roc.pro, auc.pro,
             fit.logit, post.prob.log.med, post.prob.log.mean,
             bs.log, roc.log, auc.log,
             y.i.p, y.i.o, knots,  
             s.i.o, s.i.p, timings,
             file = filename)
        
      }
    }
  }  
# }