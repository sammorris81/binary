if (cluster) {
  samp.type <- "clu"
} else {
  samp.type <- "srs"
}

if (which.y == 1) {
  this.Y <- "Y1"
} else {
  this.Y <- "Y2"
}

sample.list <- paste(samp.type, ".lst.", this.Y, ".", n, sep = "")

for (set in these.sets) {
  print(paste("Start set ", set, sep = ""))
  
  these.train <- get(sample.list)[[set]]
  y <- get(this.Y)
  # if (cluster) {
  #   samp.type <- "clu"
  #   if (which.y == 1) {
  #     these.train <- clu.lst.Y1[[set]]
  #     y <- Y1
  #   } else if (which.y == 2) {
  #     these.train <- clu.lst.Y2[[set]]
  #     y <- Y2
  #   }
  # } else {
  #   samp.type <- "srs"
  #   if (which.y == 1) {
  #     these.train <- srs.lst.Y1[[set]]
  #     y <- Y1
  #   } else if (which.y == 2) {
  #     these.train <- srs.lst.Y2[[set]]
  #     y <- Y2
  #   }
  # }
  
  upload.pre <- paste("samorris@hpc.stat.ncsu.edu:~/repos-git/rare-binary/",
                      "code/analysis/swd/ss-tables/", sep = "")
  
  table.file   <- paste("./ss-tables/", samp.type, "-", which.y, "-", n, "-", 
                        set, ".txt", sep = "")
  results.file <- paste("./ss-results/", samp.type, "-", which.y, "-", n, "-", 
                        set, ".RData", sep = "")
  fit.file     <- paste("./ss-fit/", samp.type, "-", which.y, "-", n, "-", 
                        set, "-fit.RData", sep = "")
  sample.file  <- paste("./ss-sample/", samp.type, "-", which.y, "-", n, "-", 
                        set, ".txt", sep = "")
  
  y.o <- y[these.train]
  y.p <- y[-these.train]
  # this.yp <- paste("new", this.Y, sep = "")
  # y.p <- get(this.yp)
  
  # extract info about simulation settings
  ns     <- length(y.o)
  npred  <- length(y.p)
  nt     <- 1
  
  # scale sites so in [0, 1] x [0, 1] (or close)
  s.min <- apply(s, 2, min)
  s.max <- apply(s, 2, max)
  s.range <- c(diff(range(s[, 1])), diff(range(s[, 2])))
  s.scale <- s
  s.scale[, 1] <- (s[, 1] - s.min[1]) / max(s.range)
  s.scale[, 2] <- (s[, 2] - s.min[2]) / max(s.range)
  # knots[, 1] <- (knots[, 1] - s.min[1]) / max(s.range)
  # knots[, 2] <- (knots[, 2] - s.min[2]) / max(s.range)
  
  y.o <- matrix(y.o, ns, nt)
  s.o <- s.scale[these.train, ]
  X.o <- matrix(1, nrow(s.o), 1)
  y.p <- matrix(y.p, npred, nt)
  
  s.p <- s.scale[-these.train, ]
  # s.p <- new.s
  X.p <- matrix(1, nrow(s.p), 1)
  # knots <- s.o
  knots <- as.matrix(expand.grid(seq(0.05, 0.95, length = 15),
                                 seq(0.05, 0.95, length = 15)))
  knots <- rbind(knots, s.o[y.o == 1, ])
  nknots <- nrow(knots)
  
  rho.lower <- 0.05
  rho.upper <- 1
  nknots <- nrow(knots)
  
  this.save <- cbind(y.o, s.o)
  colnames(this.save) <- c("y", "s1", "s2")
  write.table(this.save, file = sample.file)
  
  ####################################################################
  #### Start MCMC setup: Most of this is used for the spBayes package
  ####################################################################
  iters <- 25000; burn <- 20000; update <- 500; thin <- 1; iterplot <- FALSE
  # iters <- 15000; burn <- 10000; update <- 500; thin <- 1; iterplot <- TRUE
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
  
  # storage for some of the results
  # scores <- matrix(NA, 3, 4)  # place to store brier scores and auc
  # scores <- array(NA, dim = c(3, 4, length(these.sets)))
  rownames(scores) <- c("gev", "probit", "logit")
  colnames(scores) <- c("bs", "auc", "bs.1", "bs.0")
  
  timings <- rep(NA, 3)
  
  # start the simulation
  set.seed(n * 10 + set)
  
  rho.init.pcl <- 0.05
  dw2.o     <- rdist(s.o, knots)^2
  d.o       <- as.matrix(rdist(s.o))
  diag(d.o) <- 0
  max.dist  <- 1
  
  #### spatial GEV
  cat("  Start gev \n")
  
  # using pairwise estimates as starting points for rho, alpha, and beta. also
  # using the the pairwise estimate of alpha as the mean of the prior
  # distribution along with a standard deviation of 0.05 to allow for some
  # variability, but hopefully also some better convergence w.r.t. alpha.
  # we set max.dist to 0.15 in order to only consider pairs of sites that
  # are relatively close to one another.
  cat("    Start pairwise fit \n")
  fit.pcl <- tryCatch(
    fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                      alpha.init = 0.5, rho.init = rho.init.pcl,
                      xi.fix = TRUE, alpha.fix = FALSE,
                      rho.fix = FALSE, beta.fix = TRUE,
                      y = y.o, dw2 = dw2.o, d = d.o,
                      cov = X.o, method = "BFGS",
                      max.dist = max.dist,
                      alpha.min = 0.1, alpha.max = 0.9,
                      threads = 2),
    error = function(e) {
      fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                        alpha.init = 0.5, rho.init = rho.init.pcl,
                        xi.fix = TRUE, alpha.fix = FALSE,
                        rho.fix = FALSE, beta.fix = TRUE,
                        y = y.o, dw2 = dw2.o, d = d.o,
                        cov = X.o, method = "Nelder-Mead",
                        max.dist = max.dist,
                        alpha.min = 0.1, alpha.max = 0.9,
                        threads = 2)
    }
  )
  
  cat("    Finish pairwise fit \n")
  
  cat("    Start mcmc fit \n")
  mcmc.seed <- 6262
  set.seed(mcmc.seed)
  
  alpha.mn <- fit.pcl$par[1]
  alpha.sd <- 0.05
  # when alpha is close to 0, the PS random effects have a much higher 
  # variance, and when it's close to 1, then the variance will decrease
  if (alpha.mn < 0.3) {
    a.eps <- 0.5
    b.eps <- 0.1
  } else if (alpha.mn < 0.85) {
    a.eps <- 0.1
    b.eps <- 0.1
  } else {
    a.eps <- 0.05
    b.eps <- 0.1
  }
  logrho.mn <- -3
  logrho.sd <- 1
  
  # for numerical stability with the current set of starting values for the a
  # terms. if alpha is too small, the algorithm has a very hard time getting
  # started.
  if (alpha.mn < 0.3) {
    alpha.init <- 0.3
  } else {
    alpha.init <- alpha.mn
  }
  
  rho.init <- max(fit.pcl$par[2], rho.lower + 0.05)
  beta.init <- fit.pcl$beta
  
  alpha.mn <- 2 / (2 + 5)
  alpha.sd <- sqrt(2 * 5 / (49 * 8))  # beta(2, 5)
  alpha.init <- 0.5
  
  fit.gev <- spatial_GEV(y = y.o, s = s.o, x = X.o, knots = knots,
                         beta.init = log(-log(1 - mean(y.o))),
                         beta.mn = 0, beta.sd = 10,
                         beta.eps = 0.1, beta.attempts = 50,
                         xi.init = 0, xi.mn = 0, xi.sd = 0.5, xi.eps = 0.01,
                         xi.attempts = 500, xi.fix = TRUE,
                         a.init = 1, a.eps = a.eps, a.attempts = 50,
                         a.cutoff = 1, a.steps = 7,
                         b.init = 0.5, b.eps = b.eps,
                         b.attempts = 500, b.steps = 5,
                         alpha.init = alpha.init, alpha.attempts = 50,
                         alpha.mn = alpha.mn, alpha.sd = alpha.sd,
                         a.alpha.joint = FALSE, alpha.eps = 0.01,
                         rho.init = rho.init, logrho.mn = logrho.mn, 
                         logrho.sd = logrho.sd, 
                         rho.lower = rho.lower, rho.upper = rho.upper,
                         rho.eps = 0.1, rho.attempts = 50, threads = 1,
                         iters = iters, burn = burn,
                         update = update, 
                         # iterplot = iterplot,
                         # update = 10, 
                         iterplot = TRUE,
                         thin = thin, thresh = 0)
  
  cat("    Start mcmc predict \n")
  y.pred.gev <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                           s.pred = s.p, knots = knots,
                           start = 1, end = iters - burn, update = update, 
                           thin = 10)
  timings[1] <- fit.gev$minutes
  
  post.prob.gev <- apply(y.pred.gev, 2, mean)
  bs.gev        <- mean((y.p - post.prob.gev)^2)
  bs.1.gev      <- mean((y.p[y.p == 1] - post.prob.gev[y.p == 1])^2)
  bs.0.gev      <- mean((y.p[y.p == 0] - post.prob.gev[y.p == 0])^2)
  roc.gev       <- roc(y.p ~ post.prob.gev)
  auc.gev       <- roc.gev$auc
  rocs.gev[[set]] <- roc.gev

  print(bs.gev * 100)
  rm(y.pred.gev)  # to help conserve memory

  # copy table to tables folder on beowulf
  scores[1, , set] <- c(bs.gev, auc.gev, bs.1.gev, bs.0.gev)
  # write.table(scores, file = table.file)
  # if (do.upload) {
  #   upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  #   system(upload.cmd)
  # }
  
  ###### spatial probit
  cat("  Start probit \n")
  
  cat("    Start mcmc fit \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  fit.probit <- probit(Y = y.o, X = X.o, s = s.o, knots = knots,
                       logbw.mn = logrho.mn, logbw.sd = logrho.sd,
                       bw.lower = rho.lower, bw.upper = rho.upper,
                       a = 1, b = 1,
                       iters = iters, burn = burn, update = update,
                       iterplot = TRUE)
  
  cat("    Start mcmc predict \n")
  y.pred.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                            s.pred = s.p, knots = knots, thin = 10,
                            start = 1, end = iters - burn, update = update)
  timings[2] <- fit.probit$minutes

  post.prob.pro <- apply(y.pred.pro, 2, mean)
  bs.pro        <- mean((y.p - post.prob.pro)^2)
  bs.1.pro      <- mean((y.p[y.p == 1] - post.prob.pro[y.p == 1])^2)
  bs.0.pro      <- mean((y.p[y.p == 0] - post.prob.pro[y.p == 0])^2)
  roc.pro       <- roc(y.p ~ post.prob.pro)
  auc.pro       <- roc.pro$auc
  rocs.pro[[set]] <- roc.pro

  print(bs.pro * 100)
  rm(y.pred.pro)  # to help conserve memory

  # copy table to tables folder on beowulf
  scores[2, , set] <- c(bs.pro, auc.pro, bs.1.pro, bs.0.pro)
  # probit 0.02791056 0.5523248: rho.upper = 1
  
  # write.table(scores, file = table.file)
  # if (do.upload) {
  #   upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  #   system(upload.cmd)
  # }
  
  ####### spatial logit
  cat("  start logit \n")
  
  cat("    Start mcmc fit \n")
  mcmc.seed <- mcmc.seed + 1
  set.seed(mcmc.seed)
  tic       <- proc.time()[3]
  fit.logit <- spGLM(formula = y.o ~ 1, family = "binomial",
                     coords = s.o, knots = knots, starting = starting,
                     tuning = tuning, priors = priors,
                     cov.model = cov.model, n.samples = iters,
                     verbose = verbose, n.report = n.report, amcmc = amcmc)
  toc        <- proc.time()[3]
  
  print("    start mcmc predict")
  post.prob.log <- spPredict(sp.obj = fit.logit, pred.coords = s.p,
                             pred.covars = X.p, start = burn + 1,
                             end = iters, thin = 1, verbose = TRUE,
                             n.report = 500)$p.y.predictive.samples

  post.prob.log <- t(post.prob.log)
  y.pred.log <- matrix(
    rbinom(n = length(post.prob.log), size = 1, prob = post.prob.log),
    nrow = nrow(post.prob.log), ncol = ncol(post.prob.log))
  rm(post.prob.log)

  timings[3] <- toc - tic

  post.prob.log <- apply(y.pred.log, 2, mean)
  bs.log        <- mean((y.p - post.prob.log)^2)
  bs.1.log      <- mean((y.p[y.p == 1] - post.prob.log[y.p == 1])^2)
  bs.0.log      <- mean((y.p[y.p == 0] - post.prob.log[y.p == 0])^2)
  roc.log       <- roc(y.p ~ post.prob.log)
  auc.log       <- roc.log$auc
  rocs.log[[set]] <- roc.log

  print(bs.log * 100)
  rm(y.pred.log)

  # copy table to tables folder on beowscoulf
  scores[3, , set] <- c(bs.log, auc.log, bs.1.log, bs.0.log)
  # write.table(scores, file = table.file)
  # if (do.upload) {
  #   upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  #   system(upload.cmd)
  # }
  
  # if ((set - 1) %% 5 == 0) {
  #   save(fit.gev, fit.probit, fit.logit, 
  #        post.prob.gev, post.prob.pro, post.prob.log,
  #        y.o, y.p, s.o, s.p, file = fit.file)
  # }
  # save(post.prob.gev, post.prob.pro, post.prob.log, 
  #      y.o, y.p, s.o, s.p, file = results.file)
}