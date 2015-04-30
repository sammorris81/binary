#########################################################################
# A small-scale simulation study to determine when spatial GEV link
# performs better than
#
# data settings:
#   All: s in [0, 6] x [0, 6]
#   1: GEV link
#      a: alpha = 0.3, 100 knots, 1% rareness, bw = 3
#      b: alpha = 0.7, 100 knots, 1% rareness, bw = 3
#      c: alpha = 0.3, 100 knots, 5% rareness, bw = 3
#      d: alpha = 0.7, 100 knots, 5% rareness, bw = 3
#   2: Logit link
#      a: rho = 1, 100 knots, 1% rareness
#      b: rho = 3, 100 knots, 1% rareness
#      c: rho = 1, 100 knots, 5% rareness
#      d: rho = 3, 100 knots, 5% rareness
#   3: Probit link -- Hold off for now
#      a: gamma = 0.9, 100 knots, 1% rareness
#      b: gamma = 0.1, 100 knots, 1% rareness
#      c: gamma = 0.9, 100 knots, 5% rareness
#      d: gamma = 0.1, 100 knots, 5% rareness
#
# methods:
#   1: Independent probit
#   2: Independent GEV (Wang and Dey)
#   3: Spatial logit (spbayes)
#   4: Spatial probit
#   5: Spatial GEV (our method)
#   0: Independent logit - tried and MCMClogit gets stuck
#########################################################################
# NOTE: if rerunning with covariates, make sure to adjust X matrix accordingly
rm(list=ls())
load("simdata2.RData")
source("initialize.R", chdir=T)  # loads packages and sources files
library(doMC)
registerDoMC(6)

# y is ns, nt, nsets, nsettings
iters <- 40000; burn <- 30000; update <- 500; thin <- 1
nsets <- 10
nsettings <- 8
nmethods  <- 4
ntrain <- 3000
ntest  <- 1000
obs <- c(rep(T, ntrain), rep(F, ntest))

# setup for spGLM
n.report <- 500
verbose <- TRUE
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1,
               "beta"=0.1, "w"=0.1)
starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1,
                 "beta"=0, "w"=0)
priors <- list("beta.norm"=list(1, 100),
               "phi.unif"=c(0.1, 1e4), "sigma.sq.ig"=c(1, 1),
               "tau.sq.ig"=c(1, 1))
cov.model <- "exponential"

s.pred <- s[!obs, ]
X.pred <- matrix(1, ntest, nt)

results2 <- foreach (setting = 1:nsettings) %dopar% {
  brier.scores <- matrix(NA, 10, 4)
  for (set in 1:nsets) {
    y.validate   <- y[!obs, , set, setting]

    # 1: Spatial logit
    filename <- paste(setting, "-1-", set, ".RData", sep="")
    load(filename)
    yp.sp.lo <- spPredict(sp.obj = fit, pred.coords = s.pred,
                          pred.covars = X.pred, start = 30001, end = 40000,
                          thin = 1, verbose = TRUE, n.report = 500)
    post.prob <- t(yp.sp.lo$p.y.predictive.samples)
    brier.scores[set, 1] <- BrierScore(post.prob, y.validate)
    # print(paste("Spatial logit - set", set, "finished"))

    # 2: Spatial probit
    filename <- paste(setting, "-2-", set, ".RData", sep="")
    load(filename)
    post.prob <- pred.spprob(mcmcoutput = fit, X.pred = X.pred,
                             s.pred = s.pred, knots = knots,
                             start = 1, end = 10000, update = 500)
    brier.scores[set, 2] <- BrierScore(post.prob, y.validate)
    # print(paste("Spatial probit - set", set, "finished"))

    # 3: Spatial GEV
    filename <- paste(setting, "-3-", set, ".RData", sep="")
    load(filename)
    post.prob <- pred.spgev(mcmcoutput = fit, x.pred = X.pred,
                            s.pred = s.pred, knots = knots,
                            start = 1, end = 10000, update = 500)
    brier.scores[set, 3] <- BrierScore(post.prob, y.validate)
    # print(paste("Spatial GEV - set", set, "finished"))

    if (setting <= 4) {
      # 4: Spatial GEV - Fixed alpha and rho
      filename <- paste(setting, "-4-", set, ".RData", sep="")
      load(filename)
      post.prob <- pred.spgev(mcmcoutput = fit, x.pred = X.pred,
                              s.pred = s.pred, knots = knots,
                              start = 1, end = 10000, update = 500)
      brier.scores[set, 4] <- BrierScore(post.prob, y.validate)
      # print(paste("Spatial GEV - set", set, ", fixed alpha & rho, finished"))
    }
  }
  return(brier.scores)

#   save(brier.scores, file="binary-results.RData")
#   print(paste("Setting", setting, "finished"))
}

save(results, file="binary-results.RData")  # saves a list
load("binary-results.RData")
brier.scores <- array(0, dim=c(nsets, nmethods, nsettings))

for (k in 1:nsettings) {
  brier.scores[, , k] <- results[[k]]
}

apply(brier.scores, c(2, 3), mean, na.rm=T) * 1000


