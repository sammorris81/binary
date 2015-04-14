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
#   1: Independent logit
#   2: Independent probit
#   3: Independent GEV (Wang and Dey)
#   4: Spatial logit (spbayes)
#   5: Spatial probit
#   6: Spatial GEV (our method)
#########################################################################
# NOTE: if rerunning with covariates, make sure you adjust X matrix accordingly
rm(list=ls())
load("simdata1.RData")
source("initialize.R", chdir=T)  # loads packages and sources files

setting <- 1
analysis <- "indlogit"

# y is ns, nt, nsets, nsettings
iters <- 100000; burn <- 80000; update <- 1000; thin <- 1
nsets <- 5
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
               "phi.unif"=c(0.5, 1e4), "sigma.sq.ig"=c(1, 1),
               "tau.sq.ig"=c(1, 1))
cov.model <- "exponential"

for (g in 1:2) {
  fit.1 <- vector(mode = "list", length = nsets)
  fit.2 <- fit.3 <- fit.4 <- fit.5 <- fit.6 <- fit.1
  outputfile <- paste(setting, "-", analysis, "-", g, ".RData", sep="")

  start <- proc.time()
  for (d in 1:nsets) {
    dataset <- (g - 1) * 5 + d
    cat("start dataset", dataset, "\n")
    cur.seed <- setting * 100 + dataset  # need for MCMCpack

    y.d <- y[, , dataset, setting]

    if (nt == 1) {
      y.o <- matrix(y.d[obs], ntrain, 1)
      X.o <- matrix(X[obs], ntrain, 1)
    } else {
      y.o <- y.d[obs, , drop = F]
      X.o <- X[obs, , drop = F]
      s.o <- s[obs, ]

      y.validate[, , d] <- y.d[!obs, ]
      X.p <- X[!obs, ]
      s.p <- s[!obs, ]
    }

    # independent logit - sets seed inside MCMClogit
    fit.1[[d]] <- MCMClogit(formula = y.o ~ 1,
                            burnin = burn, mcmc = (iters - burn),
                            tune = 0.5,
                            verbose = update, seed = cur.seed,
                            beta.start = 0, B0 = 0.01)
    cur.seed <- cur.seed + 1

    # independent probit - sets seed inside MCMCprobit
    fit.2[[d]] <- MCMCprobit(formula = y.o ~ 1,
                             burnin = burn, mcmc = (iters - burn),
                             verbose = update, seed = cur.seed, B0 = 0.01)

    cur.seed <- cur.seed + 1

    # independent GEV
    set.seed(cur.seed)
    fit.3[[d]] <- mcmc(y = y.o, s = s.o, x = x.o, beta.init = 0, beta.m = 0,
                       beta.s = 100, xi.init = 0.1, xi.m = 0, xi.s = 0.5,
                       beta.tune = 0.01, xi.tune = 0.1, beta.attempts = 50,
                       xi.attempts = 50, spatial = FALSE, iterplot = FALSE,
                       iters = iters, burn = burn, update = update, thin = 1)
    cur.seed <- cur.seed + 1

    # spatial logit
    set.seed(cur.seed)
    fit.4[[d]] <- spGLM(formula = y.o ~ 1, family = "binomial", coords = s.o,
                        knots = knots, starting = starting, tuning = tuning,
                        priors = priors, cov.model = cov.model,
                        n.samples = iters, verbose = verbose,
                        n.report = n.report)
    cur.seed <- cur.seed + 1

    # spatial probit
    set.seed(cur.seed)
    fit.5[[d]] <- probit(Y = y.o, X = X.o, s = s.o, iters = iters, burn = burn,
                         update = update)
    cur.seed <- cur.seed + 1

    # spatial GEV
    set.seed(cur.seed)
    fit.6[[d]] <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                       beta.init = 0, beta.m = 0, beta.s = 100,
                       xi.init = 0.1, xi.m = 0, xi.s = 0.5,
                       knots = knots, beta.tune = 1, xi.tune = 1,
                       alpha.tune = 0.05, rho.tune = 0.05, A.tune = 1,
                       beta.attempts = 50, xi.attempts = 50,
                       alpha.attempts = 200, rho.attempts = 200,
                       spatial = TRUE, rho.init = 1, rho.upper = Inf,
                       alpha.init = 0.5, a.init = 1, iterplot = FALSE,
                       iters = iters, burn = burn, update = update, thin = 1)

    save(fit.1, fit.2, fit.3, fit.4, fit.5, fit.6, file=outputfile)
  }


}
