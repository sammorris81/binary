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
load("simdata1.RData")
source("initialize.R", chdir=T)  # loads packages and sources files

setting <- 8

# y is ns, nt, nsets, nsettings
iters <- 100000; burn <- 80000; update <- 500; thin <- 1
nsets <- 10
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

start <- proc.time()

# (1 - 8)-(1 - 5)-(1 - 10)
# for fit method 1 (*-1-*: independent probit), the posterior P(Y = 1) = 1 = Phi(fit)
# for fit method 2 (*-2-*: independent GEV), the posterior
  z = rep(NA, length(fit$xi))
  for (i in 1:length(fit$xi)) {
    z[i] = getZ(xi = fit$xi[i], x.beta = fit$beta[i], thresh = 0)
  }
#   
#   P(Y = 1) = exp(-1 / z)
for (d in 1:nsets) {
  cat("start dataset", d, "\n")

  y.d <- y[, , d, setting]

  if (nt == 1) {
    y.o <- matrix(y.d[obs], ntrain, 1)
    X.o <- matrix(X[obs], ntrain, 1)
  } else {
    y.o <- y.d[obs, , drop = F]
    X.o <- X[obs, , drop = F]
  }
  s.o <- s[obs, ]

  # independent probit - sets seed inside MCMCprobit
  cur.seed <- setting * 100 + d  # need for MCMCpack
  cat("Start independent probit: Set", d, "\n")
  fit <- MCMCprobit(formula = y.o ~ 1,
                    burnin = burn, mcmc = (iters - burn),
                    verbose = update, seed = cur.seed, B0 = 0.01)
  cat("End independent probit: Set", d, "\n")
  # setting-method-set.RData
  outputfile <- paste(setting, "-1-", d, ".RData", sep="")
  save(fit, file = outputfile)
  rm(fit)
  gc()


  # independent GEV
  cur.seed <- cur.seed + 1
  set.seed(cur.seed)
  cat("Start independent GEV: Set", d, "\n")
  fit <- mcmc(y = y.o, s = s.o, x = X.o, beta.init = 0, beta.m = 0,
              beta.s = 100, xi.init = 0.1, xi.m = 0, xi.s = 0.5,
              beta.tune = 0.01, xi.tune = 0.1, beta.attempts = 50,
              xi.attempts = 50, spatial = FALSE, iterplot = FALSE,
              iters = iters, burn = burn, update = update, thin = 1)
  cat("End independent GEV: Set", d, "\n")
  outputfile <- paste(setting, "-2-", d, ".RData", sep="")
  save(fit, file = outputfile)
  rm(fit)
  gc()

  # spatial logit
  cur.seed <- cur.seed + 1
  set.seed(cur.seed)
  cat("Start spatial logit: Set", d, "\n")
  fit <- spGLM(formula = y.o ~ 1, family = "binomial", coords = s.o,
               knots = knots, starting = starting, tuning = tuning,
               priors = priors, cov.model = cov.model,
               n.samples = iters, verbose = verbose,
               n.report = n.report)
  cat("End spatial logit: Set", d, "\n")
  outputfile <- paste(setting, "-3-", d, ".RData", sep="")
  save(fit, file = outputfile)
  rm(fit)
  gc()

  # spatial probit
  cur.seed <- cur.seed + 1
  set.seed(cur.seed)
  cat("Start spatial probit: Set", d, "\n")
  fit <- probit(Y = y.o, X = X.o, s = s.o, knots = knots,
                iters = iters, burn = burn, update = update)
  cat("End spatial probit: Set", d, "\n")
  outputfile <- paste(setting, "-4-", d, ".RData", sep="")
  save(fit, file = outputfile)
  rm(fit)
  gc()


  # spatial GEV
  cur.seed <- cur.seed + 1
  set.seed(cur.seed)
  cat("Start spatial GEV: Set", d, "\n")
  fit <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
              beta.init = 0, beta.m = 0, beta.s = 100,
              xi.init = 0.1, xi.m = 0, xi.s = 0.5,
              knots = knots, beta.tune = 1, xi.tune = 1,
              alpha.tune = 0.05, rho.tune = 0.05, A.tune = 1,
              beta.attempts = 50, xi.attempts = 50,
              alpha.attempts = 200, rho.attempts = 200,
              spatial = TRUE, rho.init = 1, rho.upper = 9,
              alpha.init = 0.5, a.init = 1, iterplot = FALSE,
              iters = iters, burn = burn, update = update, thin = 1)
  cat("End spatial GEV: Set", d, "\n")
  outputfile <- paste(setting, "-5-", d, ".RData", sep="")
  save(fit, file = outputfile)
  rm(fit)
  gc()
}
