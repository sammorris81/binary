rm(list=ls())

# libraries
library(fields)
library(Rcpp)
library(evd)
library(spBayes)
library(mvtnorm)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp(file = "../code/R/pairwise.cpp")

source("../code/R/auxfunctions.R", chdir = TRUE)
source("../code/R/updateModel.R")
source("../code/R/mcmc.R")
source("../code/R/probit.R", chdir=T)

# true knots
knots.t <- as.matrix(expand.grid(seq(0.00, 1.00, length=12),
                                 seq(0.00, 1.00, length=12)))

# knots used to fit model
knots <- as.matrix(expand.grid(seq(0.00, 1.00, length=15),
                               seq(0.00, 1.00, length=15)))

# grid for searching over alpha and rho
rhos   <- (knots[2, 1] - knots[1, 1]) * seq(1, 3, by = 0.5)

# get sites and distances from knots
set.seed(7483)  # sites
ns  <- 1500
s   <- cbind(runif(ns), runif(ns))
dw2 <- rdist(s, knots)
d   <- as.matrix(rdist(s))
diag(d) <- 0
x   <- matrix(1, ns, 1)

# only focusing on the case of strong dependence
alpha.t <- 0.20
xi.t    <- 0.25
rho.t   <- 0.15
prop    <- c(0.05, 0.01)

# y is ns, nt, nsets, nsettings
iters <- 15000; burn <- 10000; update <- 500; thin <- 1
# iters <- 100; burn <- 50; update <- 10; thin <- 1
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


data.seed <- 3282
set.seed(data.seed)  # data
data <- rRareBinarySpat(x, s = s, knots = knots.t, beta = 0, xi = xi.t,
                        alpha = alpha.t, rho = rho.t, prob.success = prop[1])
y <- data$y

data.seed <- data.seed + 1
set.seed(data.seed)  # data
data <- rRareBinarySpat(x, s = s, knots = knots.t, beta = 0, xi = xi.t,
                        alpha = alpha.t, rho = rho.t, prob.success = prop[2])
y <- data$y

ntrain <- floor(0.75 * ns)
ntest  <- ns - ntrain
obs <- c(rep(T, ntrain), rep(F, ntest))
y.o <- y[obs, , drop = FALSE]
X.o <- matrix(x[obs], ntrain, 1)
s.o <- s[obs, ]
y.validate <- y[!obs, , drop = FALSE]
X.p <- matrix(x[!obs, ], ntest, 1)
s.p <- s[!obs, ]

nll <- rep(9999999, length(rhos))
for (i in seq_along(rhos)) {
  temp <- fit.rarebinaryCPP(c(0.5, 0, -4), rho = rhos[i],  y = y.o,
                            dw2 = dw2, d = d, cov = X.o, threads = 6)
  nll[i] <- temp$value
  if (nll[i] == min(nll)) {
    fit <- temp
  }
}

# W <- stdW(makeW(dw2 = dw2, rho = rho.t))
# fit <- fit.rarebinaryCPP(c(0.5, 0, -4), rho=rho.t, y = y.o, dw2 = dw2, d = d, cov=X.o, threads=6)
# pairwise.rarebinary3CPP(par = c(alpha.t, xi.t, -data$thresh), rho = rho.t, 
#                         y = y.o, d = d, max.dist = max(d), W = W, 
#                         cov = X.o[, 1], threads = 1)
# pairwise.rarebinary3CPP(par = c(0.1, xi.t, -data$thresh), rho = rho.t, 
#                         y = y.o, d = d, max.dist = max(d), W = W, 
#                         cov = X.o[, 1], threads = 1)

xibeta.hat <- fit$par[2:3]
xibeta.var <- solve(fit$hessian)[2:3, 2:3]

plot(knots.t, ylim = c(0, 1), xlim = c(0, 1), 
     main = "simulated dataset", xlab="", ylab="")
train1.idx <- which(y[obs] == 1)
test1.idx <- which(y[!obs] == 1) + ntrain  # to get to the testing
points(s[train1.idx, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[test1.idx, ], pch = 21, col = "firebrick4", bg = "firebrick1")

# spatial GEV
mcmc.seed <- 1
set.seed(mcmc.seed)
fit.gev.1 <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[3], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[2], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.01, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 200, rho.attempts = 200,
                  spatial = TRUE, rho.init = fit$rho, rho.upper = 9,
                  alpha.init = 0.3, a.init = 1000, iterplot = TRUE,
                  alpha.fix = TRUE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

fit.gev.2 <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[3], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[2], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.01, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 200, rho.attempts = 200,
                  spatial = TRUE, rho.init = fit$rho, rho.upper = 9,
                  alpha.init = 0.5, a.init = 100, iterplot = TRUE,
                  alpha.fix = TRUE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

fit.gev.3 <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[3], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[2], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.01, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 200, rho.attempts = 200,
                  spatial = TRUE, rho.init = fit$rho, rho.upper = 9,
                  alpha.init = 0.7, a.init = 100, iterplot = TRUE,
                  alpha.fix = TRUE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

fit.gev.4 <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[3], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[2], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.01, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 200, rho.attempts = 200,
                  spatial = TRUE, rho.init = fit$rho, rho.upper = 9,
                  alpha.init = 0.2, a.init = 30000, iterplot = TRUE,
                  alpha.fix = TRUE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

post.prob.gev.1 <- pred.spgev(mcmcoutput = fit.gev.1, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

post.prob.gev.2 <- pred.spgev(mcmcoutput = fit.gev.2, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

post.prob.gev.3 <- pred.spgev(mcmcoutput = fit.gev.3, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

post.prob.gev.4 <- pred.spgev(mcmcoutput = fit.gev.4, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

# spatial logit
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.logit <- spGLM(formula = y.o ~ 1, family = "binomial", coords = s.o,
                   knots = knots, starting = starting, tuning = tuning,
                   priors = priors, cov.model = cov.model,
                   n.samples = iters, verbose = verbose,
                   n.report = n.report)

yp.sp.log <- spPredict(sp.obj = fit.logit, pred.coords = s.p,
                       pred.covars = X.p, start = burn + 1, end = iters,
                       thin = 1, verbose = TRUE, n.report = 500)

post.prob.log <- t(yp.sp.log$p.y.predictive.samples)

# spatial probit
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.probit <- probit(Y = y.o, X = X.o, s = s.o, knots = knots,
                     iters = iters, burn = burn, update = update)

post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                             s.pred = s.p, knots = knots,
                             start = 1, end = iters - burn, update = 500)


bs.gev.1.1 <- BrierScore(post.prob.gev.1, y.validate)  # 0.0210
bs.gev.2.1 <- BrierScore(post.prob.gev.2, y.validate)  # 0.0292
bs.gev.3.1 <- BrierScore(post.prob.gev.3, y.validate)  # 0.0315
bs.gev.4.1 <- BrierScore(post.prob.gev.4, y.validate)  # 0.0560
bs.log.1   <- BrierScore(post.prob.log, y.validate)    # 0.0602
bs.pro.1   <- BrierScore(post.prob.pro, y.validate)    # 0.0082

# a couple of thoughts:
# In the smaller dataset, the MCMC really seems to find the posterior 
# distribution very shortly after starting.
# This is partly because we've really helped stabilize the MCMC by using a 
# candidate for xi and beta that are suggested by the MLE.
# One option that might work is to basically try fixing alpha at 4 different  
# values and just running the MCMC at all of them.
# This may not be the most attractive option, but it does seem to elicit a 
# reasonable estimate for alpha.


data.seed <- data.seed + 1
set.seed(data.seed)  # data
data <- rRareBinarySpat(x, s = s, knots = knots.t, beta = 0, xi = xi.t,
                        alpha = alpha.t, rho = rho.t, prob.success = prop[2])
y <- data$y

ntrain <- floor(0.75 * ns)
ntest  <- ns - ntrain
obs <- c(rep(T, ntrain), rep(F, ntest))
y.o <- y[obs, , drop = FALSE]
X.o <- matrix(x[obs], ntrain, 1)
s.o <- s[obs, ]
y.validate <- y[!obs, , drop = FALSE]
X.p <- matrix(x[!obs, ], ntest, 1)
s.p <- s[!obs, ]

nll <- rep(9999999, length(rhos))
for (i in seq_along(rhos)) {
  temp <- fit.rarebinaryCPP(c(0.5, 0, -4), rho = rhos[i],  y = y.o,
                            dw2 = dw2, d = d, cov = X.o, threads = 6)
  nll[i] <- temp$value
  if (nll[i] == min(nll)) {
    fit <- temp
  }
}

# W <- stdW(makeW(dw2 = dw2, rho = rho.t))
# fit <- fit.rarebinaryCPP(c(0.5, 0, -4), rho=rho.t, y = y.o, dw2 = dw2, d = d, cov=X.o, threads=6)

xibeta.hat <- fit$par[2:3]
xibeta.var <- solve(fit$hessian)[2:3, 2:3]

plot(knots.t, ylim = c(0, 1), xlim = c(0, 1), 
     main = "simulated dataset", xlab="", ylab="")
train1.idx <- which(y[obs] == 1)
test1.idx <- which(y[!obs] == 1) + ntrain  # to get to the testing
points(s[train1.idx, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[test1.idx, ], pch = 21, col = "firebrick4", bg = "firebrick1")

# spatial GEV
mcmc.seed <- 1
set.seed(mcmc.seed)
fit.gev.1 <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[3], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[2], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.01, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 200, rho.attempts = 200,
                  spatial = TRUE, rho.init = fit$rho, rho.upper = 9,
                  alpha.init = 0.3, a.init = 1000, iterplot = TRUE,
                  alpha.fix = TRUE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

fit.gev.2 <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[3], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[2], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.01, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 200, rho.attempts = 200,
                  spatial = TRUE, rho.init = fit$rho, rho.upper = 9,
                  alpha.init = 0.5, a.init = 100, iterplot = TRUE,
                  alpha.fix = TRUE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

fit.gev.3 <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[3], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[2], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.01, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 200, rho.attempts = 200,
                  spatial = TRUE, rho.init = fit$rho, rho.upper = 9,
                  alpha.init = 0.7, a.init = 100, iterplot = TRUE,
                  alpha.fix = TRUE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

fit.gev.4 <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[3], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[2], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.01, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 50, xi.attempts = 50,
                  alpha.attempts = 200, rho.attempts = 200,
                  spatial = TRUE, rho.init = fit$rho, rho.upper = 9,
                  alpha.init = 0.2, a.init = 30000, iterplot = TRUE,
                  alpha.fix = TRUE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

post.prob.gev.1 <- pred.spgev(mcmcoutput = fit.gev.1, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

post.prob.gev.2 <- pred.spgev(mcmcoutput = fit.gev.2, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

post.prob.gev.3 <- pred.spgev(mcmcoutput = fit.gev.3, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

post.prob.gev.4 <- pred.spgev(mcmcoutput = fit.gev.4, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

# spatial logit
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.logit <- spGLM(formula = y.o ~ 1, family = "binomial", coords = s.o,
                   knots = knots, starting = starting, tuning = tuning,
                   priors = priors, cov.model = cov.model,
                   n.samples = iters, verbose = verbose,
                   n.report = n.report)

yp.sp.log <- spPredict(sp.obj = fit.logit, pred.coords = s.p,
                       pred.covars = X.p, start = burn + 1, end = iters,
                       thin = 1, verbose = TRUE, n.report = 500)

post.prob.log <- t(yp.sp.log$p.y.predictive.samples)

# spatial probit
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.probit <- probit(Y = y.o, X = X.o, s = s.o, knots = knots,
                     iters = iters, burn = burn, update = update)

post.prob.pro <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                             s.pred = s.p, knots = knots,
                             start = 1, end = iters - burn, update = 500)


bs.gev.1.2 <- BrierScore(post.prob.gev.1, y.validate)  # 
bs.gev.2.2 <- BrierScore(post.prob.gev.2, y.validate)  # 
bs.gev.3.2 <- BrierScore(post.prob.gev.3, y.validate)  # 
bs.gev.4.2 <- BrierScore(post.prob.gev.4, y.validate)  # 
bs.log.2   <- BrierScore(post.prob.log, y.validate)    # 
bs.pro.2   <- BrierScore(post.prob.pro, y.validate)    # 

# a couple of thoughts:
# In the smaller dataset, the MCMC really seems to find the posterior 
# distribution very shortly after starting.
# This is partly because we've really helped stabilize the MCMC by using a 
# candidate for xi and beta that are suggested by the MLE.
# One option that might work is to basically try fixing alpha at 4 different  
# values and just running the MCMC at all of them.
# This may not be the most attractive option, but it does seem to elicit a 
# reasonable estimate for alpha.