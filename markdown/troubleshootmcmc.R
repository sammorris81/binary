rm(list=ls())

####################################################################
#### loading libraries and necessary files
####################################################################
# libraries
library(fields)
library(Rcpp)
library(evd)
library(spBayes)
library(mvtnorm)
library(numDeriv)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp(file = "../code/R/pairwise.cpp")

source("../code/R/auxfunctions.R", chdir = TRUE)
source("../code/R/updateModel.R")
source("../code/R/mcmc.R")
source("../code/R/probit.R", chdir=T)


####################################################################
#### data setup
####################################################################
# true knots
knots.t <- as.matrix(expand.grid(seq(0.00, 1.00, length=12),
                                 seq(0.00, 1.00, length=12)))

# knots used to fit model
knots <- as.matrix(expand.grid(seq(0.00, 1.00, length=16),
                               seq(0.00, 1.00, length=16)))

# grid for searching over rho
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
alpha.t <- 0.35
alpha.min <- 0
alpha.max <- 1
alpha.rng <- alpha.max - alpha.min
xi.t    <- 0.25
rho.t   <- 0.15
prop    <- c(0.05, 0.01)
knots.h <- knots[2, 1] - knots[1, 1]
rho.init <- log(knots.h)

data.seed <- 3282
set.seed(data.seed)  # data
data <- rRareBinarySpat(x, s = s, knots = knots.t, beta = 0, xi = xi.t,
                        alpha = alpha.t, rho = rho.t, prob.success = prop[1])
y <- data$y

# testing vs training
ntrain <- floor(0.75 * ns)
ntest  <- ns - ntrain
obs   <- c(rep(T, ntrain), rep(F, ntest))
y.o   <- y[obs, , drop = FALSE]
X.o   <- matrix(x[obs], ntrain, 1)
s.o   <- s[obs, ]
dw2.o <- rdist(s.o, knots)
d.o   <- as.matrix(rdist(s.o))
diag(d.o) <- 0
y.validate <- y[!obs, , drop = FALSE]
X.p <- matrix(x[!obs, ], ntest, 1)
s.p <- s[!obs, ]

####################################################################
#### Start MCMC setup: Most of this is used for the spBayes package
####################################################################
iters <- 15000; burn <- 10000; update <- 500; thin <- 1

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

####################################################################
#### Get ML estimates for candidate distribution for GEV
####################################################################
# plot the data
plot(knots.t, ylim = c(0, 1), xlim = c(0, 1),
     main = "simulated dataset", xlab="", ylab="")
train1.idx <- which(y[obs] == 1)
test1.idx <- which(y[!obs] == 1) + ntrain  # to get to the testing
points(s[train1.idx, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[test1.idx, ], pch = 21, col = "firebrick4", bg = "firebrick1")

# get initial values and fit pcl. parameter order: alpha, rho, xi, beta
fit <- fit.rarebinaryCPP(c(0, rho.init, 0, -4), y = y.o, dw2 = dw2.o, d = d.o,
                         cov=X.o, alpha.min = alpha.min, alpha.max = alpha.max,
                         threads=6)

alpha.hat  <- (exp(fit$par[1]) / (1 + exp(fit$par[1]))) *
               alpha.rng + alpha.min
if (alpha.hat < 0.3) {
  alpha.hat <- 0.3
} else if (alpha.hat > 0.9) {
  alpha.hat <- 0.9
}
rho.hat    <- exp(fit$par[2])
xibeta.hat <- fit$par[3:4]
xibeta.var <- solve(fit$hessian[3:4, 3:4])

####################################################################
#### Fit MCMC
####################################################################
# spatial GEV
rho.hat <- rho.t
alpha.hat <- alpha.t
xibeta.hat <- c(xi.t, -data$thresh)
xibeta.var <- matrix(c(0.4, -0.2, -0.2, 1.4))
fit <- list(par=c(alpha.hat, rho.hat, xibeta.hat))

mcmc.seed <- 1
set.seed(mcmc.seed)

# Rprof(filename = "Rprof.out", line.profiling = TRUE)
fit.gev <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                beta.init = fit$par[4], beta.m = 0, beta.s = 100,
                xi.init = fit$par[3], xi.m = 0, xi.s = 0.5,
                knots = knots, beta.tune = 1, xi.tune = 0.1,
                alpha.tune = 0.05, rho.tune = 0.1, A.tune = 1,
                beta.attempts = 50, xi.attempts = 50,
                alpha.attempts = 300, rho.attempts = 100,
                spatial = TRUE, rho.init = rho.hat, rho.upper = 9,
                alpha.init = 0.40, a.init = 1000, iterplot = TRUE,
                alpha.fix = FALSE, rho.fix = TRUE, xibeta.joint = TRUE,
                xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                iters = iters, burn = burn, update = update, thin = 1)
# Rprof(NULL)
# summaryRprof(filename = "Rprof.out", lines = "show")

post.prob.gev.1 <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

# fixing alpha and rho to be the true values
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.gev.t <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[4], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[3], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.05, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 200, xi.attempts = 200,
                  alpha.attempts = 7500, rho.attempts = 100,
                  spatial = TRUE, rho.init = rho.t, rho.upper = 9,
                  alpha.init = alpha.t, a.init = 10000, iterplot = TRUE,
                  alpha.fix = FALSE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

post.prob.gev.1t <- pred.spgev(mcmcoutput = fit.gev.t, x.pred = X.p,
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

post.prob.log.1 <- t(yp.sp.log$p.y.predictive.samples)

# spatial probit
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.probit <- probit(Y = y.o, X = X.o, s = s.o, knots = knots,
                     iters = iters, burn = burn, update = update)

post.prob.pro.1 <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.p, knots = knots,
                               start = 1, end = iters - burn, update = 500)



####################################################################
#### Get Brier scores
####################################################################
bs.gev.1  <- BrierScore(post.prob.gev.1, y.validate)   # 0.0474 - before change
                                                       # 0.0377 - after change
bs.gev.1t <- BrierScore(post.prob.gev.1t, y.validate)  # 0.0470
bs.log.1  <- BrierScore(post.prob.log.1, y.validate)   # 0.0481
bs.pro.1  <- BrierScore(post.prob.pro.1, y.validate)   # 0.0360
dic.gev.1 <- dic.spgev(mcmcoutput = fit.gev, y = y.o, x = X.o, dw2 = dw2.o)  # -738
dic.gev.1t <- dic.spgev(mcmcoutput = fit.gev.t, y = y.o, x = X.o, dw2 = dw2.o)  # 134
dic.log.1 <- spDiag(fit.logit, start = burn + 1, end = iters)  # 446.68
dic.pro.1 <- dic.spprob(mcmcoutput = fit.probit, Y = y.o, X = X.o, s = s.o,
                        knots = knots)  # 268


####################################################################
#### Try with fewer 1s
####################################################################
data.seed <- data.seed + 1
set.seed(data.seed)  # data
data <- rRareBinarySpat(x, s = s, knots = knots.t, beta = 0, xi = xi.t,
                        alpha = alpha.t, rho = rho.t, prob.success = prop[2])
y <- data$y

# testing vs training
ntrain <- floor(0.75 * ns)
ntest  <- ns - ntrain
obs   <- c(rep(T, ntrain), rep(F, ntest))
y.o   <- y[obs, , drop = FALSE]
X.o   <- matrix(x[obs], ntrain, 1)
s.o   <- s[obs, ]
dw2.o <- rdist(s.o, knots)
d.o   <- as.matrix(rdist(s.o))
diag(d.o) <- 0
y.validate <- y[!obs, , drop = FALSE]
X.p <- matrix(x[!obs, ], ntest, 1)
s.p <- s[!obs, ]

####################################################################
#### Start MCMC setup: Most of this is used for the spBayes package
####################################################################
iters <- 15000; burn <- 10000; update <- 500; thin <- 1

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

####################################################################
#### Get ML estimates for candidate distribution for GEV
####################################################################
# plot the data
plot(knots.t, ylim = c(0, 1), xlim = c(0, 1),
     main = "simulated dataset", xlab="", ylab="")
train1.idx <- which(y[obs] == 1)
test1.idx <- which(y[!obs] == 1) + ntrain  # to get to the testing
points(s[train1.idx, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[test1.idx, ], pch = 21, col = "firebrick4", bg = "firebrick1")

# get initial values and fit pcl. parameter order: alpha, rho, xi, beta
fit <- fit.rarebinaryCPP(c(0, rho.init, 0, -4), y = y.o, dw2 = dw2.o, d = d.o,
                         cov=X.o, alpha.min = alpha.min, alpha.max = alpha.max,
                         threads=6)

alpha.hat  <- (exp(fit$par[1]) / (1 + exp(fit$par[1]))) *
               alpha.rng + alpha.min
if (alpha.hat < 0.3) {
  alpha.hat <- 0.3
} else if (alpha.hat > 0.9) {
  alpha.hat <- 0.9
}
rho.hat    <- exp(fit$par[2])
xibeta.hat <- fit$par[3:4]
xibeta.var <- solve(fit$hessian[3:4, 3:4])

####################################################################
#### Fit MCMC
####################################################################
# spatial GEV
mcmc.seed <- 1
set.seed(mcmc.seed)
fit.gev <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                beta.init = fit$par[4], beta.m = 0, beta.s = 100,
                xi.init = fit$par[3], xi.m = 0, xi.s = 0.5,
                knots = knots, beta.tune = 1, xi.tune = 0.1,
                alpha.tune = 0.05, rho.tune = 0.1, A.tune = 1,
                beta.attempts = 200, xi.attempts = 200,
                alpha.attempts = 7500, rho.attempts = 100,
                spatial = TRUE, rho.init = rho.hat, rho.upper = 9,
                alpha.init = 0.40, a.init = 10000, iterplot = TRUE,
                alpha.fix = FALSE, rho.fix = TRUE, xibeta.joint = TRUE,
                xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                iters = iters, burn = burn, update = update, thin = 1)

post.prob.gev.2 <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                              s.pred = s.p, knots = knots,
                              start = 1, end = iters - burn, update = 500)

# fixing alpha and rho to be the true values
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.gev.t <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                  beta.init = fit$par[4], beta.m = 0, beta.s = 100,
                  xi.init = fit$par[3], xi.m = 0, xi.s = 0.5,
                  knots = knots, beta.tune = 1, xi.tune = 0.1,
                  alpha.tune = 0.05, rho.tune = 0.1, A.tune = 1,
                  beta.attempts = 200, xi.attempts = 200,
                  alpha.attempts = 7500, rho.attempts = 100,
                  spatial = TRUE, rho.init = rho.t, rho.upper = 9,
                  alpha.init = alpha.t, a.init = 10000, iterplot = TRUE,
                  alpha.fix = FALSE, rho.fix = TRUE, xibeta.joint = TRUE,
                  xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                  iters = iters, burn = burn, update = update, thin = 1)

post.prob.gev.2t <- pred.spgev(mcmcoutput = fit.gev.t, x.pred = X.p,
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

post.prob.log.2 <- t(yp.sp.log$p.y.predictive.samples)

# spatial probit
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.probit <- probit(Y = y.o, X = X.o, s = s.o, knots = knots,
                     iters = iters, burn = burn, update = update)

post.prob.pro.2 <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                             s.pred = s.p, knots = knots,
                             start = 1, end = iters - burn, update = 500)

####################################################################
#### Get Brier scores
####################################################################
bs.gev.2  <- BrierScore(post.prob.gev.2, y.validate)   # 0.0187
bs.gev.2t <- BrierScore(post.prob.gev.2t, y.validate)  # 0.0187
bs.log.2  <- BrierScore(post.prob.log.2, y.validate)   # 0.0185
bs.pro.2  <- BrierScore(post.prob.pro.2, y.validate)   # 0.0147
dic.gev.2 <- dic.spgev(mcmcoutput = fit.gev, y = y.o, x = X.o, dw2 = dw2.o)  # 26.7
dic.gev.2t <- dic.spgev(mcmcoutput = fit.gev.t, y = y.o, x = X.o, dw2 = dw2.o)  # -33.4
dic.log.2 <- spDiag(fit.logit, start = burn + 1, end = iters)  # 96.68
dic.pro.2 <- dic.spprob(mcmcoutput = fit.probit, Y = y.o, X = X.o, s = s.o,
                        knots = knots)  # 60.8

####################################################################
#### Try when the occurrences are only in a certain location
####################################################################
data.seed <- data.seed + 1
set.seed(data.seed)  # data
center <- c(runif(1), runif(1))
radius <- 0.16
s.dist <- sqrt((s[, 1] - center[1])^2 + (s[, 2] - center[2])^2)
mean(s.dist < radius)
prob <- (s.dist < radius) * 0.85 + (1 - s.dist < radius) * 0.001
y    <- matrix(rbinom(n = ns, size = 1, prob = prob), ns, 1)

# testing vs training
ntrain <- floor(0.75 * ns)
ntest  <- ns - ntrain
obs   <- c(rep(T, ntrain), rep(F, ntest))
y.o   <- y[obs, , drop = FALSE]
X.o   <- matrix(x[obs], ntrain, 1)
s.o   <- s[obs, ]
dw2.o <- rdist(s.o, knots)
d.o   <- as.matrix(rdist(s.o))
diag(d.o) <- 0
y.validate <- y[!obs, , drop = FALSE]
X.p <- matrix(x[!obs, ], ntest, 1)
s.p <- s[!obs, ]

####################################################################
#### Start MCMC setup: Most of this is used for the spBayes package
####################################################################
iters <- 15000; burn <- 10000; update <- 500; thin <- 1

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

####################################################################
#### Get ML estimates for candidate distribution for GEV
####################################################################
# plot the data
plot(knots.t, ylim = c(0, 1), xlim = c(0, 1),
     main = "simulated dataset", xlab="", ylab="")
train1.idx <- which(y[obs] == 1)
test1.idx <- which(y[!obs] == 1) + ntrain  # to get to the testing
points(s[train1.idx, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[test1.idx, ], pch = 21, col = "firebrick4", bg = "firebrick1")

# get initial values and fit pcl. parameter order: alpha, rho, xi, beta
fit <- fit.rarebinaryCPP(c(0, rho.init, 0, -4), y = y.o, dw2 = dw2.o, d = d.o,
                         cov=X.o, alpha.min = alpha.min, alpha.max = alpha.max,
                         threads=6)

alpha.hat  <- (exp(fit$par[1]) / (1 + exp(fit$par[1]))) *
               alpha.rng + alpha.min
if (alpha.hat < 0.3) {
  alpha.hat <- 0.3
} else if (alpha.hat > 0.9) {
  alpha.hat <- 0.9
}
rho.hat    <- exp(fit$par[2])
xibeta.hat <- fit$par[3:4]
xibeta.var <- solve(fit$hessian[3:4, 3:4])

####################################################################
#### Fit MCMC
####################################################################
# spatial GEV
mcmc.seed <- 1
set.seed(mcmc.seed)
fit.gev <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                beta.init = fit$par[4], beta.m = 0, beta.s = 100,
                xi.init = fit$par[3], xi.m = 0, xi.s = 0.5,
                knots = knots, beta.tune = 1, xi.tune = 0.1,
                alpha.tune = 0.05, rho.tune = 0.1, A.tune = 1,
                beta.attempts = 200, xi.attempts = 200,
                alpha.attempts = 7500, rho.attempts = 100,
                spatial = TRUE, rho.init = rho.hat, rho.upper = 9,
                alpha.init = 0.40, a.init = 10000, iterplot = TRUE,
                alpha.fix = FALSE, rho.fix = TRUE, xibeta.joint = TRUE,
                xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                iters = iters, burn = burn, update = update, thin = 1)

post.prob.gev.3 <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
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

post.prob.log.3 <- t(yp.sp.log$p.y.predictive.samples)

# spatial probit
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.probit <- probit(Y = y.o, X = X.o, s = s.o, knots = knots,
                     iters = iters, burn = burn, update = update)

post.prob.pro.3 <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.p, knots = knots,
                               start = 1, end = iters - burn, update = 500)

####################################################################
#### Get Brier scores
####################################################################
bs.gev.3 <- BrierScore(post.prob.gev.3, y.validate)   # 0.0237
bs.log.3 <- BrierScore(post.prob.log.3, y.validate)   # 0.0195
bs.pro.3 <- BrierScore(post.prob.pro.3, y.validate)   # 0.0198
dic.gev.3 <- dic.spgev(mcmcoutput = fit.gev, y = y.o, x = X.o, dw2 = dw2.o)  # -7328
dic.log.3 <- spDiag(fit.logit, start = burn + 1, end = iters)  # 138.7
dic.pro.3 <- dic.spprob(mcmcoutput = fit.probit, Y = y.o, X = X.o, s = s.o,
                        knots = knots)  # 151

####################################################################
#### Try when the occurrences are only in a certain location
####################################################################
data.seed <- data.seed + 1
set.seed(data.seed)  # data
center <- c(runif(1), runif(1))
radius <- 0.13
s.dist <- sqrt((s[, 1] - center[1])^2 + (s[, 2] - center[2])^2)
mean(s.dist < radius)
prob <- (s.dist < radius) * 0.90 + (1 - s.dist < radius) * 0.001
y    <- matrix(rbinom(n = ns, size = 1, prob = prob), ns, 1)

# testing vs training
ntrain <- floor(0.75 * ns)
ntest  <- ns - ntrain
obs   <- c(rep(T, ntrain), rep(F, ntest))
y.o   <- y[obs, , drop = FALSE]
X.o   <- matrix(x[obs], ntrain, 1)
s.o   <- s[obs, ]
dw2.o <- rdist(s.o, knots)
d.o   <- as.matrix(rdist(s.o))
diag(d.o) <- 0
y.validate <- y[!obs, , drop = FALSE]
X.p <- matrix(x[!obs, ], ntest, 1)
s.p <- s[!obs, ]

####################################################################
#### Start MCMC setup: Most of this is used for the spBayes package
####################################################################
iters <- 15000; burn <- 10000; update <- 500; thin <- 1

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

####################################################################
#### Get ML estimates for candidate distribution for GEV
####################################################################
# plot the data
plot(knots.t, ylim = c(0, 1), xlim = c(0, 1),
     main = "simulated dataset", xlab="", ylab="")
train1.idx <- which(y[obs] == 1)
test1.idx <- which(y[!obs] == 1) + ntrain  # to get to the testing
points(s[train1.idx, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[test1.idx, ], pch = 21, col = "firebrick4", bg = "firebrick1")

# get initial values and fit pcl. parameter order: alpha, rho, xi, beta
fit <- fit.rarebinaryCPP(c(0, rho.init, 0, -4), y = y.o, dw2 = dw2.o, d = d.o,
                         cov=X.o, alpha.min = alpha.min, alpha.max = alpha.max,
                         threads=6)

alpha.hat  <- (exp(fit$par[1]) / (1 + exp(fit$par[1]))) *
               alpha.rng + alpha.min
if (alpha.hat < 0.3) {
  alpha.hat <- 0.3
} else if (alpha.hat > 0.9) {
  alpha.hat <- 0.9
}
rho.hat    <- exp(fit$par[2])

xibeta.hat <- fit$par[3:4]
xibeta.var <- solve(fit$hessian[3:4, 3:4])

####################################################################
#### Fit MCMC
####################################################################
# spatial GEV
mcmc.seed <- 1
set.seed(mcmc.seed)
fit.gev <- mcmc(y = y.o, s = s.o, x = X.o, s.pred = NULL, x.pred = NULL,
                beta.init = fit$par[4], beta.m = 0, beta.s = 100,
                xi.init = fit$par[3], xi.m = 0, xi.s = 0.5,
                knots = knots, beta.tune = 1, xi.tune = 0.1,
                alpha.tune = 0.05, rho.tune = 0.1, A.tune = 1,
                beta.attempts = 200, xi.attempts = 200,
                alpha.attempts = 7500, rho.attempts = 100,
                spatial = TRUE, rho.init = rho.hat, rho.upper = 9,
                alpha.init = 0.40, a.init = 10000, iterplot = TRUE,
                alpha.fix = FALSE, rho.fix = TRUE, xibeta.joint = TRUE,
                xibeta.hat = xibeta.hat, xibeta.var = xibeta.var,
                iters = iters, burn = burn, update = update, thin = 1)

post.prob.gev.4 <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
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

post.prob.log.4 <- t(yp.sp.log$p.y.predictive.samples)

# spatial probit
mcmc.seed <- mcmc.seed + 1
set.seed(mcmc.seed)
fit.probit <- probit(Y = y.o, X = X.o, s = s.o, knots = knots,
                     iters = iters, burn = burn, update = update)

post.prob.pro.4 <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p,
                               s.pred = s.p, knots = knots,
                               start = 1, end = iters - burn, update = 500)

####################################################################
#### Get Brier scores
####################################################################
bs.gev.4 <- BrierScore(post.prob.gev.4, y.validate)   # 0.0159
bs.log.4 <- BrierScore(post.prob.log.4, y.validate)   # 0.0166
bs.pro.4 <- BrierScore(post.prob.pro.4, y.validate)   # 0.0146
dic.gev.4 <- dic.spgev(mcmcoutput = fit.gev, y = y.o, x = X.o, dw2 = dw2.o)  # -4522
dic.log.4 <- spDiag(fit.logit, start = burn + 1, end = iters)  # 78.6
dic.pro.4 <- dic.spprob(mcmcoutput = fit.probit, Y = y.o, X = X.o, s = s.o,
                        knots = knots)  # 73.4

####################################################################
#### Find DIC
####################################################################
dic.log.4 <- spDiag(sp.obj = fit.logit, start = burn + 1, end = iters,
                    thin = 1, verbose = TRUE, n.report = 500)