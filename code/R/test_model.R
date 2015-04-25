# model update testing
# independence:
#   beta: PASS
#   xi: PASS
#   beta and xi: PASS
# spatial dependence:
#   beta: PASS
#   xi: PASS
#   beta and xi: PASS
#   alpha: PASS
#   rho: PASS
#   alpha and rho: PASS (both at alpha = 0.4 and alpha = 0.8)
#       moves fairly slowly in the chain though
#   A:
#   beta, xi, alpha:
#   beta, xi, rho:
#   beta, xi, alpha, rho

# currently trying to figure out what's wrong with the dependent version
# of the model...issues with beta, xi, alpha, and rho updates when data
# are created using the model with dependence.

rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)

source("auxfunctions.R")
source("updateModel.R")
source("mcmc.R")
set.seed(10)
ns   <- 2000
nt   <- 1
s    <- cbind(runif(ns, 0, 6), runif(ns, 0, 6))
x <- array(1, dim=c(ns, nt, 3))
x[, , 2] <- s[, 1]
x[, , 3] <- s[, 2]
knots.1 <- seq(0, 10, length=12)
knots.2 <- seq(0, 10, length=12)
knots <- expand.grid(knots.1, knots.2)
# knots  <- s
nknots <- nrow(knots)

set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
source("mcmc.R")
iters <- 100000
burn  <- 80000
xi.t <- 0.25
beta.t <- c(0, 0, 0)
alpha.t <- 0.7
rho.t   <- 0.5
data <- rRareBinarySpat(x, s = s, knots = knots, beta = beta.t,
                        xi = xi.t, alpha = alpha.t, rho = rho.t,
                        prob.success = 0.01)

occurs <- which(data$y == 1)
plot(s[occurs,], xlim=c(0, 6), ylim=c(0, 6))
# points(knots, col=2)

obs <- c(rep(T, 1000), rep(F, 1000))
# obs <- rep(T, 350)
y.o <- data$y[obs, , drop=F]
s.o <- s[obs, ]
x.o <- x[obs, , , drop=F]

set.seed(1)
tic <- proc.time()[3]
fit.1 <- mcmc(y=y.o, s=s.o, x=x.o, knots=knots, npts=70, rho.upper=15,
              rho.init=5, iters=iters, burn=burn,
              beta.tune=1, xi.tune=1,
              alpha.tune=0.05, rho.tune=0.05, A.tune=1,
              # beta.attempts=100, xi.attempts=100, alpha.attempts=100, rho.attempts=100,
              # beta.init=beta.t, alpha.init=alpha.t, a.init=data$a,
              update=500, iterplot=TRUE)
toc <- proc.time()[3]
toc - tic

save.image(file="datasets/rarebinary2.RData")

# testing spbays
n.report <- 500
verbose <- TRUE
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1, "nu"=0.1,
               "beta"=c(0.1, 0.1, 0.1), "w"=0.1)
starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1, "nu"=0.5,
                 "beta"=c(0, 0, 0), "w"=0)
priors <- list("beta.norm"=list(rep(0,3), diag(1000,3)),
               "phi.unif"=c(0.5, 1e4), "sigma.sq.ig"=c(1, 1),
               "tau.sq.ig"=c(1, 1), "nu.unif"=c(1e-4, 5))
cov.model <- "matern"

y.o <- data$y[obs, ]
s.o <- s[obs, ]
x.o <- x[obs, , ]

set.seed(2)
tic <- proc.time()[3]
fit.2 <- spGLM(y.o~x.o-1, family="binomial", coords=s.o, knots=c(9, 9, 0),
           starting=starting, tuning=tuning, priors=priors, cov.model=cov.model,
           n.samples=iters, verbose=verbose, n.report=n.report)
toc <- proc.time()[3]
toc - tic

save.image(file="datasets/spbayes2.RData")

load("datasets/rarebinary2.RData")
load("datasets/spbayes2.RData")
source("auxfunctions.R")
source("updateModel.R")
source("mcmc.R")

s.p <- s[!obs, ]
np <- nrow(s.p)
x.p.rb  <- x[!obs, , , drop=FALSE]
x.p.spb <- matrix(x[!obs, , ], np, 3)  # need vector for spbayes

set.seed(3)
yp.rb <- predictProb(mcmcoutput = fit.1, x.pred = x.p.rb, s.pred = s.p,
                     knots = knots, start = 1, end=20000, update=500)

set.seed(4)
yp.spbayes <- spPredict(sp.obj = fit.2, pred.coords = s.p,
                        pred.covars = x.p.spb, start = 80001, end = 100000,
                        thin = 1, verbose = TRUE, n.report = 500)
yp.spb <- t(yp.spbayes$p.y.predictive.samples)

save.image(file="datasets/predictions2.RData")
load(file="datasets/predictions2.RData")

y.validate   <- data$y[!obs, ]
success.idx  <- which(y.validate == 1)
rb.prob.med  <- apply(yp.rb, 2, quantile, probs=0.50)
spb.prob.med <- apply(yp.spb, 2, quantile, probs=0.50)

rb.bs  <- BrierScore(post.prob = yp.rb, validate = y.validate)
spb.bs <- BrierScore(post.prob = yp.spb, validate = y.validate)
# when there are 350 training sites and 650 validation sites
# xi.t <- 0.1
# beta.t <- c(0, -1, 0)
# alpha.t <- 0.5
# rb.bs  = 0.0321
# spb.bs = 0.0319
# predictions are saved as predictions1.RData

# when there are 1000 training sites and 1000 validation sites
# xi.t <- 0.1
# beta.t <- c(0, -1, 0)
# alpha.t <- 0.3
# rb.bs  = 0.0150
# spb.bs = 0.0167
# predictions are saved as predictions2.RData

# test timing for whole mcmc
source("auxfunctions.R")
source("updateModel.R")
set.seed(10)
ns   <- 200
nt   <- 5
s    <- cbind(runif(ns, 0, 10), runif(ns, 0, 10))
x <- array(1, dim=c(ns, nt, 3))
x[, , 2] <- s[, 1]
x[, , 3] <- s[, 2]
knots.1 <- seq(0, 10, length=9)
knots.2 <- seq(0, 10, length=9)
knots <- expand.grid(knots.1, knots.2)
# knots  <- s
nknots <- nrow(knots)

set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
source("mcmc.R")
iters <- 2000
burn  <- 200
xi.t <- 0.1
beta.t <- c(1, -1, 0)
alpha.t <- 0.3
rho.t   <- 3
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t)

tic <- proc.time()[3]
fit <- mcmc(y=data$y, s=s, x=x, knots=knots, npts=70, rho.upper=15,
            rho.init=5, iters=iters, burn=burn,
            beta.tune=0.05, xi.tune=0.1,
            alpha.tune=0.05, rho.tune=0.05, A.tune=1,
            alpha.attempts=50, rho.attempts=50,
            # beta.init=beta.t, alpha.init=alpha.t, a.init=data$a,
            update=100, iterplot=FALSE)
toc <- proc.time()[3]
toc - tic


# test update for beta
source("auxfunctions.R")
source("updateModel.R")
set.seed(10)
ns   <- 200
nt   <- 40
s    <- cbind(runif(ns, 0.01, 0.99), runif(ns, 0.01, 0.99))
x <- array(1, dim=c(ns, nt, 3))
x[, , 2] <- s[, 1]
x[, , 3] <- s[, 2]
knots.1 <- seq(0, 10, length=9)
knots.2 <- seq(0, 10, length=9)
knots <- expand.grid(knots.1, knots.2)
# knots  <- s
nknots <- nrow(knots)

set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 2000
burn  <- 200
xi.t <- 0.1
beta.t <- c(1, -1, 0)
theta.star.t <- 1
alpha.t <- 1
thresh <- 0
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t, thresh=thresh)

# initialization and prior distribution
beta <- c(0, 0, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(0, ns, nt)
z <- getZ(xi=xi.t, x.beta=x.beta)
z.star <- z^(1 / alpha.t)
cur.lly <- logLikeY(y=y, theta.star=theta.star.t, z.star=z.star)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(0.1, 3)

# storage
beta.keep <- matrix(NA, nreps, 3)
tic <- proc.time()[3]
for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta.star=theta.star.t, alpha=alpha.t,
                            z=z, z.star=z.star, beta=beta, beta.m=beta.m,
                            beta.s=beta.s, x.beta=x.beta, xi=xi.t, x=x,
                            cur.lly=cur.lly, acc=acc.beta, att=att.beta,
                            mh=mh.beta, thresh=thresh)
  beta     <- beta.update$beta
  x.beta   <- beta.update$x.beta
  z        <- beta.update$z
  z.star   <- beta.update$z.star
  cur.lly  <- beta.update$cur.lly
  att.beta <- beta.update$att
  acc.beta <- beta.update$acc

  beta.keep[i, ] <- beta
  if (i %% 500 == 0) {
    print(i)
  }

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta)
    acc.beta  <- mh.update$acc
    att.beta  <- mh.update$att
    mh.beta   <- mh.update$mh
  }
}
toc <- proc.time()[3]
toc - tic

par(mfrow=c(1, 3))
plot(beta.keep[, 1], type="l", main="beta 0")
plot(beta.keep[, 2], type="l", main="beta 1")
plot(beta.keep[, 3], type="l", main="beta 2")


# test update for xi
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 10000
burn  <- 7000
xi.t <- 0.1
beta.t <- c(1, -1, 0)
theta.star.t <- 1
alpha.t <- 1
thresh <- 0
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t, thresh=thresh)
x.beta.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta.t[, t] <- x[, t, ] %*% beta.t
}

# initialization and prior distribution
xi <- 0
xi.m <- 0
xi.s <- 10
z <- getZ(xi=xi, x.beta=x.beta.t)
z.star <- z^(1 / alpha.t)
cur.lly <- logLikeY(y=y, theta.star=theta.star.t, z.star=z.star)

# MH adjustments
acc.xi <- att.xi <- mh.xi <- rep(1, 3)

# storage
xi.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  xi.update <- updateXi(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z,
                        z.star=z.star, x.beta=x.beta.t, xi=xi, xi.m=xi.m,
                        xi.s=xi.s, cur.lly=cur.lly,
                        acc=acc.xi, att=att.xi, mh=mh.xi, thresh=thresh)
  xi      <- xi.update$xi
  z       <- xi.update$z
  z.star  <- xi.update$z.star
  cur.lly <- xi.update$cur.lly
  att.xi  <- xi.update$att
  acc.xi  <- xi.update$acc

  xi.keep[i] <- xi
  if (i %% 500 == 0) {
    if (i > 3000) {
      start <- i - 3000
    } else {
      start <- 1
    }
    print(i)
    plot(xi.keep[start:i], type="l")
  }

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.xi, att=att.xi, mh=mh.xi)
    acc.xi  <- mh.update$acc
    att.xi  <- mh.update$att
    mh.xi   <- mh.update$mh
  }
}

# test update for xi and beta
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 7000
burn  <- 5000
xi.t <- 0.1
beta.t <- c(1, -1, 0)
theta.star.t <- 1
alpha.t <- 1
x.beta.t <- matrix(2, ns, nt)
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t)

# initialization and prior distribution
beta <- c(1, -1, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta[, t] <- x[, t, ] %*% beta
}
xi <- 0.1
xi.m <- 0
xi.s <- 1
z <- getZ(xi=xi, x.beta=x.beta)
z.star <- z^(1 / alpha.t)
cur.lly <- logLikeY(y=y, theta.star=theta.star.t, z.star=z.star)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(0.01, 3)
acc.xi <- att.xi <- mh.xi <- 0.1

# storage
beta.keep <- matrix(NA, nreps, 3)
xi.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta.star=theta.star.t, alpha=alpha.t,
                            z=z, z.star=z.star, beta=beta, beta.m=beta.m,
                            beta.s=beta.s, xi=xi.t, x=x, cur.lly=cur.lly,
                            acc=acc.beta, att=att.beta, mh=mh.beta)
  beta     <- beta.update$beta
  x.beta   <- beta.update$x.beta
  z        <- beta.update$z
  z.star   <- beta.update$z.star
  cur.lly  <- beta.update$cur.lly
  att.beta <- beta.update$att
  acc.beta <- beta.update$acc

  beta.keep[i, ] <- beta


  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta)
    acc.beta  <- mh.update$acc
    att.beta  <- mh.update$att
    mh.beta   <- mh.update$mh
  }

  xi.update <- updateXi(y=y, theta.star=theta.star.t, alpha=alpha.t,
                        z=z, z.star=z.star, x.beta=x.beta,
                        xi=xi, xi.m=xi.m, xi.s=xi.s,
                        cur.lly=cur.lly,
                        acc=acc.xi, att=att.xi, mh=mh.xi)
  xi      <- xi.update$xi
  z       <- xi.update$z
  z.star  <- xi.update$z.star
  cur.lly <- xi.update$cur.lly
  att.xi  <- xi.update$att
  acc.xi  <- xi.update$acc

  xi.keep[i] <- xi

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.xi, att=att.xi, mh=mh.xi)
    acc.xi  <- mh.update$acc
    att.xi  <- mh.update$att
    mh.xi   <- mh.update$mh
  }

  if (i %% 500 == 0) {
    print(i)
    if (i < 3000) {
      start <- 1
    } else {
      start <- i - 3000
    }
    par(mfrow=c(2, 2))
    plot(beta.keep[start:i, 1], type="l", main="beta 0",
         ylab=round(mh.beta[1], 3))
    plot(beta.keep[start:i, 2], type="l", main="beta 1",
         ylab=round(mh.beta[2], 3))
    plot(beta.keep[start:i, 3], type="l", main="beta 2",
         ylab=round(mh.beta[3], 3))
    plot(xi.keep[start:i], type="l", main="xi",
         ylab=round(mh.xi, 3))
  }
}

# test beta update when there's spatial dependence
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps   <- 20000
burn    <- 10000
xi.t    <- 0.1
beta.t  <- c(1, -1, 0)
alpha.t <- 0.7
rho.t   <- 1
thresh  <- 0
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t, thresh=thresh)

dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, rho=rho.t))
w.star.t <- w.t^(1 / alpha.t)
theta.star.t <- getThetaStar(w.star=w.star.t, a=data$a)

# initialization and prior distribution
beta <- c(0, 0, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(0, ns, nt)
z <- getZ(xi=xi.t, x.beta=x.beta)
z.star <- z^(1 / alpha.t)
cur.lly <- logLikeY(y=data$y, theta.star=theta.star.t, z.star=z.star)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(0.1, 3)

# storage
beta.keep <- matrix(NA, nreps, 3)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=data$y, theta.star=theta.star.t, alpha=alpha.t,
                            z=z, z.star=z.star, beta=beta, beta.m=beta.m,
                            beta.s=beta.s, xi=xi.t, x=x, cur.lly=cur.lly,
                            acc=acc.beta, att=att.beta, mh=mh.beta, thresh=0)
  beta     <- beta.update$beta
  x.beta   <- beta.update$x.beta
  z        <- beta.update$z
  z.star   <- beta.update$z.star
  cur.lly  <- beta.update$cur.lly
  att.beta <- beta.update$att
  acc.beta <- beta.update$acc

  beta.keep[i, ] <- beta
  if (i %% 500 == 0) {
    print(i)
    if (i < 3000) {
      start <- 1
    } else {
      start <- i - 3000
    }
    par(mfrow=c(1, 3))
    plot(beta.keep[start:i, 1], type="l", main="beta 0")
    plot(beta.keep[start:i, 2], type="l", main="beta 1")
    plot(beta.keep[start:i, 3], type="l", main="beta 2")
  }

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta)
    acc.beta  <- mh.update$acc
    att.beta  <- mh.update$att
    mh.beta   <- mh.update$mh
  }
}

# test update for xi with spatial dependence
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps   <- 20000
burn    <- 10000
xi.t    <- 0.1
beta.t  <- c(1, -1, 0)
alpha.t <- 0.7
rho.t   <- 1
thresh  <- 0
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t)
dw2  <- as.matrix(rdist(s, knots))^2
w.t  <- stdW(makeW(dw2=dw2, rho=rho.t))
w.star.t <- w.t^(1 / alpha.t)
theta.star.t <- getThetaStar(w.star=w.star.t, a=data$a)

# initialization and prior distribution
xi   <- 0.1
xi.m <- 0
xi.s <- 1
x.beta.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta.t[, t] <- x[, t, ] %*% beta.t
}
z <- getZ(xi=xi, x.beta=x.beta.t)
z.star <- z^(1 / alpha.t)
cur.lly <- logLikeY(y=data$y, theta.star=theta.star.t, z.star=z.star)

# MH adjustments
acc.xi <- att.xi <- mh.xi <- 0.01

# storage
xi.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  xi.update <- updateXi(y=data$y, theta.star=theta.star.t, alpha=alpha.t,
                        z=z, z.star=z.star, x.beta=x.beta.t, xi=xi, xi.m=xi.m,
                        xi.s=xi.s, cur.lly=cur.lly,
                        acc=acc.xi, att=att.xi, mh=mh.xi, thresh=thresh)
  xi      <- xi.update$xi
  z       <- xi.update$z
  z.star  <- xi.update$z.star
  cur.lly <- xi.update$cur.lly
  att.xi  <- xi.update$att
  acc.xi  <- xi.update$acc

  xi.keep[i] <- xi
  if (i %% 1000 == 0) {
    print(i)
    if (i < 3000) {
      start <- 1
    } else {
      start <- i - 3000
    }
    plot(xi.keep[start:i], type="l")
  }

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.xi, att=att.xi, mh=mh.xi)
    acc.xi  <- mh.update$acc
    att.xi  <- mh.update$att
    mh.xi   <- mh.update$mh
  }
}

# test beta and xi with spatial dependence
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps   <- 20000
burn    <- 10000
xi.t    <- 0.1
beta.t  <- c(1, -1, 0)
alpha.t <- 0.7
rho.t   <- 1
thresh  <- 0
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t, thresh=thresh)

dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, rho=rho.t))
w.star.t <- w.t^(1 / alpha.t)
theta.star.t <- getThetaStar(w.star=w.star.t, a=data$a)

# initialization and prior distribution
beta   <- c(0, 0, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(0, ns, nt)
xi     <- 0
xi.m   <- 0
xi.s   <- 1
z <- getZ(xi=xi, x.beta=x.beta)
z.star <- z^(1 / alpha.t)
cur.lly <- logLikeY(y=data$y, theta.star=theta.star.t, z.star=z.star)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(0.1, 3)
acc.xi <- att.xi <- mh.xi <- 0.1

# storage
beta.keep <- matrix(NA, nreps, 3)
xi.keep   <- rep(NA, nreps)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=data$y, theta.star=theta.star.t, alpha=alpha.t,
                            z=z, z.star=z.star, beta=beta, beta.m=beta.m,
                            beta.s=beta.s, xi=xi, x=x, cur.lly=cur.lly,
                            acc=acc.beta, att=att.beta, mh=mh.beta,
                            thresh=thresh)
  beta     <- beta.update$beta
  x.beta   <- beta.update$x.beta
  z        <- beta.update$z
  z.star   <- beta.update$z.star
  cur.lly  <- beta.update$cur.lly
  att.beta <- beta.update$att
  acc.beta <- beta.update$acc

  beta.keep[i, ] <- beta

  xi.update <- updateXi(y=data$y, theta.star=theta.star.t, alpha=alpha.t,
                        z=z, z.star=z.star, x.beta=x.beta, xi=xi, xi.m=xi.m,
                        xi.s=xi.s, cur.lly=cur.lly,
                        acc=acc.xi, att=att.xi, mh=mh.xi, thresh=thresh)
  xi      <- xi.update$xi
  z       <- xi.update$z
  z.star  <- xi.update$z.star
  cur.lly <- xi.update$cur.lly
  att.xi  <- xi.update$att
  acc.xi  <- xi.update$acc

  xi.keep[i] <- xi

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta)
    acc.beta  <- mh.update$acc
    att.beta  <- mh.update$att
    mh.beta   <- mh.update$mh

    mh.update <- mhUpdate(acc=acc.xi, att=att.xi, mh=mh.xi)
    acc.xi    <- mh.update$acc
    att.xi    <- mh.update$att
    mh.xi     <- mh.update$mh
  }

  if (i %% 500 == 0) {
    print(i)
    if (i < 3000) {
      start <- 1
    } else {
      start <- i - 3000
    }
    par(mfrow=c(2, 2))
    plot(beta.keep[start:i, 1], type="l", main="beta 0")
    plot(beta.keep[start:i, 2], type="l", main="beta 1")
    plot(beta.keep[start:i, 3], type="l", main="beta 2")
    plot(xi.keep[start:i], type="l", main="xi")
  }

}

# test update for alpha
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps   <- 20000
burn    <- 10000
xi.t    <- 0.1
beta.t  <- c(1, -1, 0)
alpha.t <- 0.8
rho.t   <- 1
thresh  <- 0
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t, thresh=thresh)
x.beta.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta.t[, t] <- x[, t, ] %*% beta.t
}
z.t <- getZ(xi=xi.t, x.beta=x.beta.t, thresh=thresh)

# initialization and prior distribution
npts <- 70
u.beta     <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
bin.width  <- u.beta[-1] - u.beta[-(npts + 1)]
dw2        <- as.matrix(rdist(s, knots))^2
w.t        <- stdW(makeW(dw2=dw2, rho=rho.t))
alpha      <- 0.5
z.star     <- z.t^(1 / alpha)
w.star     <- w.t^(1 / alpha)
theta.star <- getThetaStar(w.star=w.star.t, a=data$a)

cur.lly  <- logLikeY(y=data$y, theta.star=theta.star, z=z.star)
cur.llps <- dPS.Rcpp(a=data$a, alpha=alpha, mid.points=mid.points,
                     bin.width=bin.width)

# tic <- proc.time()[3]
# for (i in 1:100) {
#   cur.llps <- dPS.Rcpp(a=data$a, alpha=alpha, mid.points=mid.points,
#                        bin.width=bin.width)
# }
# toc <- proc.time()[3]
# toc - tic
# MH adjustments
acc.alpha <- att.alpha <- mh.alpha <- 0.1

# storage
alpha.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  alpha.update <- updateAlpha(y=data$y, theta.star=theta.star, a=data$a,
                              alpha=alpha, cur.lly=cur.lly, cur.llps=cur.llps,
                              z=z.t, z.star=z.star, w=w.t, w.star=w.star,
                              mid.points=mid.points, bin.width=bin.width,
                              acc=acc.alpha, att=att.alpha, mh=mh.alpha)

  alpha      <- alpha.update$alpha
  w.star     <- alpha.update$w.star
  z.star     <- alpha.update$z.star
  theta.star <- alpha.update$theta.star
  cur.lly    <- alpha.update$cur.lly
  cur.llps   <- alpha.update$cur.llps
  att.alpha  <- alpha.update$att
  acc.alpha  <- alpha.update$acc

  alpha.keep[i] <- alpha

  if ((i %% 100) == 0) {
    cat("Iter", i, "\n")
    if (i > 3000) {
      start <- i - 3000
    } else {
      start <- 1
    }
    plot(alpha.keep[start:i], type="l",
         ylab=bquote(paste(alpha, " mh =", .(round(mh.alpha, 4)))),
         main=bquote(paste(alpha, " true =", .(alpha.t))))
  }

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.alpha, att=att.alpha, mh=mh.alpha)
    acc.alpha <- mh.update$acc
    att.alpha <- mh.update$att
    mh.alpha  <- mh.update$mh
  }
}

# test update for rho
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps   <- 20000
burn    <- 10000
xi.t    <- 0.1
beta.t  <- c(1, -1, 0)
alpha.t <- 0.8
rho.t   <- 3
thresh  <- 0
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t, thresh=thresh)
x.beta.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta.t[, t] <- x[, t, ] %*% beta.t
}
z.t <- getZ(xi=xi.t, x.beta=x.beta.t)

# initialization and prior distribution
rho <- 5
dw2 <- as.matrix(rdist(s, knots))^2
w   <- stdW(makeW(dw2=dw2, rho=rho))
w.star <- w^(1 / alpha.t)
z.star <- z.t^(1 / alpha.t)
theta.star <- getThetaStar(w.star=w.star, a=data$a)

cur.lly <- logLikeY(y=data$y, theta.star=theta.star, z.star=z.star)

# MH adjustments
acc.rho <- att.rho <- mh.rho <- 0.1

# storage
rho.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  rho.update <- updateRho(y=data$y, theta.star=theta.star, a=data$a,
                          alpha=alpha.t, cur.lly=cur.lly, z.star=z.star,
                          w=w, w.star=w.star, dw2=dw2,
                          rho=rho, rho.upper=15,
                          acc=acc.rho, att=att.rho, mh=mh.rho)

  rho        <- rho.update$rho
  w          <- rho.update$w
  w.star     <- rho.update$w.star
  theta.star <- rho.update$theta.star
  cur.lly    <- rho.update$cur.lly
  att.rho    <- rho.update$att
  acc.rho    <- rho.update$acc

  rho.keep[i] <- rho

  if ((i %% 500) == 0) {
    cat("Iter", i, "\n")
    if (i > 3000) {
      start <- i - 3000
    } else {
      start <- 1
    }
    plot(rho.keep[start:i], type="l",
         ylab=bquote(paste(rho, " mh =", .(round(mh.rho, 4)))),
         main=bquote(paste(rho, " true =", .(rho.t))))
  }

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh
  }
}

# test update for alpha and rho
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps   <- 20000
burn    <- 10000
xi.t    <- 0.1
beta.t  <- c(1, -1, 0)
alpha.t <- 0.8
rho.t   <- 1
thresh  <- 0
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t, thresh=thresh)
x.beta.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta.t[, t] <- x[, t, ] %*% beta.t
}
z.t <- getZ(xi=xi.t, x.beta=x.beta.t, thresh=thresh)

# initialization and prior distribution
npts <- 70
u.beta     <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
bin.width  <- u.beta[-1] - u.beta[-(npts + 1)]
alpha      <- 0.5
rho        <- 5
dw2        <- as.matrix(rdist(s, knots))^2
w          <- stdW(makeW(dw2=dw2, rho=rho))
w.star     <- w^(1 / alpha)
z.star     <- z.t^(1 / alpha)
theta.star <- getThetaStar(w.star=w.star, a=data$a)

cur.lly  <- logLikeY(y=data$y, theta.star=theta.star, z.star=z.star)
cur.llps <- dPS.Rcpp(a=data$a, alpha=alpha, mid.points=mid.points,
                     bin.width=bin.width)
# cur.llps.2 <- matrix(NA, nknots, nt)
# for (t in 1:nt) {
#   for (k in 1:nknots) {
#     cur.llps.2[k, t] <- dPS(a=data$a[k, t], alpha=alpha,
#                           mid.points=mid.points, bin.width=bin.width)
#   }
# }

# MH adjustments
acc.alpha <- att.alpha <- mh.alpha <- 0.1
acc.rho   <- att.rho   <- mh.rho   <- 0.1

# storage
alpha.keep <- rep(NA, nreps)
rho.keep   <- rep(NA, nreps)

for (i in 1:nreps) {
  alpha.update <- updateAlpha(y=data$y, theta.star=theta.star, a=data$a,
                              alpha=alpha, cur.lly=cur.lly, cur.llps=cur.llps,
                              z=z.t, z.star=z.star, w=w, w.star=w.star,
                              mid.points=mid.points, bin.width=bin.width,
                              acc=acc.alpha, att=att.alpha, mh=mh.alpha)

  alpha      <- alpha.update$alpha
  w.star     <- alpha.update$w.star
  z.star     <- alpha.update$z.star
  theta.star <- alpha.update$theta.star
  cur.lly    <- alpha.update$cur.lly
  cur.llps   <- alpha.update$cur.llps
  att.alpha  <- alpha.update$att
  acc.alpha  <- alpha.update$acc

  alpha.keep[i] <- alpha

  rho.update <- updateRho(y=data$y, theta.star=theta.star, a=data$a,
                          alpha=alpha, cur.lly=cur.lly, z.star=z.star,
                          w=w, w.star=w.star,
                          dw2=dw2, rho=rho, rho.upper=15,
                          acc=acc.rho, att=att.rho, mh=mh.rho)

  rho        <- rho.update$rho
  w          <- rho.update$w
  w.star     <- rho.update$w.star
  theta.star <- rho.update$theta.star
  cur.lly    <- rho.update$cur.lly
  att.rho    <- rho.update$att
  acc.rho    <- rho.update$acc

  rho.keep[i] <- rho

  if ((i %% 100) == 0) {
    cat("Iter", i, "\n")
    if (i > 3000) {
      start <- i - 3000
    } else {
      start <- 1
    }
    par(mfrow=c(1, 2))
    plot(rho.keep[start:i], type="l",
         ylab=bquote(paste(rho, " mh =", .(round(mh.rho, 4)))),
         main=bquote(paste(rho, " true =", .(rho.t))))
    plot(alpha.keep[start:i], type="l",
         ylab=bquote(paste(alpha, " mh =", .(round(mh.alpha, 4)))),
         main=bquote(paste(alpha, " true =", .(alpha.t))))
  }

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.alpha, att=att.alpha, mh=mh.alpha)
    acc.alpha <- mh.update$acc
    att.alpha <- mh.update$att
    mh.alpha  <- mh.update$mh

    mh.update <- mhUpdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho <- mh.update$acc
    att.rho <- mh.update$att
    mh.rho  <- mh.update$mh
  }
}

# test update for alpha and a
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps   <- 20000
burn    <- 10000
xi.t    <- 0.1
beta.t  <- c(1, -1, 0)
alpha.t <- 0.8
rho.t   <- 1
thresh  <- 0
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t, thresh=thresh)
x.beta.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta.t[, t] <- x[, t, ] %*% beta.t
}
z.t <- getZ(xi=xi.t, x.beta=x.beta.t, thresh=thresh)
z.star.t <- z.t^(1 / alpha.t)

# initialization and prior distribution
npts <- 70
u.beta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
bin.width <- u.beta[-1] - u.beta[-(npts + 1)]
dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, rho=rho.t))
w.star.t <- w.t^(1 / alpha.t)
a <- matrix(5, nknots, nt)
theta.star <- getThetaStar(w.star=w.star.t, a=a)
theta.star.t <- getThetaStar(w.star=w.star.t, a=data$a)

cur.lly    <- logLikeY(y=data$y, theta.star=theta.star, z.star=z.star.t)
cur.lly.t  <- logLikeY(y=data$y, theta.star=theta.star.t, z.star=z.star.t)
cur.llps   <- dPS.Rcpp(a, alpha=alpha.t, mid.points=mid.points,
                     bin.width=bin.width)
cur.llps.t <- dPS.Rcpp(a=data$a, alpha=alpha.t, mid.points=mid.points,
                       bin.width=bin.width)

# MH adjustments
cuts  <- exp(c(-1, 0, 1, 2, 5, 10))
mh.a  <- rep(0.5, 100)
acc.a <- att.a <- 0 * mh.a

# storage
a.keep     <- array(NA, dim=c(nreps, nknots, nt))

# Rprof(filename = "Rprof.out",line.profiling = TRUE)
for (i in 1:nreps) {
  old.a <- a
  a.update <- updateA(y=data$y, theta.star=theta.star, a=a, alpha=alpha.t,
                      cur.lly=cur.lly, cur.llps=cur.llps, z.star=z.star.t,
                      w.star=w.star.t, mid.points=mid.points,
                      bin.width=bin.width, mh=mh.a, cuts=cuts)

  a          <- a.update$a
  theta.star <- a.update$theta.star
  cur.lly    <- a.update$cur.lly
  cur.llps   <- a.update$cur.llps

  a.keep[i, , ] <- a

  if ((i %% 50) == 0) {
    par(mfrow=c(3, 4))
    cat("Iter", i, "\n")
    if (i > 3000) {
      start <- i - 3000
    } else {
      start <- 1
    }
    plot(a.keep[start:i, 1, 1], type="l",
         main=paste("a 1, 1 true:", round(data$a[1, 1], 3)))
    plot(a.keep[start:i, 3, 1], type="l",
         main=paste("a 3, 1 true:", round(data$a[3, 1], 3)))
    plot(a.keep[start:i, 5, 1], type="l",
         main=paste("a 5, 1 true:", round(data$a[5, 1], 3)))
    plot(a.keep[start:i, 7, 1], type="l",
         main=paste("a 7, 1 true:", round(data$a[7, 1], 3)))
    plot(a.keep[start:i, 1, 5], type="l",
         main=paste("a 1, 5 true:", round(data$a[1, 5], 3)))
    plot(a.keep[start:i, 3, 5], type="l",
         main=paste("a 3, 5 true:", round(data$a[3, 5], 3)))
    plot(a.keep[start:i, 5, 5], type="l",
         main=paste("a 5, 5 true:", round(data$a[5, 5], 3)))
    plot(a.keep[start:i, 7, 5], type="l",
         main=paste("a 7, 5 true:", round(data$a[7, 5], 3)))
    plot(a.keep[start:i, 1, 10], type="l",
         main=paste("a 1, 10 true:", round(data$a[1, 10], 3)))
    plot(a.keep[start:i, 3, 10], type="l",
         main=paste("a 3, 10 true:", round(data$a[3, 10], 3)))
    plot(a.keep[start:i, 5, 10], type="l",
         main=paste("a 5, 10 true:", round(data$a[5, 10], 3)))
    plot(a.keep[start:i, 7, 10], type="l",
         main=paste("a 7, 10 true:", round(data$a[7, 10], 3)))
  }
}

# test update for a
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps   <- 20000
burn    <- 10000
xi.t    <- 0.1
beta.t  <- c(1, -1, 0)
alpha.t <- 0.8
rho.t   <- 1
thresh  <- 0
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t, thresh=thresh)
x.beta.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta.t[, t] <- x[, t, ] %*% beta.t
}
z.t <- getZ(xi=xi.t, x.beta=x.beta.t, thresh=thresh)
z.star.t <- z.t^(1 / alpha.t)

# initialization and prior distribution
npts <- 70
u.beta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
bin.width <- u.beta[-1] - u.beta[-(npts + 1)]
dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, rho=rho.t))
w.star.t <- w.t^(1 / alpha.t)
a <- matrix(5, nknots, nt)
theta.star <- getThetaStar(w.star=w.star.t, a=a)
theta.star.t <- getThetaStar(w.star=w.star.t, a=data$a)

cur.lly    <- logLikeY(y=data$y, theta.star=theta.star, z.star=z.star.t)
cur.lly.t  <- logLikeY(y=data$y, theta.star=theta.star.t, z.star=z.star.t)
cur.llps   <- dPS.Rcpp(a, alpha=alpha.t, mid.points=mid.points,
                     bin.width=bin.width)
cur.llps.t <- dPS.Rcpp(a=data$a, alpha=alpha.t, mid.points=mid.points,
                       bin.width=bin.width)

# MH adjustments
cuts  <- exp(c(-1, 0, 1, 2, 5, 10))
mh.a  <- rep(0.5, 100)
acc.a <- att.a <- 0 * mh.a

# storage
a.keep     <- array(NA, dim=c(nreps, nknots, nt))

# Rprof(filename = "Rprof.out",line.profiling = TRUE)
for (i in 1:nreps) {
  old.a <- a
  a.update <- updateA(y=data$y, theta.star=theta.star, a=a, alpha=alpha.t,
                      cur.lly=cur.lly, cur.llps=cur.llps, z.star=z.star.t,
                      w.star=w.star.t, mid.points=mid.points,
                      bin.width=bin.width, mh=mh.a, cuts=cuts)

  a          <- a.update$a
  theta.star <- a.update$theta.star
  cur.lly    <- a.update$cur.lly
  cur.llps   <- a.update$cur.llps

  a.keep[i, , ] <- a

  if ((i %% 50) == 0) {
    par(mfrow=c(3, 4))
    cat("Iter", i, "\n")
    if (i > 3000) {
      start <- i - 3000
    } else {
      start <- 1
    }
    plot(a.keep[start:i, 1, 1], type="l",
         main=paste("a 1, 1 true:", round(data$a[1, 1], 3)))
    plot(a.keep[start:i, 3, 1], type="l",
         main=paste("a 3, 1 true:", round(data$a[3, 1], 3)))
    plot(a.keep[start:i, 5, 1], type="l",
         main=paste("a 5, 1 true:", round(data$a[5, 1], 3)))
    plot(a.keep[start:i, 7, 1], type="l",
         main=paste("a 7, 1 true:", round(data$a[7, 1], 3)))
    plot(a.keep[start:i, 1, 5], type="l",
         main=paste("a 1, 5 true:", round(data$a[1, 5], 3)))
    plot(a.keep[start:i, 3, 5], type="l",
         main=paste("a 3, 5 true:", round(data$a[3, 5], 3)))
    plot(a.keep[start:i, 5, 5], type="l",
         main=paste("a 5, 5 true:", round(data$a[5, 5], 3)))
    plot(a.keep[start:i, 7, 5], type="l",
         main=paste("a 7, 5 true:", round(data$a[7, 5], 3)))
    plot(a.keep[start:i, 1, 10], type="l",
         main=paste("a 1, 10 true:", round(data$a[1, 10], 3)))
    plot(a.keep[start:i, 3, 10], type="l",
         main=paste("a 3, 10 true:", round(data$a[3, 10], 3)))
    plot(a.keep[start:i, 5, 10], type="l",
         main=paste("a 5, 10 true:", round(data$a[5, 10], 3)))
    plot(a.keep[start:i, 7, 10], type="l",
         main=paste("a 7, 10 true:", round(data$a[7, 10], 3)))
  }
}




# test update for alpha and a
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 20000
burn  <- 10000
xi.t <- 0.1
beta.t <- c(2, 0, 0)
alpha.t <- 0.8
rho.t   <- 1
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t)

# initialization and prior distribution
alpha <- 0.5
npts <- 50
u.beta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
bin.width <- u.beta[-1] - u.beta[-(npts + 1)]
x.beta <- matrix(2, ns, nt)
z.t <- getZ(xi=xi.t, x.beta=x.beta)
dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, rho=rho.t))
w.star <- w.t^(1 / alpha)
a <- matrix(1, nknots, nt)
theta.star <- getThetaStar(w.star=w.star, a=a)

cur.lly <- logLikeY(y=data$y, theta.star=theta.star, alpha=alpha, z=z.t)
cur.llps <- matrix(NA, nknots, nt)
for (t in 1:nt) {
  for (k in 1:nknots) {
    cur.llps[k, t] <- dPS(a=a[k, t], alpha=alpha,
                          mid.points=mid.points, bin.width=bin.width)
  }
}

# MH adjustments
cuts      <- exp(c(-1, 0, 1, 2, 5, 10))
mh.a      <- rep(1, 100)
acc.a     <- att.a     <- 0 * mh.a
acc.alpha <- att.alpha <- mh.alpha <- 1

# storage
a.keep     <- array(NA, dim=c(nreps, nknots, nt))
alpha.keep <- rep(NA, nreps)

Rprof("rprof.out", line.profiling=TRUE)
for (i in 1:2000) {
  old.a <- a
  a.update <- updateA(y=data$y, theta.star=theta.star, a=a, alpha=alpha,
                      cur.lly=cur.lly, cur.llps=cur.llps, z=z.t,
                      w.star=w.star.t, mid.points=mid.points,
                      bin.width=bin.width, mh=mh.a, cuts=cuts)

  a          <- a.update$a
  theta.star <- a.update$theta.star
  cur.lly    <- a.update$cur.lly
  cur.llps   <- a.update$cur.llps

  a.keep[i, , ] <- a

  alpha.update <- updateAlpha(y=data$y, theta.star=theta.star, a=a,
                              alpha=alpha, cur.lly=cur.lly, cur.llps=cur.llps,
                              z=z.t, w.star=w.star.t,
                              mid.points=mid.points, bin.width=bin.width,
                              acc=acc.alpha, att=att.alpha, mh=mh.alpha)

  alpha      <- alpha.update$alpha
  theta.star <- alpha.update$theta.star
  cur.lly    <- alpha.update$cur.lly
  cur.llps   <- alpha.update$cur.llps
  att.alpha  <- alpha.update$att
  acc.alpha  <- alpha.update$acc

  alpha.keep[i] <- alpha

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.alpha, att=att.alpha, mh=mh.alpha)
    acc.alpha <- mh.update$acc
    att.alpha <- mh.update$att
    mh.alpha  <- mh.update$mh
  }

  if ((i %% 100) == 0) {
    cat("Iter", i, "\n")
    if (i > 3000) {
      start <- i - 3000
    } else {
      start <- 1
    }
    par(mfrow=c(3, 4))
    plot(alpha.keep[start:i], type="l",
         main=bquote(paste(alpha, " mh =", .(round(mh.alpha, 4)))))
    plot(a.keep[start:i, 1, 1], type="l",
         main=paste("a 1, 1 true:", round(data$a[1, 1], 2)))
    plot(a.keep[start:i, 3, 1], type="l",
         main=paste("a 10, 1 true:", round(data$a[3, 1], 2)))
    plot(a.keep[start:i, 5, 1], type="l",
         main=paste("a 20, 1 true:", round(data$a[5, 1], 2)))
    plot(a.keep[start:i, 7, 1], type="l",
         main=paste("a 30, 1 true:", round(data$a[5, 1], 2)))
    plot(a.keep[start:i, 1, 5], type="l",
         main=paste("a 1, 1 true:", round(data$a[1, 5], 2)))
    plot(a.keep[start:i, 3, 5], type="l",
         main=paste("a 10, 1 true:", round(data$a[3, 5], 2)))
    plot(a.keep[start:i, 5, 5], type="l",
         main=paste("a 20, 1 true:", round(data$a[5, 5], 2)))
    plot(a.keep[start:i, 7, 5], type="l",
         main=paste("a 30, 1 true:", round(data$a[7, 5], 2)))
    plot(a.keep[start:i, 1, 10], type="l",
         main=paste("a 1, 1 true:", round(data$a[1, 10], 2)))
    plot(a.keep[start:i, 3, 10], type="l",
         main=paste("a 10, 1 true:", round(data$a[3, 10], 2)))
    plot(a.keep[start:i, 5, 10], type="l",
         main=paste("a 20, 1 true:", round(data$a[5, 10], 2)))
  }


}
Rprof(NULL)
summaryRprof(filename="rprof.out", lines="show")

# test update for rho
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 10000
burn  <- 7000
xi.t <- 0.1
beta.t <- c(2, 0, 0)
alpha.t <- 0.2
rho.t   <- 1
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t)

# initialization and prior distribution
rho <- 5
x.beta <- matrix(2, ns, nt)
z.t <- getZ(xi=xi.t, x.beta=x.beta)
dw2 <- as.matrix(rdist(s, knots))^2
w <- stdW(makeW(dw2=dw2, rho=rho))
w.star <- w^(1 / alpha.t)
theta.star <- getThetaStar(w.star=w.star, a=data$a)
cur.lly <- logLikeY(y=data$y, theta.star=theta.star, alpha=alpha.t, z=z.t)

# MH adjustments
acc.rho <- att.rho <- mh.rho <- 0.1

# storage
rho.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  rho.update <- updateRho(y=data$y, theta.star=theta.star, a=data$a,
                          alpha=alpha.t, cur.lly=cur.lly, z=z.t, w=w,
                          w.star=w.star, dw2=dw2, rho=rho, rho.upper=15,
                          acc=acc.rho, att=att.rho, mh=mh.rho)
  rho        <- rho.update$rho
  w          <- rho.update$w
  w.star     <- rho.update$w.star
  theta.star <- rho.update$theta.star
  cur.lly    <- rho.update$cur.lly
  att.rho    <- rho.update$att
  acc.rho    <- rho.update$acc

  rho.keep[i] <- rho

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.rho, att=att.rho, mh=mh.rho)
    acc.rho   <- mh.update$acc
    att.rho   <- mh.update$att
    mh.rho    <- mh.update$mh
  }
  if (i %% 500 == 0) {
    if (i < 3000) {
      start <- 1
    } else {
      start <- i - 3000
    }
    plot(rho.keep[start:i], type="l")
    cat("iter ", i, "\n")
  }
}
plot(rho.keep, type="l")


# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent=FALSE) {
    vars <- codetools::findGlobals(f)
    found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
    names <- names(found)[found]

    if ((length(names) > 0)) {
      sum.nfncs <- 0
      for (i in 1:length(names)) {
        if(!is.function(eval(parse(text=names[i])))) {sum.nfncs <- sum.nfncs + 1}
      }
      if (sum.nfncs > 0) {
        warning("global variables used: ", paste(names(found)[found], collapse=', '))
        return(invisible(FALSE))
      }
    }

    !any(found)
}

checkStrict(updateBeta)
checkStrict(updateXi)
checkStrict(updateA)
checkStrict(updateAlpha)
checkStrict(updateRho)
checkStrict(getZ)
checkStrict(getThetaStar)
checkStrict(makeW)
checkStrict(stdW)
checkStrict(logLikeY)
checkStrict(dPS)
checkStrict(ld)
checkStrict(rPS)
checkStrict(rRareBinaryInd)
checkStrict(rRareBinarySpat)
checkStrict(mhUpdate)
checkStrict(ld2)
checkStrict(dlognormal)
checkStrict(get.level)
checkStrict(mcmc)

# timing tests with ijulia
# Looking at timing
beta <- rep(0, 3)
beta.m <- rep(0, 3)
beta.s <- rep(0, 3)
alpha <- 0.5
y <- matrix(rbinom(ns * nt, size=1, prob=0.5), ns, nt)
theta.star <- matrix(rchisq(ns * nt, df=3), ns, nt)
z <- matrix(rchisq(ns * nt, df=2), ns, nt)
x <- array(runif(ns * nt * 3), dim=c(ns, nt, 3))
x[, , 1] <- 1
x.beta <- matrix(rnorm(ns * nt), ns, nt)
nreps <- 5000

cur.lly <- logLikeY(y=y, theta.star=theta.star, alpha=alpha, z=z)
tic <- proc.time()[3]
for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta.star=theta.star, alpha=alpha, z=z,
                            beta=beta, beta.m=beta.m, beta.s=beta.s,
                            xi=xi.t, x=x, cur.lly=cur.lly,
                            acc=acc.beta, att=att.beta, mh=mh.beta)
  beta     <- beta.update$beta
  x.beta   <- beta.update$x.beta
  z        <- beta.update$z
  cur.lly  <- beta.update$cur.lly
  att.beta <- beta.update$att
  acc.beta <- beta.update$acc

  beta.keep[i, ] <- beta
  if (i %% 500 == 0) {
    print(i)
  }

  # if (i < burn / 2) {
  #   mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta)
  #   acc.beta  <- mh.update$acc
  #   att.beta  <- mh.update$att
  #   mh.beta   <- mh.update$mh
  # }
}
toc <- proc.time()[3]
toc - tic


can.llps <- matrix(0.0, nknots, nt)
tic <- proc.time()
  for(t in 1:nt) {
    for (k in 1:nknots) {
      can.llps[k, t] <- dPS(data$a[k, t], alpha, mid.points, bin.width)
    }
  }
toc <- proc.time()

# timing updates agains julia
tic <- proc.time()
for (i in 1:50) {
  alpha.update <- updateAlpha(y=data$y, theta.star=theta.star, a=data$a,
                              alpha=alpha, cur.lly=cur.lly, cur.llps=cur.llps,
                              z=z.t, w.star=w.star.t,
                              mid.points=mid.points, bin.width=bin.width,
                              acc=acc.alpha, att=att.alpha, mh=mh.alpha)

  alpha      <- alpha.update$alpha
  theta.star <- alpha.update$theta.star
  cur.lly    <- alpha.update$cur.lly
  cur.llps   <- alpha.update$cur.llps
  att.alpha  <- alpha.update$att
  acc.alpha  <- alpha.update$acc

  alpha.keep[i] <- alpha

  # if ((i %% 100) == 0) {
  #   cat("Iter", i, "\n")
  #   if (i > 3000) {
  #     start <- i - 3000
  #   } else {
  #     start <- 1
  #   }
  #   plot(alpha.keep[start:i], type="l",
  #        ylab=bquote(paste(alpha, " mh =", .(round(mh.alpha, 4)))),
  #        main=bquote(paste(alpha, " true =", .(alpha.t))))
  # }

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.alpha, att=att.alpha, mh=mh.alpha)
    acc.alpha <- mh.update$acc
    att.alpha <- mh.update$att
    mh.alpha  <- mh.update$mh
  }
}
toc <- proc.time()
toc - tic