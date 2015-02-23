rm(list=ls())
options(warn=2)
library(fields)
source("auxfunctions.R")
source("updateModel.R")
set.seed(10)
ns   <- 50
nt   <- 20
s    <- cbind(runif(ns), runif(ns))
x <- array(1, dim=c(ns, nt, 3))
x[, , 2] <- s[, 1]
x[, , 3] <- s[, 2]
knots.1 <- seq(0, 1, length=10)
knots.2 <- seq(0, 1, length=10)
knots <- expand.grid(knots.1, knots.2)

# test update for beta
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 20000
burn  <- 10000
xi.t <- 0.1
beta.t <- c(2, 0, 0)
theta.t <- 1
alpha.t <- 1
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t)

# initialization and prior distribution
beta <- c(0, 0, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(1, ns, nt)
z <- getZ(xi=xi.t, x.beta=x.beta)
cur.lly <- logLikeY(y=y, theta=theta.t, alpha=alpha.t, z=z)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(1, 3)

# storage
beta.keep <- matrix(NA, nreps, 3)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta=theta.t, alpha=alpha.t, z=z,
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


  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta)
    acc.beta  <- mh.update$acc
    att.beta  <- mh.update$att
    mh.beta   <- mh.update$mh
  }
}

par(mfrow=c(1, 3))
plot(beta.keep[, 1], type="l", main="beta 0")
plot(beta.keep[, 2], type="l", main="beta 1")
plot(beta.keep[, 3], type="l", main="beta 2")

# test update for xi
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 20000
burn  <- 10000
xi.t <- 0.3
beta.t <- c(2, 0, 0)
theta.t <- 1
alpha.t <- 1
x.beta.t <- matrix(2, ns, nt)
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t)

# initialization and prior distribution
xi <- 0
xi.m <- 0
xi.s <- 10
z <- (1 + xi * x.beta.t)^(1 / xi)
cur.lly <- logLikeY(y=y, theta=theta.t, alpha=alpha.t, z=z)

# MH adjustments
acc.xi <- att.xi <- mh.xi <- rep(1, 3)

# storage
xi.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  xi.update <- updateXi(y=y, theta=theta.t, alpha=alpha.t, z=z, x.beta=x.beta.t,
                        xi=xi, xi.m=xi.m, xi.s=xi.s,
                        cur.lly=cur.lly,
                        acc=acc.xi, att=att.xi, mh=mh.xi)
  xi      <- xi.update$xi
  z       <- xi.update$z
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
}

plot(xi.keep, type="l")

# test update for xi and beta (just intercept)
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 20000
burn  <- 10000
xi.t <- 0.3
beta.t <- c(2, 0, 0)
theta.t <- 1
alpha.t <- 1
x.beta.t <- matrix(2, ns, nt)
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t)

# initialization and prior distribution
beta <- c(0, 0, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(1, ns, nt)
xi <- 0
xi.m <- 0
xi.s <- 10
z <- getZ(xi=xi, x.beta=x.beta)
cur.lly <- logLikeY(y=y, theta=theta.t, alpha=alpha.t, z=z)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(1, 3)
acc.xi <- att.xi <- mh.xi <- rep(1, 3)

# storage
beta.keep <- matrix(NA, nreps, 3)
xi.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta=theta.t, alpha=alpha.t, z=z,
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


  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta)
    acc.beta  <- mh.update$acc
    att.beta  <- mh.update$att
    mh.beta   <- mh.update$mh
  }

  xi.update <- updateXi(y=y, theta=theta.t, alpha=alpha.t, z=z, x.beta=x.beta,
                        xi=xi, xi.m=xi.m, xi.s=xi.s,
                        cur.lly=cur.lly,
                        acc=acc.xi, att=att.xi, mh=mh.xi)
  xi      <- xi.update$xi
  z       <- xi.update$z
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
}

par(mfrow=c(2, 2))
plot(beta.keep[, 1], type="l", main="beta 0")
plot(beta.keep[, 2], type="l", main="beta 1")
plot(beta.keep[, 3], type="l", main="beta 2")
plot(xi.keep, type="l")

# test update for xi and beta (intercept and x)
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 20000
burn  <- 10000
xi.t <- 0.1
beta.t <- c(2, 1, 0)
theta.t <- 1
alpha.t <- 1
x.beta.t <- matrix(2, ns, nt)
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t)

# initialization and prior distribution
beta <- c(0, 0, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(1, ns, nt)
xi <- 0
xi.m <- 0
xi.s <- 10
z <- getZ(xi=xi, x.beta=x.beta)
cur.lly <- logLikeY(y=y, theta=theta.t, alpha=alpha.t, z=z)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(1, 3)
acc.xi <- att.xi <- mh.xi <- rep(1, 3)

# storage
beta.keep <- matrix(NA, nreps, 3)
xi.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta=theta.t, alpha=alpha.t, z=z,
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


  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta)
    acc.beta  <- mh.update$acc
    att.beta  <- mh.update$att
    mh.beta   <- mh.update$mh
  }

  xi.update <- updateXi(y=y, theta=theta.t, alpha=alpha.t, z=z, x.beta=x.beta,
                        xi=xi, xi.m=xi.m, xi.s=xi.s,
                        cur.lly=cur.lly,
                        acc=acc.xi, att=att.xi, mh=mh.xi)
  xi      <- xi.update$xi
  z       <- xi.update$z
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
}

par(mfrow=c(2, 2))
plot(beta.keep[, 1], type="l", main="beta 0")
plot(beta.keep[, 2], type="l", main="beta 1")
plot(beta.keep[, 3], type="l", main="beta 2")
plot(xi.keep, type="l")

# test update for alpha
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 20000
burn  <- 10000
xi.t <- 0.1
beta.t <- c(2, 0, 0)
alpha.t <- 0.8
rho.t   <- 0.1
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t)

# initialization and prior distribution
x.beta <- matrix(2, ns, nt)
z.t <- getZ(xi=xi.t, x.beta=x.beta)
dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, logrho=log(rho.t)))
alpha <- 0.5
theta <- getTheta(w=w.t, a=data$a, alpha=alpha)

cur.lly <- logLikeY(y=y, theta=theta.t, alpha=alpha, z=z)
cur.llps <- matrix(NA, nknots, nt)
for (t in 1:nt) {
  for (k in 1:nknots) {
    cur.llps[k, t] <- dPS(a, alpha)
  }
}

# MH adjustments
acc.alpha <- att.alpha <- mh.alpha <- 1

# storage
alpha.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  alpha.update <- updateAlpha(y=y, theta=theta, a=data$a, alpha=alpha,
                              cur.lly=cur.lly, cur.llps=cur.llps, z=z.t, w=w.t,
                              acc=acc.alpha, att=att.alpha, mh=mh.alpha)

  alpha     <- alpha.update$alpha
  theta     <- alpha.update$theta
  cur.lly   <- alpha.update$cur.lly
  cur.llps  <- alpha.update$cur.llps
  att.alpha <- alpha.update$att
  acc.alpha <- alpha.update$acc

  alpha.keep[i] <- alpha

  if (i < burn / 2) {
    mh.update <- mhUpdate(acc=acc.alpha, att=att.alpha, mh=mh.alpha)
    acc.alpha <- mh.update$acc
    att.alpha <- mh.update$att
    mh.alpha  <- mh.update$mh
  }
}

