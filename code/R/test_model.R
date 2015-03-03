# model update testing
# independence:
#   beta: PASS
#   xi: PASS
#   beta and xi: PASS
# spatial dependence:
#   beta: PASS
#   xi: 
#   beta and xi:
#   alpha:
#   rho:
#   alpha and rho:
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
source("auxfunctions.R")
source("updateModel.R")
set.seed(10)
ns   <- 200
nt   <- 40
s    <- cbind(runif(ns, 0, 10), runif(ns, 0, 10))
x <- array(1, dim=c(ns, nt, 3))
x[, , 2] <- s[, 1]
x[, , 3] <- s[, 2]
knots.1 <- seq(1, 10, length=9)
knots.2 <- seq(1, 10, length=9)
knots <- expand.grid(knots.1, knots.2)
nknots <- nrow(knots)

# test update for beta
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 5000
burn  <- 2000
xi.t <- 0.1
beta.t <- c(2, 0, 0)
theta.star.t <- 1
alpha.t <- 1
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t)

# initialization and prior distribution
beta <- c(0, 0, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(1, ns, nt)
z <- getZ(xi=xi.t, x.beta=x.beta)
cur.lly <- logLikeY(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(0.1, 3)

# storage
beta.keep <- matrix(NA, nreps, 3)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z,
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

# test update for beta
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 10000
burn  <- 7000
xi.t <- 0.1
beta.t <- c(1, -1, 0)
theta.star.t <- 1
alpha.t <- 1
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t)

# initialization and prior distribution
beta <- c(0, 0, 0)
beta.m <- 0
beta.s <- 5
x.beta <- matrix(0, ns, nt)
z <- getZ(xi=xi.t, x.beta=x.beta)
cur.lly <- logLikeY(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(0.1, 3)

# storage
beta.keep <- matrix(NA, nreps, 3)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z,
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
nreps <- 10000
burn  <- 7000
xi.t <- 0.1
beta.t <- c(1, -1, 0)
theta.star.t <- 1
alpha.t <- 1
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t)
x.beta.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta.t[, t] <- x[, t, ] %*% beta.t
}

# initialization and prior distribution
xi <- 0
xi.m <- 0
xi.s <- 10
z <- getZ(xi=xi, x.beta=x.beta.t)
cur.lly <- logLikeY(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z)

# MH adjustments
acc.xi <- att.xi <- mh.xi <- rep(1, 3)

# storage
xi.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  xi.update <- updateXi(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z,
                        x.beta=x.beta.t, xi=xi, xi.m=xi.m, xi.s=xi.s,
                        cur.lly=cur.lly,
                        acc=acc.xi, att=att.xi, mh=mh.xi)
  xi      <- xi.update$xi
  z       <- xi.update$z
  cur.lly <- xi.update$cur.lly
  att.xi  <- xi.update$att
  acc.xi  <- xi.update$acc

  xi.keep[i] <- xi
  if (i %% 500 == 0) {
    print(i)
  }

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
nreps <- 7000
burn  <- 5000
xi.t <- 0.1
beta.t <- c(1, 0, 0)
theta.star.t <- 1
alpha.t <- 1
y <- rRareBinaryInd(x, beta=beta.t, xi=xi.t)

# initialization and prior distribution
beta <- c(0, 0, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(0, ns, nt)
xi <- 0.1
xi.m <- 0
xi.s <- 1
z <- getZ(xi=xi, x.beta=x.beta)
cur.lly <- logLikeY(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(0.1, 3)
acc.xi <- att.xi <- mh.xi <- rep(1, 3)

# storage
beta.keep <- matrix(NA, nreps, 3)
xi.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z,
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

  xi.update <- updateXi(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z,
                        x.beta=x.beta, xi=xi, xi.m=xi.m, xi.s=xi.s,
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
    plot(xi.keep[start:i], type="l")
  }
}

# test update for xi and beta
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 15000
burn  <- 12000
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
cur.lly <- logLikeY(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(0.01, 3)
acc.xi <- att.xi <- mh.xi <- 0.1

# storage
beta.keep <- matrix(NA, nreps, 3)
xi.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z,
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

  xi.update <- updateXi(y=y, theta.star=theta.star.t, alpha=alpha.t, z=z,
                        x.beta=x.beta, xi=xi, xi.m=xi.m, xi.s=xi.s,
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
# seems to be working
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 20000
burn  <- 10000
xi.t <- 0.1
beta.t <- c(1, -1, 0)
alpha.t <- 0.7
rho.t   <- 1
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t)

dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, rho=rho.t))
theta.star.t <- getThetaStar(w=w.t, a=data$a, alpha=alpha.t)

# initialization and prior distribution
beta <- c(0, 0, 0)
beta.m <- 0
beta.s <- 10
x.beta <- matrix(0, ns, nt)
z <- getZ(xi=xi.t, x.beta=x.beta)
cur.lly <- logLikeY(y=data$y, theta.star=theta.star.t, alpha=alpha.t, z=z)

# MH adjustments
acc.beta <- att.beta <- mh.beta <- rep(0.1, 3)

# storage
beta.keep <- matrix(NA, nreps, 3)

for (i in 1:nreps) {
  beta.update <- updateBeta(y=data$y, theta.star=theta.star.t, alpha=alpha.t,
                            z=z, beta=beta, beta.m=beta.m, beta.s=beta.s,
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













#### there's a problem in the updates for the dependence structure
#### rho and alpha both are estimated poorly
#### as alpha -> 1, it does better, but anywhere below 0.8, alpha starts to
#### come back way too low.

#### need to go back through the likelihood function and all other functions
#### that are associated with alpha to make sure that they are all correct

# test update for alpha
set.seed(15)
source("auxfunctions.R")
source("updateModel.R")
nreps <- 20000
burn  <- 10000
xi.t <- 0.1
beta.t <- c(1, -1, 0)
alpha.t <- 0.7
rho.t   <- 1
data <- rRareBinarySpat(x, s=s, knots=knots, beta=beta.t,
                        xi=xi.t, alpha=alpha.t, rho=rho.t)
x.beta.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  x.beta.t[, t] <- x[, t, ] %*% beta.t
}
z.t <- getZ(xi=xi.t, x.beta=x.beta.t)

# initialization and prior distribution
npts <- 70
u.beta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
bin.width <- u.beta[-1] - u.beta[-(npts + 1)]
dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, rho=rho.t))
alpha <- 0.7
theta.star <- getThetaStar(w=w.t, a=data$a, alpha=alpha)

cur.lly <- logLikeY(y=data$y, theta.star=theta.star, alpha=alpha, z=z.t)
cur.llps <- matrix(NA, nknots, nt)
for (t in 1:nt) {
  for (k in 1:nknots) {
    cur.llps[k, t] <- dPS(a=data$a[k, t], alpha=alpha,
                          mid.points=mid.points, bin.width=bin.width)
  }
}

# MH adjustments
acc.alpha <- att.alpha <- mh.alpha <- 0.01

# storage
alpha.keep <- rep(NA, nreps)

for (i in 1:5000) {
  alpha.update <- updateAlpha(y=data$y, theta.star=theta.star, a=data$a,
                              alpha=alpha, cur.lly=cur.lly, cur.llps=cur.llps,
                              z=z.t, w=w.t,
                              mid.points=mid.points, bin.width=bin.width,
                              acc=acc.alpha, att=att.alpha, mh=mh.alpha)

  alpha      <- alpha.update$alpha
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



















# test update for a
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
npts <- 100
u.beta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
bin.width <- u.beta[-1] - u.beta[-(npts + 1)]
x.beta <- matrix(2, ns, nt)
z.t <- getZ(xi=xi.t, x.beta=x.beta)
dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, rho=rho.t))
a <- matrix(5, nknots, nt)
theta.star <- getThetaStar(w=w.t, a=a, alpha=alpha.t)
theta.star.t <- getThetaStar(w=w.t, a=data$a, alpha=alpha.t)

cur.lly <- logLikeY(y=data$y, theta.star=theta.star, alpha=alpha.t, z=z.t)
cur.lly.t <- logLikeY(y=data$y, theta.star=theta.star.t, alpha=alpha.t, z=z.t)
cur.llps <- matrix(NA, nknots, nt)
for (t in 1:nt) {
  for (k in 1:nknots) {
    cur.llps[k, t] <- dPS(a=a[k, t], alpha=alpha.t,
                          mid.points=mid.points, bin.width=bin.width)
  }
}

cur.llps.t <- matrix(NA, nknots, nt)
for (t in 1:nt) {
  for (k in 1:nknots) {
    cur.llps.t[k, t] <- dPS(a=data$a[k, t], alpha=alpha.t,
                          mid.points=mid.points, bin.width=bin.width)
  }
}

# MH adjustments
cuts  <- exp(c(-1, 0, 1, 2, 5, 10))
mh.a  <- rep(0.5, 100)
acc.a <- att.a <- 0 * mh.a

# storage
a.keep     <- array(NA, dim=c(nreps, nknots, nt))

for (i in 1:2000) {
  old.a <- a
  a.update <- updateA(y=data$y, theta.star=theta.star, a=a, alpha=alpha.t,
                      cur.lly=cur.lly, cur.llps=cur.llps, z=z.t, w=w.t,
                      mid.points=mid.points, bin.width=bin.width,
                      mh=mh.a, cuts=cuts)
  level <- get.level(old.a, cuts)
  for (j in 1:length(mh.a)) {
    acc.a[j] <- acc.a[j] + sum(old.a[level == j] != a[level == j])
    att.a[j] <- att.a[j] + sum(level == j)
    if ((i < burn / 2) & (att.a[j] > 200)) {
      if (acc.a[j] / att.a[j] < 0.3) { mh.a[j] <- mh.a[j] * 0.9 }
      if (acc.a[j] / att.a[j] > 0.6) { mh.a[j] <- mh.a[j] * 1.1 }
      acc.a[j] <- att.a[j] <- 0
    }
  }

  a          <- a.update$a
  theta.star <- a.update$theta.star
  cur.lly    <- a.update$cur.lly
  cur.llps   <- a.update$cur.llps

  a.keep[i, , ] <- a

  if ((i %% 100) == 0) {
    par(mfrow=c(3, 4))
    cat("Iter", i, "\n")
    if (i > 3000) {
      start <- i - 3000
    } else {
      start <- 1
    }
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
    plot(a.keep[start:i, 7, 10], type="l",
         main=paste("a 30, 1 true:", round(data$a[7, 10], 2)))
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
npts <- 50
u.beta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
bin.width <- u.beta[-1] - u.beta[-(npts + 1)]
x.beta <- matrix(2, ns, nt)
z.t <- getZ(xi=xi.t, x.beta=x.beta)
dw2 <- as.matrix(rdist(s, knots))^2
w.t <- stdW(makeW(dw2=dw2, rho=rho.t))
alpha <- 0.5
a <- matrix(1, nknots, nt)
theta.star <- getThetaStar(w=w.t, a=a, alpha=alpha)

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
                      cur.lly=cur.lly, cur.llps=cur.llps, z=z.t, w=w.t,
                      mid.points=mid.points, bin.width=bin.width,
                      mh=mh.a, cuts=cuts)

  a          <- a.update$a
  theta.star <- a.update$theta.star
  cur.lly    <- a.update$cur.lly
  cur.llps   <- a.update$cur.llps

  a.keep[i, , ] <- a

  alpha.update <- updateAlpha(y=data$y, theta.star=theta.star, a=a,
                              alpha=alpha, cur.lly=cur.lly, cur.llps=cur.llps,
                              z=z.t, w=w.t,
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
theta.star <- getThetaStar(w=w, a=data$a, alpha=alpha.t)
cur.lly <- logLikeY(y=data$y, theta.star=theta.star, alpha=alpha.t, z=z.t)

# MH adjustments
acc.rho <- att.rho <- mh.rho <- 0.1

# storage
rho.keep <- rep(NA, nreps)

for (i in 1:nreps) {
  rho.update <- updateRho(y=data$y, theta.star=theta.star, a=data$a,
                          alpha=alpha.t, cur.lly=cur.lly, z=z.t, w=w, dw2=dw2,
                          rho=rho, rho.upper=15,
                          acc=acc.rho, att=att.rho, mh=mh.rho)
  rho        <- rho.update$rho
  w          <- rho.update$w
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