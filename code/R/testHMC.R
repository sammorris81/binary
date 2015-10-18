rm(list=ls())
source("./hmc_aux.R")

options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1500
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
rho.t <- 0.25
alpha.t <- 0.25
s <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0, 1, length = nknotsx), seq(0, 1, length = nknotsy))
dw2 <- rdist(s, knots)
w.t <- stdW(makeW(dw2 = dw2, rho = rho.t))
z.t <- matrix(rgev(n = ns * nt, 1, 1, 1), ns, nt)
a.t <- matrix(rPS(nknots * nt, alpha = alpha.t), nknots, nt)
wz.t <- getwzCPP(z = z.t, w = w.t)
theta.t <- getThetaCPP(wz= wz.t, a_star = a.t^alpha.t, alpha = alpha.t)
y.t <- matrix(rbinom(ns * nt, size = 1, prob = -expm1(-theta.t)), ns, nt) 

niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
q.a <- matrix(log(10), nknots, nt)
q.b <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = q.b, a = q.a)

set.seed(200)
tic.1 <- proc.time()
for (i in 1:niters) {
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q.a, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    q.a      <- HMCout$q
    others$a <- exp(HMCout$q)
  }
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q.b, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    others$b <- transform$inv.logit(HMCout$q)
    q.b <- HMCout$q
  }
  storage.a[i, , ] <- others$a
  storage.b[i, , ] <- others$b
  if (i %% 100 == 0) {
    par(mfrow=c(3, 4))
    plot.idx <- seq(1, 23, by = 2)
    for (idx in plot.idx){
      plot(log(storage.a[1:i, idx, 1]), type = "l")
    }
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()


# try with the update for a and alpha
rm(list=ls())
source("./hmc_aux.R")

options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1500
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
nkt <- nknots * nt
rho.t <- 0.25
alpha.t <- 0.25
s <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0, 1, length = nknotsx), seq(0, 1, length = nknotsy))
dw2 <- rdist(s, knots)
w.t <- stdW(makeW(dw2 = dw2, rho = rho.t))
z.t <- matrix(rgev(n = ns * nt, 1, 1, 1), ns, nt)
a.t <- matrix(rPS(nknots * nt, alpha = alpha.t), nknots, nt)
wz.t <- getwzCPP(z = z.t, w = w.t)
theta.t <- getThetaCPP(wz= wz.t, a_star = a.t^alpha.t, alpha = alpha.t)
y.t <- matrix(rbinom(ns * nt, size = 1, prob = -expm1(-theta.t)), ns, nt) 

niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
q.a <- matrix(1, nknots, nt)
q.b <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = q.b, a = q.a)

q.ps <- c(as.vector(q.a), 0)

set.seed(200)
tic.1 <- proc.time()
for (i in 1:niters) {
  HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q.ps, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    q.ps     <- HMCout$q
    q.a      <- matrix(q.ps[1:nkt], nknots, nt)
    q.alpha  <- tail(q.ps, 1)
    others$a <- exp(q.a)
    others$alpha <- transform$inv.logit(q.alpha)
  }
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q.b, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    others$b <- transform$inv.logit(HMCout$q)
    q.b <- HMCout$q
  }
  storage.a[i, , ] <- others$a
  storage.b[i, , ] <- others$b
  storage.alpha[i] <- others$alpha
  if (i %% 500 == 0) {
    par(mfrow=c(3, 4))
    plot.idx <- seq(1, 6)
    for (idx in plot.idx){
      plot(log(storage.a[1:i, idx, 1]), type = "l", main = round(log(a.t[idx, 1]), 2))
    }
    plot.idx <- seq(1, 5)
    for (idx in plot.idx){
      plot(storage.b[1:i, idx, 1], type = "l")
    }
    plot(storage.alpha[1:i], type = "l")
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()


# try with the update for a and alpha
rm(list=ls())
source("./hmc_aux.R")

options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1500
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
nkt <- nknots * nt
rho.t <- 0.25
alpha.t <- 0.8
s <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0, 1, length = nknotsx), seq(0, 1, length = nknotsy))
dw2 <- rdist(s, knots)
w.t <- stdW(makeW(dw2 = dw2, rho = rho.t))
z.t <- matrix(rgev(n = ns * nt, 1, 1, 1), ns, nt)
a.t <- matrix(rPS(nknots * nt, alpha = alpha.t), nknots, nt)
wz.t <- getwzCPP(z = z.t, w = w.t)
theta.t <- getThetaCPP(wz= wz.t, a_star = a.t^alpha.t, alpha = alpha.t)
y.t <- matrix(rbinom(ns * nt, size = 1, prob = -expm1(-theta.t)), ns, nt) 

niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
q.a <- matrix(1, nknots, nt)
q.b <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = q.b, a = q.a)

q.ps <- c(as.vector(q.a), 0)

set.seed(200)
tic.1 <- proc.time()
for (i in 1:niters) {
  HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q.ps, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    q.ps     <- HMCout$q
    q.a      <- matrix(q.ps[1:nkt], nknots, nt)
    q.alpha  <- tail(q.ps, 1)
    others$a <- exp(q.a)
    others$alpha <- transform$inv.logit(q.alpha)
  }
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q.b, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    others$b <- transform$inv.logit(HMCout$q)
    q.b <- HMCout$q
  }
  storage.a[i, , ] <- others$a
  storage.b[i, , ] <- others$b
  storage.alpha[i] <- others$alpha
  if (i %% 500 == 0) {
    par(mfrow=c(3, 4))
    plot.idx <- seq(1, 6)
    for (idx in plot.idx){
      plot(log(storage.a[1:i, idx, 1]), type = "l", main = round(log(a.t[idx, 1]), 2))
    }
    plot.idx <- seq(1, 5)
    for (idx in plot.idx){
      plot(storage.b[1:i, idx, 1], type = "l")
    }
    plot(storage.alpha[1:i], type = "l")
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()
