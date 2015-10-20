rm(list=ls())
source("./hmc_aux.R")
source("./auxfunctions.R")
options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 10
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
rho.t <- 0.25
alpha.t <- 0.25
x <- matrix(1, ns, nt)
s <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0, 1, length = nknotsx), seq(0, 1, length = nknotsy))
dw2 <- rdist(s, knots)

others <- list(A.cutoff = max(sqrt(dw2)), thresh = 0)
data <- list(x = x, s = s, knots = knots, dw2 = dw2)
params.t <- list(rho = rho.t, alpha = alpha.t)
calc.t <- list()
calc.t$w <- getW(d = data, p = params.t, c = calc.t, o = others)
calc.t$z <- matrix(rgev(n = ns * nt, 1, 1, 1), ns, nt)
params.t$a <- matrix(rPS(nknots * nt, alpha = alpha.t), nknots, nt)
calc.t$aw <- getAW(d = data, p = params.t, c = calc.t, o = others)
calc.t$theta <- getTheta(d = data, p = params.t, c = calc.t, o = others)

y <- matrix(rbinom(ns * nt, size = 1, prob = -expm1(-calc.t$theta)), ns, nt) 

# create lists for MCMC
data   <- list(y = y, x = x, s = s, knots = knots, dw2 = dw2)
params <- list(beta = 0, xi = 0, rho = rho.t, alpha = alpha.t)  # still need a and b
calc   <- list(z = calc.t$z, x.beta = 0, w = calc.t$w)  # need aw, theta
priors <- list(beta.mn = 0, beta.sd = 100)

# initial values
a <- matrix(log(10), nknots, nt)
b <- matrix(0, nknots, nt)

params$a <- a
params$b <- b

calc$aw <- getAW(d = data, p = params, c = calc)
calc$theta <- getTheta(d = data, p = params, c = calc)

niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
set.seed(200)
tic.1 <- proc.time()
for (i in 1:niters) {
  q <- log(params$a)
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q, epsilon = 0.01, L = 10, 
                 dothers)
  if (HMCout$accept) {
    a   <- exp(q.a)
  }
  
  others <- list(y = y)
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
      plot(log(storage.a[1:i, idx, 1]), type = "l", main = round(log(a.t[idx, 1]), 2))
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


# try with the update for a and alpha
rm(list=ls())
library(Rcpp)
source("./auxfunctions.R")
source("./hmc_aux.R")
source("./updateModel.R")

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
theta.t <- getThetaCPP(wz = wz.t, a_star = a.t^alpha.t, alpha = alpha.t)
y   <- matrix(rbinom(ns * nt, size = 1, prob = -expm1(-theta.t)), ns, nt) 

niters <- 10000
burn   <- 100
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
a <- matrix(1, nknots, nt)
q.a <- log(a)
b <- matrix(0.5, nknots, nt)
q.b <- transform$logit(b)
q.ps <- c(as.vector(q.a), 0)
q.beta <- 0

x <- matrix(rep(1, ns * nt), nrow = ns, ncol = nt)
beta <- 0
xi <- 0
thresh <- 0
alpha <- 0.5
x.beta <- getXBeta(x = x, ns = ns, nt = nt, beta = beta)
z <- getZ(xi = xi, x.beta = x.beta, thresh = thresh)
wz <- getwzCPP(z = z, w = w.t)
att.beta <- acc.beta <- mh.beta <- 0.05
beta.m <- 0
beta.s <- 10
beta.attempts <- 50
cur.lly <- logLikeY(y = y, theta = theta.t)
storage.beta <- rep(NA, niters)


theta  <- theta.t
set.seed(200)
tic.1 <- proc.time()
Rprof(filename = "Rprof.out", line.profiling = TRUE)
for (iter in 1:niters) {
  aw <- getawCPP(a_star = a^alpha, w = w.t, alpha = alpha)
  others <- list(y = y, x = x, xi = xi, aw = aw, alpha = alpha, 
                 pri.mn = 0, pri.sd = 10)
  HMCout <- HMC(neg_log_post_beta, neg_log_post_grad_beta, q.beta, 
                epsilon = 0.005, L = 10, others)
  if (HMCout$accept) {
    q.beta <- HMCout$q
    beta   <- q.beta
    x.beta <- getXBeta(x = x, ns = ns, nt = nt, beta = beta)
    z <- getZ(xi = xi, x.beta = x.beta, thresh = 0)
  }
  
  wz <- getwzCPP(z = z, w = w.t)
  others <- list(y = y, alpha = alpha, b = b, a = a, wz = wz)
  HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q.ps, 
                 epsilon = 0.001, L = 10, others)
  if (HMCout$accept) {
    q.ps     <- HMCout$q
    q.a      <- matrix(q.ps[1:nkt], nknots, nt)
    q.alpha  <- tail(q.ps, 1)
    a <- exp(q.a)
    alpha <- transform$inv.logit(q.alpha)
  }
  
  others <- list(a = a, alpha = alpha)
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q.b, epsilon = 0.005, 
                 L = 10, others)
  if (HMCout$accept) {
    q.b <- HMCout$q
    b   <- transform$inv.logit(q.b)
  }
  storage.a[iter, , ] <- a
  storage.b[iter, , ] <- b
  storage.alpha[iter] <- alpha
  storage.beta[iter]  <- beta
  if (iter %% 500 == 0) {
    par(mfrow=c(3, 4))
    plot.idx <- seq(1, 5)
    for (idx in plot.idx){
      plot(log(storage.a[1:iter, idx, 1]), type = "l", main = round(log(a.t[idx, 1]), 2))
    }
    plot.idx <- seq(1, 5)
    for (idx in plot.idx){
      plot(storage.b[1:iter, idx, 1], type = "l")
    }
    plot(storage.beta[1:iter], type = "l")
    plot(storage.alpha[1:iter], type = "l")
    print(paste("iter:", iter, "of", niters, sep=" "))
  }
}
Rprof(filename = NULL)
summaryRprof(filename = "Rprof.out", lines = "show")
toc.1 <- proc.time()


# check to make sure they come back the same
#   y:      data
#   x:      covariates
#   xi:     xi
#   w:      w
#   alpha:  spatial dependence
#   wz:     kernel weights
#   a:      positive stable random effects
#   b:      auxiliary random variable
#   pri.mn: prior mean
#   pri.sd: prior standard deviation
xi.t <- 0
others1 <- list(y = y, x = x, xi = xi.t, w = w.t, alpha = alpha.t, wz = wz.t, a = a.t, b = b, pri.mn = 0, pri.sd = 10)
neg_log_post_beta(0.5, others = others1)

#   y(ns, nt):      data
#   x(ns, nt * np): covariates
#   xi(1):          xi
#   aw(ns, nt):     sum_l (A_l * w_li^(1 / alpha))
#   alpha(1):       spatial dependence
#   pri.mn(1):      prior mean
#   pri.sd(1):      prior standard deviation
aw <- getawCPP(a_star = a.t^alpha.t, w = w.t, alpha = alpha.t)
others2 <- list(y = y, x = x, xi = xi.t, aw = aw, alpha = alpha.t, pri.mn = 0, pri.sd = 10)
neg_log_post_beta_2(0.5, others = others2)

library(microbenchmark)
microbenchmark(neg_log_post_beta(0, others = others1), neg_log_post_beta_2(0, others = others2))
