rm(list=ls())
source("./hmc_aux.R")
source("./auxfunctions.R")
options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1000
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
prior  <- list(beta.mn = 0, beta.sd = 100)

# initial values
a <- matrix(1, nknots, nt)
b <- matrix(0.5, nknots, nt)

params$a <- a
params$b <- b

calc$aw <- getAW(d = data, p = params, c = calc)
calc$theta <- getTheta(d = data, p = params, c = calc)


library(numDeriv)
neg_log_post_a(q = log(params$a), d = data, p = params, c = calc, o = others, 
               prior = prior)
neg_log_post_grad_a(q = log(params$a), d = data, p = params, c = calc, 
                    o = others, prior = prior)
grad(func = neg_log_post_a, x = log(params$a), d = data, p = params, c = calc, 
     o = others, prior = prior)

neg_log_post_b(q = transform$logit(params$b), d = data, p = params, c = calc,
               o = others, prior = prior)
neg_log_post_grad_b(q = transform$logit(params$b), d = data, p = params, 
                    c = calc, o = others, prior = prior)
grad(func = neg_log_post_b, x = transform$logit(params$b), d = data, p = params, 
     c = calc, o = others, prior = prior)

niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
set.seed(200)

tic.1 <- proc.time()
for (i in 1:niters) {
  q <- log(params$a)
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q, epsilon = 0.005, L = 10, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior)
  if (HMCout$accept) {
    params$a <- exp(HMCout$q)
    calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
  q <- transform$logit(params$b)
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q, epsilon = 0.01, L = 10, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior)
  if (HMCout$accept) {
    params$b <- transform$inv.logit(HMCout$q)
  }
  storage.a[i, , ] <- params$a
  storage.b[i, , ] <- params$b
  if (i %% 500 == 0) {
    par(mfrow=c(3, 4))
    plot.idx <- seq(1, 23, by = 2)
    for (idx in plot.idx){
      plot(log(storage.a[1:i, idx, 1]), type = "l", 
           main = round(log(params.t$a[idx, 1]), 2))
    }
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()


# try with the update for a and alpha
rm(list=ls())
source("./hmc_aux.R")
source("./auxfunctions.R")
options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1000
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
nkt <- nknots * nt
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
prior  <- list(beta.mn = 0, beta.sd = 2)

# initial values
a <- matrix(1, nknots, nt)
b <- matrix(0.5, nknots, nt)
alpha <- 0.5

params$a <- a
params$b <- b
params$alpha <- 0.5

calc$aw <- getAW(d = data, p = params, c = calc)
calc$theta <- getTheta(d = data, p = params, c = calc)

library(numDeriv)
neg_log_post_a_alpha(q = as.vector(c(log(params$a), transform$logit(params$alpha))), 
                     d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a_alpha(q = as.vector(c(log(params$a), transform$logit(params$alpha))), 
                          d = data, p = params, c = calc, o = others, prior = prior)
grad(func = neg_log_post_a_alpha, x = as.vector(c(log(params$a), transform$logit(params$alpha))), 
     d = data, p = params, c = calc, o = others, prior = prior)
library(microbenchmark)
microbenchmark(neg_log_post_grad_a_alpha(q = as.vector(c(log(params$a), transform$logit(params$alpha))), 
                                         d = data, p = params, c = calc, o = others, prior = prior))

neg_log_post_b(q = transform$logit(params$b), d = data, p = params, c = calc,
               o = others, prior = prior)
neg_log_post_grad_b(q = transform$logit(params$b), d = data, p = params, 
                    c = calc, o = others, prior = prior)
grad(func = neg_log_post_b, x = transform$logit(params$b), d = data, p = params, 
     c = calc, o = others, prior = prior)

niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
set.seed(200)

for (i in 1:niters) {
  q <- as.vector(c(log(params$a), transform$logit(params$alpha)))
  HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q, epsilon = 0.005, L = 10, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior)
  if (HMCout$accept) {
    params$a <- matrix(exp(HMCout$q[1:nkt]), nknots, nt)
    params$alpha <- transform$inv.logit(tail(HMCout$q, 1))
    calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
  q <- transform$logit(params$b)
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q, epsilon = 0.01, L = 10, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior)
  if (HMCout$accept) {
    params$b <- transform$inv.logit(HMCout$q)
  }
  storage.a[i, , ] <- params$a
  storage.b[i, , ] <- params$b
  storage.alpha[i] <- params$alpha
  
  if (i %% 500 == 0) {
    par(mfrow=c(3, 4))
    plot.idx <- seq(1, 6)
    for (idx in plot.idx){
      plot(log(storage.a[1:i, idx, 1]), type = "l", 
           main = round(log(params.t$a[idx, 1]), 2))
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


# try with the update for beta, a, alpha
rm(list=ls())
source("./hmc_aux.R")
source("./auxfunctions.R")
options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1000
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
nkt <- nknots * nt
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

gen <- rRareBinarySpat(x = x, s = s, knots = knots, beta = 0, xi = 0, alpha = alpha.t,
                     rho = rho.t, nt = 1, prob.success = 0.05, dw2 = dw2)

# create lists for MCMC
data   <- list(y = gen$y, x = x, s = s, knots = knots, dw2 = dw2)
params <- list(beta = 0, xi = 0, rho = rho.t, alpha = alpha.t)  # still need a and b
calc   <- list(z = calc.t$z, x.beta = 0, w = calc.t$w)  # need aw, theta
prior  <- list(beta.mn = 0, beta.sd = 10)

# initial values
a <- matrix(1, nknots, nt)
b <- matrix(0.5, nknots, nt)
alpha <- 0.5
beta  <- 0

params$a <- a
params$b <- b
params$alpha <- alpha
params$beta <- beta

calc$x.beta <- getXBeta(d = data, p = params, c = calc, o = others)
calc$z      <- getZ(d = data, p = params, c = calc, o = others)
calc$aw     <- getAW(d = data, p = params, c = calc, o = others)
calc$theta  <- getTheta(d = data, p = params, c = calc, o = others)

library(numDeriv)
neg_log_post_a(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a_alpha(q = as.vector(c(log(params$a), transform$logit(params$alpha))), 
                          d = data, p = params, c = calc, o = others, prior = prior)
grad(func = neg_log_post_a_alpha, x = as.vector(c(log(params$a), transform$logit(params$alpha))), 
     d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_alpha(q = transform$logit(params$alpha), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_alpha(q = transform$logit(params$alpha), d = data, p = params, c = calc, o = others, prior = prior)
grad(func = neg_log_post_alpha, x = transform$logit(params$alpha), d = data, p = params, c = calc, o = others, prior = prior)

neg_log_post_beta(q = params$beta, d = data, p = params, c = calc, o = others, 
                  prior = prior)
neg_log_post_grad_beta(q = params$beta, d = data, p = params, c = calc, o = others, 
                       prior = prior)
grad(func = neg_log_post_beta, x = params$beta, d = data, p = params, c = calc, o = others, 
     prior = prior)


niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
storage.beta  <- rep(NA, niters)
set.seed(200)

for (i in 1:niters) {
  q <- params$beta
  HMCout <- HMC(neg_log_post_beta, neg_log_post_grad_beta, q, epsilon = 0.01, 
                L = 10, data = data, params = params, calc = calc, others = others, 
                prior = prior)
  if (HMCout$accept) {
    params$beta <- HMCout$q
    calc$x.beta <- getXBeta(d = data, p = params, c = calc, o = others)
    calc$z <- getZ(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
  q <- log(params$a)
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q, epsilon = 0.005, L = 10, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior)
  if (HMCout$accept) {
    params$a <- exp(HMCout$q)
    calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
  q <- transform$logit(params$alpha)
  HMCout  <- HMC(neg_log_post_alpha, neg_log_post_grad_alpha, q, epsilon = 0.0005, L = 10, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior)
  if (HMCout$accept) {
    alpha <- transform$inv.logit(HMCout$q)
    if (alpha > 0.1 & alpha < 0.9) {
      params$alpha <- transform$inv.logit(HMCout$q)
      calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
      calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
    }
  }
  
  q <- transform$logit(params$b)
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q, epsilon = 0.01, L = 10, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior)
  if (HMCout$accept) {
    params$b <- transform$inv.logit(HMCout$q)
  }
  
  storage.a[i, , ] <- params$a
  storage.b[i, , ] <- params$b
  storage.alpha[i] <- params$alpha
  storage.beta[i]  <- params$beta
  
  if (i %% 500 == 0) {
    par(mfrow=c(3, 4))
    plot.idx <- seq(1, 5)
    for (idx in plot.idx){
      plot(log(storage.a[1:i, idx, 1]), type = "l", 
           main = round(log(gen$a[idx, 1]), 2))
    }
    plot.idx <- seq(1, 5)
    for (idx in plot.idx){
      plot(storage.b[1:i, idx, 1], type = "l")
    }
    plot(storage.beta[1:i], type = "l")
    plot(storage.alpha[1:i], type = "l")
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()

# try with the update for beta, a, alpha
rm(list=ls())
source("./hmc_aux.R")
source("./auxfunctions.R")
options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1000
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
nkt <- nknots * nt
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

gen <- rRareBinarySpat(x = x, s = s, knots = knots, beta = 0, xi = 0, alpha = alpha.t,
                       rho = rho.t, nt = 1, prob.success = 0.05, dw2 = dw2)

# create lists for MCMC
data   <- list(y = gen$y, x = x, s = s, knots = knots, dw2 = dw2)
params <- list(beta = 0, xi = 0, rho = rho.t, alpha = alpha.t)  # still need a and b
calc   <- list(z = calc.t$z, x.beta = 0, w = calc.t$w)  # need aw, theta
prior  <- list(beta.mn = 0, beta.sd = 10)

# initial values
a <- matrix(100, nknots, nt)
b <- matrix(0.5, nknots, nt)
alpha <- 0.5
beta  <- 0

params$a <- a
params$b <- b
params$alpha <- alpha
params$beta <- beta

calc$x.beta <- getXBeta(d = data, p = params, c = calc, o = others)
calc$z      <- getZ(d = data, p = params, c = calc, o = others)
calc$aw     <- getAW(d = data, p = params, c = calc, o = others)
calc$theta  <- getTheta(d = data, p = params, c = calc, o = others)

library(numDeriv)
neg_log_post_a(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a_alpha(q = as.vector(c(log(params$a), transform$logit(params$alpha))), 
                          d = data, p = params, c = calc, o = others, prior = prior)
grad(func = neg_log_post_a_alpha, x = as.vector(c(log(params$a), transform$logit(params$alpha))), 
     d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_alpha(q = transform$logit(params$alpha), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_alpha(q = transform$logit(params$alpha), d = data, p = params, c = calc, o = others, prior = prior)
grad(func = neg_log_post_alpha, x = transform$logit(params$alpha), d = data, p = params, c = calc, o = others, prior = prior)

neg_log_post_beta(q = params$beta, d = data, p = params, c = calc, o = others, 
                  prior = prior)
neg_log_post_grad_beta(q = params$beta, d = data, p = params, c = calc, o = others, 
                       prior = prior)
grad(func = neg_log_post_beta, x = params$beta, d = data, p = params, c = calc, o = others, 
     prior = prior)


niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
storage.beta  <- rep(NA, niters)
storage.prob  <- array(NA, dim = c(niters, ns, nt))
set.seed(200)

beta.att  <- beta.acc  <- beta.eps  <- 0.01
a.att     <- a.acc     <- a.eps     <- 0.05       
alpha.att <- alpha.acc <- alpha.eps <- 0.01
b.att     <- b.acc     <- b.eps     <- 0.05

for (i in 1:niters) {
  beta.att <- beta.att + 1
  q <- params$beta
  HMCout <- HMC(neg_log_post_beta, neg_log_post_grad_beta, q, epsilon = beta.eps, 
                L = 20, data = data, params = params, calc = calc, 
                others = others, prior = prior, this.param = "beta")
  if (HMCout$accept) {
    beta.acc <- beta.acc + 1
    params$beta <- HMCout$q
    calc$x.beta <- getXBeta(d = data, p = params, c = calc, o = others)
    calc$z <- getZ(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
  #   if (beta.att > 500) {
  #     beta.rate <- beta.acc / beta.att
  #     if (beta.rate < 0.20) {
  #       beta.eps <- beta.eps * 0.8
  #     } else if (beta.rate > 0.60) {
  #       beta.eps <- beta.eps * 1.2
  #     }
  #     beta.acc <- beta.att <- 0
  #   }
  
  a.att <- a.att + 1
  q <- log(params$a)
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q, 
                 epsilon = a.eps, L = 30, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior, this.param = "a")
  if (HMCout$accept) {
    a.acc <- a.acc + 1
    params$a <- exp(HMCout$q)
    calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
  #   if (a.att > 100) {
  #     a.rate <- a.acc / a.att
  #     if (a.rate < 0.20) {
  #       a.eps <- a.eps * 0.8
  #     } else if (a.rate > 0.60) {
  #       a.eps <- a.eps * 1.2
  #     }
  #     a.acc <- a.att <- 0
  #   }
  
  alpha.att <- alpha.att + 1
  q <- transform$logit(params$alpha)
  HMCout  <- HMC(neg_log_post_alpha, neg_log_post_grad_alpha, q, 
                 epsilon = alpha.eps, L = 30, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior, this.param = "alpha")
  if (HMCout$accept) {
    alpha.acc <- alpha.acc + 1
    params$alpha <- transform$inv.logit(HMCout$q)
    calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
  #   if (alpha.att > 500) {
  #     alpha.rate <- alpha.acc / alpha.att
  #     if (alpha.rate < 0.20) {
  #       alpha.eps <- alpha.eps * 0.8
  #     } else if (alpha.rate > 0.60) {
  #       alpha.eps <- alpha.eps * 1.2
  #     }
  #     alpha.acc <- alpha.att <- 0
  #   }
  
  #   a_alpha.att <- a_alpha.att + 1
  #   q <- as.vector(c(log(params$a), transform$logit(params$alpha)))
  #   HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q, 
  #                  epsilon = a_alpha.eps, L = 30, 
  #                  data = data, params = params, calc = calc, others = others, 
  #                  prior = prior, this.param = "a and alpha")
  #   if (HMCout$accept) {
  #     a_alpha.acc <- a_alpha.acc + 1
  #     params$a <- matrix(exp(HMCout$q[1:nkt]), nknots, nt)
  #     params$alpha <- transform$inv.logit(tail(HMCout$q, 1))
  #     calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
  #     calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  #   }
  
  #   if (a_alpha.att > 100) {
  #     a_alpha.rate <- a_alpha.acc / a_alpha.att
  #     if (a_alpha.rate < 0.20) {
  #       a_alpha.eps <- a_alpha.eps * 0.8
  #     } else if (a_alpha.rate > 0.60) {
  #       a_alpha.eps <- a_alpha.eps * 1.2
  #     }
  #     a_alpha.acc <- a_alpha.att <- 0
  #   }
  
  q <- transform$logit(params$b)
  b.att <- b.att + 1
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q, epsilon = b.eps, 
                 L = 10, data = data, params = params, calc = calc, 
                 others = others, prior = prior, this.param = "b")
  if (HMCout$accept) {
    b.acc <- b.acc + 1
    params$b <- transform$inv.logit(HMCout$q)
  }
  
  #   if (b.att > 100) {
  #     b.rate <- b.acc / b.att
  #     if (b.rate < 0.20) {
  #       b.eps <- b.eps * 0.8
  #     } else if (b.rate > 0.60) {
  #       b.eps <- b.eps * 1.2
  #     }
  #     b.acc <- b.att <- 0
  #   }
  
  storage.a[i, , ] <- params$a
  storage.b[i, , ] <- params$b
  storage.alpha[i] <- params$alpha
  storage.beta[i]  <- params$beta
  storage.prob[i, , ] <- 1 - exp(-calc$theta)
  
  if (i %% 500 == 0) {
    start <- max(i - 5000, 1)
    end   <- i
    par(mfrow=c(4, 5))
    plot.idx <- seq(1, 18, by = 2)
    for (idx in plot.idx){
      plot(log(storage.a[1:i, idx, 1]), type = "l", 
           main = round(log(gen$a[idx, 1]), 2), 
           xlab = round(a.acc / a.att, 3))
    }
    plot.idx <- seq(1, 18, by = 2)
    for (idx in plot.idx){
      plot(storage.b[1:i, idx, 1], type = "l", 
           xlab = round(b.acc / b.att, 3))
    }
    #     plot.idx <- 1:18
    #     for (idx in plot.idx){
    #       plot(storage.prob[start:end, idx, 1], type = "l")
    #     }
    plot(storage.beta[start:end], type = "l", 
         xlab = round(beta.acc / beta.att, 2))
    plot(storage.alpha[start:end], type = "l",
         xlab = round(alpha.acc / alpha.att, 2))
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()

par(mfrow = c(5, 5))
plot.idx <- seq(1, 46, by = 2)
for (idx in plot.idx) {
  plot(storage.prob[5001:30000, idx, 1], type = "l",
       main = paste("Site ", idx), ylab = "P(Y = 1)")
}
plot(storage.alpha[5001:30000], type = "l", main = bquote(alpha))
plot(storage.beta[5001:30000], type = "l", main = bquote(beta))



# try with the update for beta, a, alpha (higher rate of occurrence)
rm(list=ls())
source("./hmc_aux.R")
source("./auxfunctions.R")
options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1000
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
nkt <- nknots * nt
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

gen <- rRareBinarySpat(x = x, s = s, knots = knots, beta = 0, xi = 0, alpha = alpha.t,
                       rho = rho.t, nt = 1, prob.success = 0.5, dw2 = dw2)

# create lists for MCMC
data   <- list(y = gen$y, x = x, s = s, knots = knots, dw2 = dw2)
params <- list(beta = 0, xi = 0, rho = rho.t, alpha = alpha.t)  # still need a and b
calc   <- list(z = calc.t$z, x.beta = 0, w = calc.t$w)  # need aw, theta
prior  <- list(beta.mn = 0, beta.sd = 10)

# initial values
a <- matrix(10, nknots, nt)
b <- matrix(0.5, nknots, nt)
alpha <- 0.5
beta  <- 0

params$a <- a
params$b <- b
params$alpha <- alpha
params$beta <- beta

calc$x.beta <- getXBeta(d = data, p = params, c = calc, o = others)
calc$z      <- getZ(d = data, p = params, c = calc, o = others)
calc$aw     <- getAW(d = data, p = params, c = calc, o = others)
calc$theta  <- getTheta(d = data, p = params, c = calc, o = others)

library(numDeriv)
neg_log_post_a(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a_alpha(q = as.vector(c(log(params$a), transform$logit(params$alpha))), 
                          d = data, p = params, c = calc, o = others, prior = prior)
grad(func = neg_log_post_a_alpha, x = as.vector(c(log(params$a), transform$logit(params$alpha))), 
     d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_alpha(q = transform$logit(params$alpha), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_alpha(q = transform$logit(params$alpha), d = data, p = params, c = calc, o = others, prior = prior)
grad(func = neg_log_post_alpha, x = transform$logit(params$alpha), d = data, p = params, c = calc, o = others, prior = prior)

neg_log_post_beta(q = params$beta, d = data, p = params, c = calc, o = others, 
                  prior = prior)
neg_log_post_grad_beta(q = params$beta, d = data, p = params, c = calc, o = others, 
                       prior = prior)
grad(func = neg_log_post_beta, x = params$beta, d = data, p = params, c = calc, o = others, 
     prior = prior)


niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
storage.beta  <- rep(NA, niters)
storage.prob  <- array(NA, dim = c(niters, ns, nt))
set.seed(200)

beta.att  <- beta.acc  <- beta.eps  <- 0.01
a.att     <- a.acc     <- a.eps     <- 0.1       
alpha.att <- alpha.acc <- alpha.eps <- 0.005
b.att     <- b.acc     <- b.eps     <- 0.3

for (i in 1:niters) {
  beta.att <- beta.att + 1
  q <- params$beta
  HMCout <- HMC(neg_log_post_beta, neg_log_post_grad_beta, q, epsilon = beta.eps, 
                L = 20, data = data, params = params, calc = calc, 
                others = others, prior = prior, this.param = "beta")
  if (HMCout$accept) {
    beta.acc <- beta.acc + 1
    params$beta <- HMCout$q
    calc$x.beta <- getXBeta(d = data, p = params, c = calc, o = others)
    calc$z <- getZ(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
#   if (beta.att > 500) {
#     beta.rate <- beta.acc / beta.att
#     if (beta.rate < 0.20) {
#       beta.eps <- beta.eps * 0.8
#     } else if (beta.rate > 0.60) {
#       beta.eps <- beta.eps * 1.2
#     }
#     beta.acc <- beta.att <- 0
#   }
  
  a.att <- a.att + 1
  q <- log(params$a)
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q, 
                 epsilon = a.eps, L = 30, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior, this.param = "a")
  if (HMCout$accept) {
    a.acc <- a.acc + 1
    params$a <- exp(HMCout$q)
    calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
#   if (a.att > 100) {
#     a.rate <- a.acc / a.att
#     if (a.rate < 0.20) {
#       a.eps <- a.eps * 0.8
#     } else if (a.rate > 0.60) {
#       a.eps <- a.eps * 1.2
#     }
#     a.acc <- a.att <- 0
#   }
  
  alpha.att <- alpha.att + 1
  q <- transform$logit(params$alpha)
  HMCout  <- HMC(neg_log_post_alpha, neg_log_post_grad_alpha, q, 
                 epsilon = alpha.eps, L = 30, 
                 data = data, params = params, calc = calc, others = others, 
                 prior = prior, this.param = "alpha")
  if (HMCout$accept) {
    alpha.acc <- alpha.acc + 1
    params$alpha <- transform$inv.logit(HMCout$q)
    calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
    calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  }
  
#   if (alpha.att > 500) {
#     alpha.rate <- alpha.acc / alpha.att
#     if (alpha.rate < 0.20) {
#       alpha.eps <- alpha.eps * 0.8
#     } else if (alpha.rate > 0.60) {
#       alpha.eps <- alpha.eps * 1.2
#     }
#     alpha.acc <- alpha.att <- 0
#   }
  
#   a_alpha.att <- a_alpha.att + 1
#   q <- as.vector(c(log(params$a), transform$logit(params$alpha)))
#   HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q, 
#                  epsilon = a_alpha.eps, L = 30, 
#                  data = data, params = params, calc = calc, others = others, 
#                  prior = prior, this.param = "a and alpha")
#   if (HMCout$accept) {
#     a_alpha.acc <- a_alpha.acc + 1
#     params$a <- matrix(exp(HMCout$q[1:nkt]), nknots, nt)
#     params$alpha <- transform$inv.logit(tail(HMCout$q, 1))
#     calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
#     calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
#   }
  
#   if (a_alpha.att > 100) {
#     a_alpha.rate <- a_alpha.acc / a_alpha.att
#     if (a_alpha.rate < 0.20) {
#       a_alpha.eps <- a_alpha.eps * 0.8
#     } else if (a_alpha.rate > 0.60) {
#       a_alpha.eps <- a_alpha.eps * 1.2
#     }
#     a_alpha.acc <- a_alpha.att <- 0
#   }
  
  q <- transform$logit(params$b)
  b.att <- b.att + 1
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q, epsilon = b.eps, 
                 L = 10, data = data, params = params, calc = calc, 
                 others = others, prior = prior, this.param = "b")
  if (HMCout$accept) {
    b.acc <- b.acc + 1
    params$b <- transform$inv.logit(HMCout$q)
  }
  
#   if (b.att > 100) {
#     b.rate <- b.acc / b.att
#     if (b.rate < 0.20) {
#       b.eps <- b.eps * 0.8
#     } else if (b.rate > 0.60) {
#       b.eps <- b.eps * 1.2
#     }
#     b.acc <- b.att <- 0
#   }
  
  storage.a[i, , ] <- params$a
  storage.b[i, , ] <- params$b
  storage.alpha[i] <- params$alpha
  storage.beta[i]  <- params$beta
  storage.prob[i, , ] <- 1 - exp(-calc$theta)
  
  if (i %% 500 == 0) {
    start <- max(i - 5000, 1)
    end   <- i
    par(mfrow=c(4, 5))
        plot.idx <- seq(1, 18, by = 2)
        for (idx in plot.idx){
          plot(log(storage.a[1:i, idx, 1]), type = "l", 
               main = round(log(gen$a[idx, 1]), 2), 
               xlab = round(a.acc / a.att, 3))
        }
        plot.idx <- seq(1, 18, by = 2)
        for (idx in plot.idx){
          plot(storage.b[1:i, idx, 1], type = "l", 
               xlab = round(b.acc / b.att, 3))
        }
#     plot.idx <- 1:18
#     for (idx in plot.idx){
#       plot(storage.prob[start:end, idx, 1], type = "l")
#     }
    plot(storage.beta[start:end], type = "l", 
         xlab = round(beta.acc / beta.att, 2))
    plot(storage.alpha[start:end], type = "l",
         xlab = round(alpha.acc / alpha.att, 2))
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()

par(mfrow = c(5, 5))
plot.idx <- seq(1, 46, by = 2)
for (idx in plot.idx) {
  plot(storage.prob[5001:30000, idx, 1], type = "l",
       main = paste("Site ", idx), ylab = "P(Y = 1)")
}
plot(storage.alpha[5001:30000], type = "l", main = bquote(alpha))
plot(storage.beta[5001:30000], type = "l", main = bquote(beta))


# try with the update for beta (MH), a, alpha (higher rate of occurrence)
rm(list=ls())
source("./hmc_aux.R")
source("./updateModel.R")
source("./auxfunctions.R")
options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1000
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
nkt <- nknots * nt
rho.t <- list(cur = 0.25)
alpha.t <- list(cur = 0.25)
x <- matrix(1, ns, nt)
s <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0, 1, length = nknotsx), seq(0, 1, length = nknotsy))
dw2 <- rdist(s, knots)

others <- list(A.cutoff = max(sqrt(dw2)), thresh = 0, dw2 = dw2)
data <- list(x = x, s = s, knots = knots)
calc.t <- list()
calc.t$w <- getW(rho = rho.t, others = others)
calc.t$z <- matrix(rgev(n = ns * nt, 1, 1, 1), ns, nt)
a <- list(cur = matrix(rPS(nknots * nt, alpha = alpha.t$cur), nknots, nt))
calc.t$aw <- getAW(alpha = alpha.t, a = a, calc = calc.t)
calc.t$theta <- getTheta(alpha = alpha.t, calc = calc.t)

gen <- rRareBinarySpat(x = x, s = s, knots = knots, beta = 0, xi = 0, alpha = alpha.t,
                       rho = rho.t, nt = 1, prob.success = 0.05, dw2 = dw2)

# create lists for MCMC
data   <- list(y = gen$y, x = x, s = s, knots = knots)
calc   <- list(w = calc.t$w)  # need aw, theta

# initial values
beta.init <- -log(-log(mean(data$y)))
beta  <- list(cur = beta.init, att = 0, acc = 0, eps = 0.5, mn = 0, sd = 100)
xi    <- list(cur = 0, att = 0, acc = 0, eps = 0.01, mn = 0, sd = 0.5)
a     <- list(cur = matrix(10, nknots, nt), att = 0, acc = 0, eps = 0.3)
b     <- list(cur = matrix(0.5, nknots, nt), att = 0, acc = 0, eps = 0.3)
alpha <- list(cur = 0.5, att = 0, acc = 0, eps = 0.005)
rho   <- list(cur = 0.1, att = 0, acc = 0, eps = 0.01, mn = 0, sd = 1)

calc$x.beta <- getXBeta(data = data, beta = beta)
calc$z      <- getZ(xi = xi, calc = calc, others = others)
calc$aw     <- getAW(alpha = alpha, a = a, calc = calc)
calc$theta  <- getTheta(alpha = alpha, calc = calc)

library(numDeriv)
neg_log_post_a(q = log(a$cur), data = data, beta = beta, xi = xi, 
               a = a, b = b, alpha = alpha, rho = rho, calc = calc, others = others)
neg_log_post_grad_a(q = log(a$cur), data = data, beta = beta, xi = xi, a = a, 
                    b = b, alpha = alpha, rho = rho, calc = calc, others = others)
neg_log_post_grad_a_alpha(q = as.vector(c(log(a$cur), transform$logit(alpha$cur))),
                          data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha,
                          rho = rho, calc = calc, others = others, eps = 0.0001)

grad(func = neg_log_post_a_alpha, 
     x = as.vector(c(log(a$cur), transform$logit(alpha$cur))), 
     data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, rho = rho, calc = calc,
     others = others)
neg_log_post_alpha(q = transform$logit(alpha$cur), data = data, beta = beta, xi = xi, 
                   a = a, b = b, alpha = alpha, rho = rho, calc = calc, others = others)
neg_log_post_grad_alpha(q = transform$logit(alpha$cur), data = data, beta = beta, xi = xi, 
                        a = a, b = b, alpha = alpha, rho = rho, calc = calc, others = others, 
                        eps = 0.0001)
grad(func = neg_log_post_alpha, x = transform$logit(alpha$cur), 
     data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, rho = rho, 
     calc = calc, others = others)

neg_log_post_beta(q = beta$cur, data = data, beta = beta, xi = xi, a = a, 
                  b = b, alpha = alpha, rho = rho, calc = calc, others = others)
neg_log_post_grad_beta(q = beta$cur, data = data, beta = beta, xi = xi, a = a,
                       b = b, alpha = alpha, rho = rho, calc = calc, 
                       others = others, eps = 0.0001)
grad(func = neg_log_post_beta, x = beta$cur, data = data, beta = beta, xi = xi,
     a = a, b = b, alpha = alpha, rho = rho, calc = calc, others = others)


niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
storage.beta  <- rep(NA, niters)
storage.prob  <- array(NA, dim = c(niters, ns, nt))

set.seed(200)
for (i in 1:niters) {
  beta$att <- beta$att + 1
  q <- beta$cur
  MHout <- updateBeta(data = data, beta = beta, xi = xi, alpha = alpha, 
                      calc = calc, others = others)
  if (MHout$accept) {
    beta$acc    <- beta$acc + 1
    beta$cur    <- MHout$q
    calc$x.beta <- getXBeta(data = data, beta = beta)
    calc$z      <- getZ(xi = xi, calc = calc, others = others)
    calc$theta  <- getTheta(alpha = alpha, calc = calc)
  }
  
  if (beta$att > 100) {
    beta.rate <- beta$acc / beta$att
    if (beta.rate < 0.20) {
      beta$eps <- beta$eps * 0.8
    } else if (beta.rate > 0.60) {
      beta$eps <- beta$eps * 1.2
    }
    beta$acc <- beta$att <- 0
  }
  
  a$att <- a$att + 1
  q <- log(a$cur)
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q, 
                 epsilon = a$eps, L = 30, 
                 data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                 rho = rho, calc = calc, others = others, this.param = "a")
  if (HMCout$accept) {
    a$acc <- a$acc + 1
    a$cur <- exp(HMCout$q)
    calc$aw  <- getAW(alpha = alpha, a = a, calc = calc)
    calc$theta <- getTheta(alpha = alpha, calc = calc)
  }
  
  #   if (a.att > 100) {
  #     a.rate <- a.acc / a.att
  #     if (a.rate < 0.20) {
  #       a.eps <- a.eps * 0.8
  #     } else if (a.rate > 0.60) {
  #       a.eps <- a.eps * 1.2
  #     }
  #     a.acc <- a.att <- 0
  #   }
  
  alpha$att <- alpha$att + 1
  q <- transform$logit(alpha$cur)
  HMCout  <- HMC(neg_log_post_alpha, neg_log_post_grad_alpha, q, 
                 epsilon = alpha$eps, L = 30, 
                 data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                 rho = rho, calc = calc, others = others, this.param = "alpha")
  if (HMCout$accept) {
    alpha$acc <- alpha$acc + 1
    alpha$cur <- transform$inv.logit(HMCout$q)
    calc$aw  <- getAW(alpha = alpha, a = a, calc = calc)
    calc$theta <- getTheta(alpha = alpha, calc = calc)
  }
  
  #   if (alpha.att > 500) {
  #     alpha.rate <- alpha.acc / alpha.att
  #     if (alpha.rate < 0.20) {
  #       alpha.eps <- alpha.eps * 0.8
  #     } else if (alpha.rate > 0.60) {
  #       alpha.eps <- alpha.eps * 1.2
  #     }
  #     alpha.acc <- alpha.att <- 0
  #   }
  
  #   a_alpha.att <- a_alpha.att + 1
  #   q <- as.vector(c(log(params$a), transform$logit(params$alpha)))
  #   HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q, 
  #                  epsilon = a_alpha.eps, L = 30, 
  #                  data = data, params = params, calc = calc, others = others, 
  #                  prior = prior, this.param = "a and alpha")
  #   if (HMCout$accept) {
  #     a_alpha.acc <- a_alpha.acc + 1
  #     params$a <- matrix(exp(HMCout$q[1:nkt]), nknots, nt)
  #     params$alpha <- transform$inv.logit(tail(HMCout$q, 1))
  #     calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
  #     calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  #   }
  
  #   if (a_alpha.att > 100) {
  #     a_alpha.rate <- a_alpha.acc / a_alpha.att
  #     if (a_alpha.rate < 0.20) {
  #       a_alpha.eps <- a_alpha.eps * 0.8
  #     } else if (a_alpha.rate > 0.60) {
  #       a_alpha.eps <- a_alpha.eps * 1.2
  #     }
  #     a_alpha.acc <- a_alpha.att <- 0
  #   }
  
  q <- transform$logit(b$cur)
  b$att <- b$att + 1
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q, epsilon = b$eps, 
                 L = 10, 
                 data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                 rho = rho, calc = calc, others = others, this.para = "b")
  if (HMCout$accept) {
    b$acc <- b$acc + 1
    b$cur <- transform$inv.logit(HMCout$q)
  }
  
  #   if (b.att > 100) {
  #     b.rate <- b.acc / b.att
  #     if (b.rate < 0.20) {
  #       b.eps <- b.eps * 0.8
  #     } else if (b.rate > 0.60) {
  #       b.eps <- b.eps * 1.2
  #     }
  #     b.acc <- b.att <- 0
  #   }
  
  storage.a[i, , ] <- a$cur
  storage.b[i, , ] <- b$cur
  storage.alpha[i] <- alpha$cur
  storage.beta[i]  <- beta$cur
  storage.prob[i, , ] <- 1 - exp(-calc$theta)
  
  if (i %% 500 == 0) {
    start <- max(i - 5000, 1)
    end   <- i
    par(mfrow=c(4, 5))
    plot.idx <- seq(1, 18, by = 2)
    for (idx in plot.idx){
      plot(log(storage.a[1:i, idx, 1]), type = "l", 
           main = round(log(gen$a[idx, 1]), 2), 
           xlab = round(a$acc / a$att, 3))
    }
    plot.idx <- seq(1, 18, by = 2)
    for (idx in plot.idx){
      plot(storage.b[1:i, idx, 1], type = "l", 
           xlab = round(b$acc / b$att, 3))
    }
    #     plot.idx <- 1:18
    #     for (idx in plot.idx){
    #       plot(storage.prob[start:end, idx, 1], type = "l")
    #     }
    plot(storage.beta[start:end], type = "l", 
         xlab = round(beta$acc / beta$att, 2), main = round(beta$eps, 3))
    plot(storage.alpha[start:end], type = "l",
         xlab = round(alpha$acc / alpha$att, 2))
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()


# try with the update for beta (HMC), a, alpha (higher rate of occurrence)
rm(list=ls())
source("./hmc_aux.R")
source("./updateModel.R")
source("./auxfunctions.R")
options(warn = 2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1000
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
nkt <- nknots * nt
rho.t <- list(cur = 0.25)
alpha.t <- list(cur = 0.25)
x <- matrix(1, ns, nt)
s <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0, 1, length = nknotsx), seq(0, 1, length = nknotsy))
dw2 <- rdist(s, knots)

others <- list(A.cutoff = max(sqrt(dw2)), thresh = 0, dw2 = dw2)
data <- list(x = x, s = s, knots = knots)
calc.t <- list()
calc.t$w <- getW(rho = rho.t, others = others)
calc.t$z <- matrix(rgev(n = ns * nt, 1, 1, 1), ns, nt)
a <- list(cur = matrix(rPS(nknots * nt, alpha = alpha.t$cur), nknots, nt))
calc.t$aw <- getAW(alpha = alpha.t, a = a, calc = calc.t)
calc.t$theta <- getTheta(alpha = alpha.t, calc = calc.t)

gen <- rRareBinarySpat(x = x, s = s, knots = knots, beta = 0, xi = 0, alpha = alpha.t,
                       rho = rho.t, nt = 1, prob.success = 0.05, dw2 = dw2)

# create lists for MCMC
data   <- list(y = gen$y, x = x, s = s, knots = knots)
calc   <- list(w = calc.t$w)  # need aw, theta

# initial values
beta.init <- -log(-log(mean(data$y)))
beta  <- list(cur = beta.init, att = 0, acc = 0, eps = 0.1, mn = 0, sd = 100)
xi    <- list(cur = 0, att = 0, acc = 0, eps = 0.01, mn = 0, sd = 0.5)
a     <- list(cur = matrix(10, nknots, nt), att = 0, acc = 0, eps = 0.3)
b     <- list(cur = matrix(0.5, nknots, nt), att = 0, acc = 0, eps = 0.3)
alpha <- list(cur = 0.5, att = 0, acc = 0, eps = 0.005)
rho   <- list(cur = 0.1, att = 0, acc = 0, eps = 0.01, mn = 0, sd = 1)

calc$x.beta <- getXBeta(data = data, beta = beta)
calc$z      <- getZ(xi = xi, calc = calc, others = others)
calc$aw     <- getAW(alpha = alpha, a = a, calc = calc)
calc$theta  <- getTheta(alpha = alpha, calc = calc)

library(numDeriv)
neg_log_post_a(q = log(a$cur), data = data, beta = beta, xi = xi, 
               a = a, b = b, alpha = alpha, rho = rho, calc = calc, others = others)
neg_log_post_grad_a(q = log(a$cur), data = data, beta = beta, xi = xi, a = a, 
                    b = b, alpha = alpha, rho = rho, calc = calc, others = others)
neg_log_post_grad_a_alpha(q = as.vector(c(log(a$cur), transform$logit(alpha$cur))),
                          data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha,
                          rho = rho, calc = calc, others = others, eps = 0.0001)

grad(func = neg_log_post_a_alpha, 
     x = as.vector(c(log(a$cur), transform$logit(alpha$cur))), 
     data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, rho = rho, calc = calc,
     others = others)
neg_log_post_alpha(q = transform$logit(alpha$cur), data = data, beta = beta, xi = xi, 
                   a = a, b = b, alpha = alpha, rho = rho, calc = calc, others = others)
neg_log_post_grad_alpha(q = transform$logit(alpha$cur), data = data, beta = beta, xi = xi, 
                        a = a, b = b, alpha = alpha, rho = rho, calc = calc, others = others, 
                        eps = 0.0001)
grad(func = neg_log_post_alpha, x = transform$logit(alpha$cur), 
     data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, rho = rho, 
     calc = calc, others = others)

neg_log_post_beta(q = beta$cur, data = data, beta = beta, xi = xi, a = a, 
                  b = b, alpha = alpha, rho = rho, calc = calc, others = others)
neg_log_post_grad_beta(q = beta$cur, data = data, beta = beta, xi = xi, a = a,
                       b = b, alpha = alpha, rho = rho, calc = calc, 
                       others = others, eps = 0.0001)
grad(func = neg_log_post_beta, x = beta$cur, data = data, beta = beta, xi = xi,
     a = a, b = b, alpha = alpha, rho = rho, calc = calc, others = others)


niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
storage.beta  <- rep(NA, niters)
storage.prob  <- array(NA, dim = c(niters, ns, nt))

set.seed(200)
for (i in 1:niters) {
  beta$att <- beta$att + 1
  q <- beta$cur
  HMCout <- HMC(neg_log_post_beta, neg_log_post_grad_beta, q, 
                epsilon = beta$eps, L = 10, 
                data = data, beta = beta, xi = xi, a = a, b =b, alpha = alpha, 
                calc = calc, others = others, this.param = "beta")
  if (HMCout$accept) {
    beta$acc    <- beta$acc + 1
    beta$cur    <- HMCout$q
    calc$x.beta <- getXBeta(data = data, beta = beta)
    calc$z      <- getZ(xi = xi, calc = calc, others = others)
    calc$theta  <- getTheta(alpha = alpha, calc = calc)
  }
  
#   if (beta$att > 100) {
#     beta.rate <- beta$acc / beta$att
#     if (beta.rate < 0.20) {
#       beta$eps <- beta$eps * 0.8
#     } else if (beta.rate > 0.60) {
#       beta$eps <- beta$eps * 1.2
#     }
#     beta$acc <- beta$att <- 0
#   }
  
  a$att <- a$att + 1
  q <- log(a$cur)
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q, 
                 epsilon = a$eps, L = 30, 
                 data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                 rho = rho, calc = calc, others = others, this.param = "a")
  if (HMCout$accept) {
    a$acc <- a$acc + 1
    a$cur <- exp(HMCout$q)
    calc$aw  <- getAW(alpha = alpha, a = a, calc = calc)
    calc$theta <- getTheta(alpha = alpha, calc = calc)
  }
  
  #   if (a.att > 100) {
  #     a.rate <- a.acc / a.att
  #     if (a.rate < 0.20) {
  #       a.eps <- a.eps * 0.8
  #     } else if (a.rate > 0.60) {
  #       a.eps <- a.eps * 1.2
  #     }
  #     a.acc <- a.att <- 0
  #   }
  
  alpha$att <- alpha$att + 1
  q <- transform$logit(alpha$cur)
  HMCout  <- HMC(neg_log_post_alpha, neg_log_post_grad_alpha, q, 
                 epsilon = alpha$eps, L = 30, 
                 data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                 rho = rho, calc = calc, others = others, this.param = "alpha")
  if (HMCout$accept) {
    alpha$acc <- alpha$acc + 1
    alpha$cur <- transform$inv.logit(HMCout$q)
    calc$aw  <- getAW(alpha = alpha, a = a, calc = calc)
    calc$theta <- getTheta(alpha = alpha, calc = calc)
  }
  
  #   if (alpha.att > 500) {
  #     alpha.rate <- alpha.acc / alpha.att
  #     if (alpha.rate < 0.20) {
  #       alpha.eps <- alpha.eps * 0.8
  #     } else if (alpha.rate > 0.60) {
  #       alpha.eps <- alpha.eps * 1.2
  #     }
  #     alpha.acc <- alpha.att <- 0
  #   }
  
  #   a_alpha.att <- a_alpha.att + 1
  #   q <- as.vector(c(log(params$a), transform$logit(params$alpha)))
  #   HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q, 
  #                  epsilon = a_alpha.eps, L = 30, 
  #                  data = data, params = params, calc = calc, others = others, 
  #                  prior = prior, this.param = "a and alpha")
  #   if (HMCout$accept) {
  #     a_alpha.acc <- a_alpha.acc + 1
  #     params$a <- matrix(exp(HMCout$q[1:nkt]), nknots, nt)
  #     params$alpha <- transform$inv.logit(tail(HMCout$q, 1))
  #     calc$aw  <- getAW(d = data, p = params, c = calc, o = others)
  #     calc$theta <- getTheta(d = data, p = params, c = calc, o = others)
  #   }
  
  #   if (a_alpha.att > 100) {
  #     a_alpha.rate <- a_alpha.acc / a_alpha.att
  #     if (a_alpha.rate < 0.20) {
  #       a_alpha.eps <- a_alpha.eps * 0.8
  #     } else if (a_alpha.rate > 0.60) {
  #       a_alpha.eps <- a_alpha.eps * 1.2
  #     }
  #     a_alpha.acc <- a_alpha.att <- 0
  #   }
  
  q <- transform$logit(b$cur)
  b$att <- b$att + 1
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q, epsilon = b$eps, 
                 L = 10, 
                 data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                 rho = rho, calc = calc, others = others, this.para = "b")
  if (HMCout$accept) {
    b$acc <- b$acc + 1
    b$cur <- transform$inv.logit(HMCout$q)
  }
  
  #   if (b.att > 100) {
  #     b.rate <- b.acc / b.att
  #     if (b.rate < 0.20) {
  #       b.eps <- b.eps * 0.8
  #     } else if (b.rate > 0.60) {
  #       b.eps <- b.eps * 1.2
  #     }
  #     b.acc <- b.att <- 0
  #   }
  
  storage.a[i, , ] <- a$cur
  storage.b[i, , ] <- b$cur
  storage.alpha[i] <- alpha$cur
  storage.beta[i]  <- beta$cur
  storage.prob[i, , ] <- 1 - exp(-calc$theta)
  
  if (i %% 500 == 0) {
    start <- max(i - 5000, 1)
    end   <- i
    par(mfrow=c(4, 5))
    plot.idx <- seq(1, 18, by = 2)
    for (idx in plot.idx){
      plot(log(storage.a[1:i, idx, 1]), type = "l", 
           main = round(log(gen$a[idx, 1]), 2), 
           xlab = round(a$acc / a$att, 3))
    }
    plot.idx <- seq(1, 18, by = 2)
    for (idx in plot.idx){
      plot(storage.b[1:i, idx, 1], type = "l", 
           xlab = round(b$acc / b$att, 3))
    }
    #     plot.idx <- 1:18
    #     for (idx in plot.idx){
    #       plot(storage.prob[start:end, idx, 1], type = "l")
    #     }
    plot(storage.beta[start:end], type = "l", 
         xlab = round(beta$acc / beta$att, 2), main = round(beta$eps, 3))
    plot(storage.alpha[start:end], type = "l",
         xlab = round(alpha$acc / alpha$att, 2))
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()


# checking to make sure that the gradients are coming back near what we expect from the 
# full log posterior

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
nknotsx <- 2
nknotsy <- 3
nknots <- nknotsx * nknotsy
nkt <- nknots * nt
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

gen <- rRareBinarySpat(x = x, s = s, knots = knots, beta = 0, xi = 0, alpha = alpha.t,
                       rho = rho.t, nt = 1, prob.success = 0.05, dw2 = dw2)

# create lists for MCMC
data   <- list(y = gen$y, x = x, s = s, knots = knots, dw2 = dw2)
params <- list(beta = 0, xi = 0, rho = rho.t, alpha = alpha.t)  # still need a and b
calc   <- list(z = calc.t$z, x.beta = 0, w = calc.t$w)  # need aw, theta
prior  <- list(beta.mn = 0, beta.sd = 100)

# initial values
a <- matrix(1, nknots, nt)
b <- matrix(0.5, nknots, nt)
alpha <- 0.5
beta  <- 0

params$a <- a
params$b <- b
params$alpha <- alpha
params$beta <- beta

calc$x.beta <- getXBeta(d = data, p = params, c = calc, o = others)
calc$z      <- getZ(d = data, p = params, c = calc, o = others)
calc$aw     <- getAW(d = data, p = params, c = calc, o = others)
calc$theta  <- getTheta(d = data, p = params, c = calc, o = others)


neg_log_post_beta_full(q = 0, d = data, p = params, c = calc, o = others, prior = prior)
grad(neg_log_post_beta_full, x = 0, d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_beta(q = 0, d = data, p = params, c = calc, o = others, prior = prior)

neg_log_post_a_full(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
grad(neg_log_post_a_full, x = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a2(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)

neg_log_post_grad_a(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_a2(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)

microbenchmark(neg_log_post_grad_a(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior), 
               neg_log_post_grad_a2(q = log(params$a), d = data, p = params, c = calc, o = others, prior = prior))

grad(neg_log_post_a, x = log(params$a), d = data, p = params, c = calc, o = others, prior = prior)

neg_log_post_alpha_full(q = 0, d = data, p = params, c = calc, o = others, prior = prior)
grad(neg_log_post_alpha_full, x = 0, d = data, p = params, c = calc, o = others, prior = prior)
neg_log_post_grad_alpha(q = 0, d = data, p = params, c = calc, o = others, prior = prior)

microbenchmark(grad(neg_log_post_alpha_full, x = 0, d = data, p = params, c = calc, o = others, prior = prior),
               neg_log_post_grad_alpha(q = 0, d = data, p = params, c = calc, o = others, prior = prior))