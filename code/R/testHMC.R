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
  
  q <- as.vector(c(log(params$a), transform$logit(params$alpha)))
  HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q, epsilon = 0.01, L = 10, 
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