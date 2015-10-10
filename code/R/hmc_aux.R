rm(list=ls())
library(Rcpp)
source("./HMC.R")
source("./auxfunctions.R")

neg_log_post_a <- function(q, others) {
  # q: log(a)
  # others is a list
  #   y: data
  #   alpha: spatial dependence
  #   wz: kernel weights
  
  # extract from the list
  y <- others$y
  alpha <- others$alpha
  
  nt <- ncol(y)
  nknots <- length(q)
  a  <- exp(q)
  b  <- transform$inv.logit(others$b)
  
  alpha1m <- 1 - alpha
  
  # start with the log prior
  # Remember: q = log(a)
  lc <- logc(b = b, alpha = alpha)
  ll <- sum(alpha / alpha1m * q - exp(lc) * a^(-alpha / alpha1m))
  
  # log prob of success for likelihood
  theta <- getThetaCPP(wz = others$wz, a_star = a^alpha, alpha = alpha)
  
  # add in the log likelihood
  ll <- ll + sum(dbinom(others$y, size = 1, prob = -expm1(-theta), log = TRUE))
  
  return (-ll)
}


neg_log_post_b <- function(q, others) {
  # q: logit(b)
  # others is a list
  #   y: data
  #   alpha: spatial dependence
  #   wz: kernel weights
  
  # extract from the list
  y <- others$y
  alpha <- others$alpha
  
  nknots <- length(q)
  a  <- exp(others$a)
  b  <- transform$inv.logit(q)
  
  alpha1m <- 1 - alpha
  
  # start with the log prior
  lc <- logc(b = b, alpha = alpha)
  ll <- sum(lc - exp(lc) * a^(-alpha / alpha1m))
  
  return (-ll)
}

# separate a and b
neg_log_post <- function(q, others) {
  # q contains the current values for both a and b
  #   a: log(positive stable random effects)
  #   b: unif(0, 1) aux variable
  #   first half of q is a, second half is b
  # others is a list
  #   y: data
  #   alpha: spatial dependence
  #   wz: kernel weights
  
  # extract from the list
  y <- others$y
  alpha <- others$alpha
  
  nt <- ncol(y)
  nknots <- length(q) / (2 * nt)
  la <- q[1:nknots, , drop = FALSE]
  a  <- exp(la)
  bstar <- q[(nknots + 1):(2 * nknots), drop = FALSE]
  b  <- transform$inv.logit(bstar)
  
  alpha1m <- 1 - alpha
  
  # start with the log prior
  lc <- logc(b = b, alpha = alpha)
  ll <- sum(alpha / alpha1m * la + lc - exp(lc) * a^(-alpha / alpha1m))
  
  # log prob of success for likelihood
  theta <- getThetaCPP(wz = others$wz, a_star = a^alpha, alpha = alpha)
  
  # add in the log likelihood
  ll <- ll + sum(dbinom(others$y, size = 1, prob = -expm1(-theta), log = TRUE))
  
  return (-ll)
}

logc <- function(b, alpha) {
  alpha1m <- 1 - alpha  # used a few times below
  results <- alpha * log(sin(alpha * pi * b)) / (alpha1m) - 
             log(sin(pi * b)) / (alpha1m) + 
             log(sin(alpha1m) * pi * b)
  
  return (results)
}

neg_log_post_grad_a <- function(q, others) {
  # q: log(a)
  # others is a list
  #   y: data
  #   alpha: spatial dependence
  #   wz: kernel weights
  
  # extract things from others list if they're used more than once
  alpha <- others$alpha
  wz    <- others$wz
  y     <- others$y
  
  nt <- ncol(y)
  nknots <- length(q)
  a  <- exp(q)
  
  # get calculated values
  lc      <- logc(b = others$b, alpha = alpha)
  theta   <- getThetaCPP(wz = wz, a_star = a^alpha, alpha = alpha)
  alpha1m <- 1 - alpha
  
  # grad wrt log(a)
  grad <- matrix(0, nrow(a), ncol(a))
  for (t in 1:nt) {
    for (l in 1:nknots) {
      these <- which(y[, t] == 1)
      grad[l, t] <- a[l, t] * (sum((y[!these, t] - 1) * wz[!these, l, t]) + 
                                   sum(y[these, t] * wz[these, l, t] / expm1(-theta[these, t])))
    }
  }
  grad <- grad + alpha / alpha1m * (1 + exp(lc) * a^(-alpha / alpha1m))
  
  return(-grad)
}

neg_log_post_grad_b <- function(q, others) {
  # q: logit(b)
  # others is a list
  #   y: data
  #   alpha: spatial dependence
  #   wz: kernel weights
  
  # extract things from others list if they're used more than once
  alpha <- others$alpha
  a  <- exp(others$a)
  b  <- transform$inv.logit(q)
  
  # get calculated values
  lc      <- logc(b = b, alpha = alpha)
  alpha1m <- 1 - alpha  # used a lot in gradiant for b
  
  # grad wrt b
  apb   <- alpha * pi * b
  a1mpb <- alpha1m * pi * b
  pb    <- pi * b
  
  grad <- pi * (alpha^2 / (alpha1m * tan(apb)) - 1 / (alpha1m * tan(pb)) + 
          (alpha1m) / tan(a1mpb)) - a^(-alpha / alpha1m) * 
          (sin(apb)^(alpha / alpha1m) * sin(pb)^(-1 / alpha1m) * 
            (alpha1m * pi * cos(a1mpb) - pi / alpha1m * sin(a1mpb) / tan(pb)) + 
            sin(pb)^(-1 / alpha1m) * sin(a1mpb) * alpha^2 * pi / alpha1m * 
            cos(apb) * sin(apb)^(-1 / alpha1m)
          )
  
  
  return(-grad)
}

neg_log_post_grad <- function(q, others) {
  # q contains the current values in unconstrained space for both a and b
  #   a: log(positive stable random effects)
  #   b: logit(unif(0, 1)) aux variable
  #   first half of q is a, second half is b
  # others is a list
  #   y: data
  #   alpha: spatial dependence
  #   wz: kernel weights
  
  # extract things from others list if they're used more than once
  alpha <- others$alpha
  wz    <- others$wz
  y     <- others$y
  
  nt <- ncol(y)
  nknots <- length(q) / (2 * nt)
  la <- q[1:nknots, , drop = FALSE]
  a  <- exp(la)
  bstar <- q[(nknots + 1):(2 * nknots), , drop = FALSE]
  b  <- transform$inv.logit(bstar)
  
  # get calculated values
  lc      <- logc(b = b, alpha = alpha)
  theta   <- getThetaCPP(wz = wz, a_star = a^alpha, alpha = alpha)
  alpha1m <- 1 - alpha  # used a lot in gradiant for b
  
  # storage
  grad <- matrix(NA, nrow(q), ncol(q))
  
  # grad wrt a
  grad_a <- matrix(0, nrow(a), ncol(a))
  for (t in 1:nt) {
    for (l in 1:nknots) {
      these <- which(y[, t] == 1)
      grad_a[l, t] <- a[l, t] * (sum((y[!these, t] - 1) * wz[!these, l, t]) + 
                                 sum(y[these, t] * wz[these, l, t] / expm1(-theta[these, t])))
    }
  }
  grad_a <- grad_a + alpha / alpha1m * (1 + exp(lc) * a^(-alpha / alpha1m))
  
  # grad wrt b
  apb   <- alpha * pi * b
  a1mpb <- alpha1m * pi * b
  pb    <- pi * b
  
  grad_b <- pi * (alpha^2 / (alpha1m * tan(apb)) - 1 / (alpha1m * tan(pb)) + 
            (alpha1m) / tan(a1mpb)) - a^(-alpha / alpha1m) * 
            (sin(apb)^(alpha / alpha1m) * sin(pb)^(-1 / alpha1m) * 
              (alpha1m * pi * cos(a1mpb) - pi / alpha1m * sin(a1mpb) / tan(pb)) + 
              sin(pb)^(-1 / alpha1m) * sin(a1mpb) * alpha^2 * pi / alpha1m * 
              cos(apb) * sin(apb)^(-1 / alpha1m)
            )
  
  grad[1:nknots, ] <- -grad_a
  grad[(nknots + 1):(2 * nknots), ] <- -grad_b
  
  return(grad)
}

# Test out the functions
library(fields)
library(evd)
ns <- 20
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
rho.t <- 0.25
alpha.t <- 0.3
s <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0, 1, length = nknotsx), seq(0, 1, length = nknotsy))
dw2 <- rdist(s, knots)
w.t <- stdW(makeW(dw2 = dw2, rho = rho.t))
z.t <- matrix(rgev(n = ns * nt, 1, 1, 1), ns, nt)
a.t <- matrix(rPS(nknots * nt, alpha = alpha.t), nknots, nt)
wz.t <- getwzCPP(z = z.t, w = w.t)
theta.t <- getThetaCPP(wz= wz.t, a_star = a.t^alpha.t, alpha = alpha.t)
y.t <- matrix(rbinom(ns * nt, size = 1, prob = -expm1(-theta.t)), ns, nt) 

b.t <- matrix(runif(nknots * nt), nknots, nt)
q <- rbind(log(a.t), transform$logit(b.t))

others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = b.t, a = a.t)
neg_log_post(q, others)

q.a <- log(a.t)
neg_log_post_grad_a(q.a, others)

q.b <- transform$logit(b.t)
neg_log_post_grad_b(q.b, others)

niters <- 5000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
q.a <- matrix(log(100), nknots, nt)
q.b <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = q.b, a = q.a)
for (i in 1:niters) {
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q.a, epsilon=0.001, L=10, others)
  if (HMCout$accept) {
    q.a      <- HMCout$q
    others$a <- HMCout$q
  }
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q.b, epsilon=0.001, L=10, others)
  if (HMCout$accept) {
    others.b <- HMCout$q
    q.b <- HMCout$q
  }
  storage.a[i, , ] <- others$a
  storage.b[i, , ] <- others$b
  if (i %% 500 == 0) {
    par(mfrow=c(3, 4))
    plot(storage.b[1:i, 1, 1], type = "l")
    plot(storage.b[1:i, 3, 1], type = "l")
    plot(storage.b[1:i, 5, 1], type = "l")
    plot(storage.b[1:i, 7, 1], type = "l")
    plot(storage.b[1:i, 9, 1], type = "l")
    plot(storage.b[1:i, 11, 1], type = "l")
    plot(storage.b[1:i, 13, 1], type = "l")
    plot(storage.b[1:i, 15, 1], type = "l")
    plot(storage.b[1:i, 17, 1], type = "l")
    plot(storage.b[1:i, 21, 1], type = "l")
    plot(storage.b[1:i, 23, 1], type = "l")
    plot(storage.b[1:i, 25, 1], type = "l")
  }
}