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
  ll <- sum(-alpha / alpha1m * q - exp(lc) * a^(-alpha / alpha1m))
  
  # log prob of success for likelihood
  theta <- getThetaCPP(wz = others$wz, a_star = a^alpha, alpha = alpha)
  
  # add in the log likelihood
  ll <- ll + sum(logLikeY(y = y, theta = theta))
  
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
  alpha1m <- 1 - alpha
  results <- alpha * log(sin(alpha * pi * b)) / alpha1m - 
             log(sin(pi * b)) / alpha1m + log(sin(alpha1m * pi * b))
  
  return(results)
}

neg_log_post_grad_a <- function(q, others) {
  # q: log(a)
  # others is a list
  #   y: data
  #   alpha: spatial dependence
  #   wz: kernel weights
  #   b: logit(b)
  
  # extract things from others list if they're used more than once
  alpha <- others$alpha
  wz    <- others$wz
  y     <- others$y
  b     <- transform$inv.logit(others$b)
  
  nt <- ncol(y)
  nknots <- length(q)
  a  <- exp(q)
  
  # get calculated values
  lc      <- logc(b = b, alpha = alpha)
  theta   <- getThetaCPP(wz = wz, a_star = a^alpha, alpha = alpha)
  alpha1m <- 1 - alpha
  
  # grad wrt log(a)
  grad <- matrix(0, nrow(a), ncol(a))
  wz.star <- wz^(1 / alpha)
  
  # the likelihood component should be
  # For non-observances:
  #   -(1 - y.i) * wz.star.lt * a.lt
  # For observances:
  #   -y.i * (wz.star.lt) * 1 / (expm1(theta.t)) * a.lt
  for (t in 1:nt) {
    for (l in 1:nknots) {
      these <- which(y[, t] == 1)
      grad[l, t] <- -a[l, t] * (sum((1 - y[!these, t]) * wz.star[!these, l, t]) + 
                                sum(y[these, t] * wz.star[these, l, t] / expm1(theta[these, t])))
    }
  }
  grad <- grad - alpha / alpha1m * (1 + exp(lc) * a^(-alpha / alpha1m))
  
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
  alpha1m <- 1 - alpha  # used a lot in grad
  
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

