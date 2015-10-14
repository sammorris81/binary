library(Rcpp)
source("./HMC.R")
source("./auxfunctions.R")

neg_log_post_a <- function(q, others) {
  # q: log(a)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  a <- exp(q)
  
  # extract from the list and get calculated quantities
  y       <- others$y
  alpha   <- others$alpha
  b       <- others$b
  alpha1m <- 1 - alpha
  
  # log prior
  # Remember: q = log(a) and jacobian is included
  lc <- logc(b = b, alpha = alpha)
  ll <- sum(alpha / alpha1m * q - exp(lc) * a^(-alpha / alpha1m))
  
  # data and log likelihood
  theta <- getThetaCPP(wz = others$wz, a_star = a^alpha, alpha = alpha)
  ll <- ll + sum(logLikeY(y = y, theta = theta))
  
  return (-ll)
}


neg_log_post_b <- function(q, others) {
  # q: logit(b)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  b  <- transform$inv.logit(q)
  
  # extract from the list and get calculated quantities
  a       <- others$a
  alpha   <- others$alpha
  alpha1m <- 1 - alpha
  
  # log prior
  lc <- logc(b = b, alpha = alpha)
  ll <- sum(lc - exp(lc) * a^(-alpha / alpha1m) + 
        2 * log(b) - q)  # jacobian
  
  # b does not have an impact on the likelihood
  
  return (-ll)
}

neg_log_post_alpha <- function(q, others) {
  # q: logit(alpha)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  alpha   <- transform$inv.logit(q)
  alpha1m <- 1 - alpha
  
  # extract from the list and get calculated quantities
  y  <- others$y
  a  <- others$a
  b  <- others$b
  lc <- logc(b = b, alpha = alpha)
  
  nt     <- ncol(a)
  nknots <- nrow(a)
  
  # log density for PS random variable
  ll <- nknots * (log(alpha) - log(alpha1m)) + 
        2 * log(alpha) - q  # Jacobian
  ll <- ll + sum(alpha / alpha1m * log(a) + lc - 
                 exp(lc) * a^(-alpha / alpha1m))
  
  # data and log likelihood
  theta <- getThetaCPP(wz = others$wz, a_star = a^alpha, alpha = alpha)
  ll    <- ll + sum(logLikeY(y = y, theta = theta))
  
  return(-ll)
}

neg_log_post_a_alpha <- function(q, others) {
  # q: c(log(a), logit(alpha))
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  nknots  <- nrow(others$a)
  nt      <- ncol(others$a)
  q.a     <- matrix(q[1:(nknots * nt)], nknots, nt)
  q.alpha <- tail(q, 1)
  a       <- exp(q.a)
  alpha   <- transform$inv.logit(q.alpha)
  alpha1m <- 1 - alpha
  
  # extract from the list and get calculated quantities
  y       <- others$y
  b       <- others$b
  lc      <- logc(b = b, alpha = alpha)
  
  # log prior
  # Remember: q = log(a) and jacobian is included
  
  ll <- sum(alpha / alpha1m * q.a + lc - exp(lc) * a^(-alpha / alpha1m)) +
        nknots * (log(alpha) - log(alpha1m)) + 2 * log(alpha) - q.alpha
  
  # data and log likelihood
  theta <- getThetaCPP(wz = others$wz, a_star = a^alpha, alpha = alpha)
  ll <- ll + sum(logLikeY(y = y, theta = theta))

  return (-ll)
}

logc <- function(b, alpha) {
  alpha1m <- 1 - alpha
  apb <- alpha * pi * b
  pb  <- pi * b
  a1mpb <- alpha1m * pi * b 
  u <- sin(apb) / sin(pb)
  results <- log(u) / alpha1m + log(sin(a1mpb)) - log(sin(apb))
  
  return(results)
}

neg_log_post_grad_a <- function(q, others) {
  # q: log(a)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  a      <- exp(q)
  
  nt     <- ncol(a)
  nknots <- nrow(a)
  
  # extract from the list and get calculated quantities
  alpha <- others$alpha
  wz    <- others$wz
  y     <- others$y
  b     <- others$b
  
  # get calculated values
  lc      <- logc(b = b, alpha = alpha)
  theta   <- getThetaCPP(wz = wz, a_star = a^alpha, alpha = alpha)
  alpha1m <- 1 - alpha
  wz.star <- wz^(1 / alpha)
  
  # grad wrt log(a)
  grad <- matrix(0, nknots, nt)
  
  # the likelihood component should be
  # For non-observances:
  #   a.lt * wz.star.lt
  # For observances:
  #   a.lt * wz.star.lt / expm1(theta.t)
  for (t in 1:nt) {
    for (l in 1:nknots) {
      these <- which(y[, t] == 1)
      grad[l, t] <- a[l, t] * (sum(wz.star[y == 0, l, t]) - 
                               sum(wz.star[y == 1, l, t] / expm1(theta[y == 1, t])))
    }
  }
  grad <- grad - alpha / alpha1m * (1 + exp(lc) * a^(-alpha / alpha1m))

  return(grad)
}

neg_log_post_grad_b <- function(q, others) {
  # q: logit(b)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  b  <- transform$inv.logit(q)
  
  # extract from the list and get calculated quantities
  alpha   <- others$alpha
  alpha1m <- 1 - alpha  # used a lot in gradiant for b
  a       <- others$a
  apb     <- alpha * pi * b
  a1mpb   <- alpha1m * pi * b
  pb      <- pi * b
  u       <- (sin(apb) / sin(pb))
  
  # jacobian of b wrt q
  dbdq    <- b^2 * exp(-q)
  
  # dUdb
  grad <- alpha * pi / (alpha1m * tan(apb)) - u * pi * cos(pb) / (alpha1m * sin(apb)) +
          alpha1m * pi / tan(a1mpb) - alpha * pi / tan(apb)
  
  grad <- -grad + u^(1 / alpha1m) * a^(-alpha / alpha1m) * (
     (alpha * pi * cos(apb) - u * pi * cos(pb)) * sin(a1mpb) / (alpha1m * sin(apb)) +
       alpha1m * pi * cos(a1mpb) - alpha * pi * sin(a1mpb) / tan(apb)
   ) / sin(apb)
  
  # dUdq = dUdb * dbdq + d(log(dbdq)) / dq
  grad <- grad * dbdq - (2 * b * exp(-q) - 1)
  
  return(grad)
}

neg_log_post_grad_alpha <- function(q, others, eps = 0.0001) {
  # q: logit(alpha)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # using a quick and easy numerical approximation
  # might be good to eventually get the actual value, but this is pretty reliable
  
  grad <- (neg_log_post_alpha(q = q + 0.0001, others = others) - 
           neg_log_post_alpha(q = q, others = others)) / 0.0001
  
  return(grad)
}

neg_log_post_grad_a_alpha <- function(q, others) {
  # q is a list: 
  #   a:     log(a)
  #   alpha: logit(alpha))
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  nknots  <- length(q) - 1
  nt      <- ncol(others$y)
  q.a     <- matrix(q[1:nknots], nknots, nt)  # a needs to be a matrix
  q.alpha <- tail(q, 1)
  
  grad <- rep(NA, length(q))
  
  others$a <- exp(q.a)
  others$alpha <- transform$inv.logit(q.alpha)
  
  grad[1:nknots] <- neg_log_post_grad_a(q = q.a, others = others)
  grad[(nknots + 1)] <- neg_log_post_grad_alpha(q = q.alpha, others = others)
  
  return(grad)
}
