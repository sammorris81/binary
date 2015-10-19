library(Rcpp)
source("./HMC.R")

neg_log_post_beta <- function(q, others) {
  # q: beta
  # others is a list
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
  
  # extract from the list and get calculated quantities
  y     <- others$y
  alpha <- others$alpha
  w     <- others$w
  x     <- others$x
  ns <- nrow(y)
  nt <- ncol(y)
  
  x.beta <- getXBeta(x = x, ns = ns, nt = nt, beta = q)
  z      <- getZ(xi = xi, x.beta = x.beta, thresh = 0)
  wz     <- getwzCPP(z = z, w = w)
  
  # log prior
  ll <- dnorm(x = q, mean = others$pri.mn, sd = others$pri.sd, log = TRUE)
  
  theta <- getThetaCPP(wz = wz, a_star = a^alpha, alpha = alpha)
  ll    <- ll + sum(logLikeY(y = y, theta = theta))
  
  return(-ll)
}

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


llps <- function(q, others, addup = TRUE) {
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
  ll <- lc - exp(lc) * a^(-alpha / alpha1m) + log(b) + log(1 - b)
  
  if (addup) {
    ll <- sum(ll)
  }
  
  # b does not have an impact on the likelihood
}

# aux functions for b
#   c(b) = u * p * q
#   u: sin(alpha * pi * b)^(alpha / (1 - alpha))
#   p: sin(pi * b)^(-1 / (1 - alpha))
#   q: sin((1 - alpha) * pi * b)

neg_log_post_b <- function(q, others, addup = TRUE) {
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
  
  # prior is U(0, 1)
  lc <- logc(b = b, alpha = alpha)
  ll <- lc - exp(lc) * a^(-alpha / alpha1m) + 
        log(b) + log(1 - b)  # jacobian
  
  if (addup) {
    ll <- sum(ll)
  }
  
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
  ll <- nknots * (log(alpha) - log(alpha1m))
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
        nknots * (log(alpha) - log(alpha1m))
  
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

  results <- alpha * log(sin(apb)) / alpha1m - log(sin(pb)) / alpha1m +
             log(sin(a1mpb))
  
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
  c       <- exp(logc(b = b, alpha = alpha))
  
  # chain rule work
  u <- sin(apb)^(alpha / alpha1m)
  v <- sin(pb)^(-1 / alpha1m)
  w <- sin(a1mpb)
  
  # chain rule derivatives
  dbdq <- b * (1 - b)  # Jacobian
  dudb <- alpha^2 * pi / alpha1m * sin(apb)^(alpha / alpha1m - 1) * cos(apb)
  dvdb <- - pi / alpha1m * sin(pb)^(-1 / alpha1m - 1) * cos(pb)
  dwdb <- cos(a1mpb) * alpha1m * pi
  dcdb <- u * v * dwdb + u * dvdb * w + dudb * v * w
  dcdq <- dcdb * dbdq
  
  # dcdb
  grad <- dcdq * (-1 / c + a^(-alpha / alpha1m)) + 2 * b - 1
  
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
  
  grad <- (neg_log_post_alpha(q = q + eps, others = others) - 
           neg_log_post_alpha(q = q, others = others)) / eps
  
  return(grad)
}

neg_log_post_grad_beta <- function(q, others, eps = 0.0001) {
  # q: beta
  # others is a list
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
  
  # using a quick and easy numerical approximation
  # might be good to eventually get the actual value, but this is pretty reliable
  
  grad <- (neg_log_post_beta(q = q + eps, others = others) - 
           neg_log_post_beta(q = q, others = others)) / eps
  
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
