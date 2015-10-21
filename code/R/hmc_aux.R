library(Rcpp)
source("./HMC.R")

log_ps_full <- function(d, p, c, o, prior) {
  alpha1m <- 1 - p$alpha
  lc <- logc(b = p$b, alpha = p$alpha) 
  ll <- sum(log(p$alpha) - log(alpha1m) - 1 / alpha1m * log(p$a) + 
              lc - exp(lc) * p$a^(-p$alpha / alpha1m) + 
              log(p$a) + log(p$b) + log(1 - p$b))  # Jacobian
  
  return (ll)
}

log_post_full <- function(d, p, c, o, prior) {
  alpha1m <- 1 - p$alpha
  lc <- logc(b = p$b, alpha = p$alpha)
  
  # redo calculated
  c$x.beta <- getXBeta(d = d, p = p, c = c, o = o)
  c$z      <- getZ(d = d, p = p, c = c, o = o)
  c$w      <- getW(d = d, p = p, c = c, o = o)
  c$aw     <- getAW(d = d, p = p, c = c, o = o)
  c$theta  <- getTheta(d = d, p = p, c = c, o = o)
  
#   ll <- log_ps_full(d = d, p = p, c = c, o = o) + 
#         sum(logLikeY(d = d, p = p, c = c, o = o)) +
#         sum(dnorm(p$beta, prior$beta.mn, prior$beta.sd, log = TRUE))
  
  ll <- sum(log(p$alpha) - log(alpha1m) - log(p$a) / alpha1m + lc - exp(lc) * p$a^(-p$alpha / alpha1m) +
              log(p$a) + log(p$b) + log(1 - p$b))
  ll <- ll + sum(-c$theta[d$y == 0]) + sum(log(1 - exp(-c$theta[d$y == 1])))
  ll <- ll + sum(dnorm(p$beta, prior$beta.mn, prior$beta.sd, log = TRUE))
  
  return(ll)
}

neg_log_post_beta_full <- function(q, d, p, c, o, prior) {
  p$beta <- q
  c$x.beta <- getXBeta(d = d, p = p, c = c, o = o)
  c$z <- getZ(d = d, p = p, c = c, o = o)
  c$theta <- getTheta(d = d, p = p, c = c, o = o)
  ll <- log_post_full(d = d, p = p, c = c, o = o, prior = prior)
  return (-ll)
}

neg_log_post_a_full <- function(q, d, p, c, o, prior) {
  p$a <- exp(q)
  ll <- log_post_full(d = d, p = p, c = c, o = o, prior = prior)
  return (-ll)
}

neg_log_post_alpha_full <- function(q, d, p, c, o, prior) {
  p$alpha <- transform$inv.logit(q)
  ll <- log_post_full(d = d, p = p, c = c, o = o, prior = prior)
  return(-ll)
}

neg_log_post_beta <- function(q, d, p, c, o, prior) {
  # q: beta
  # others is a list
  #   y(ns, nt):      data
  #   x(ns, nt * np): covariates
  #   xi(1):          xi
  #   aw(ns, nt):     sum_l (A_l * w_li^(1 / alpha))
  #   alpha(1):       spatial dependence
  #   pri.mn(1):      prior mean
  #   pri.sd(1):      prior standard deviation
  
  # extract from the list and get calculated quantities
  p$beta <- q
  
  c$x.beta <- getXBeta(d = d, p = p, c = c, o = o)
  c$z      <- getZ(d = d, p = p, c = c, o = o)
  c$theta <- getTheta(d = d, p = p, c = c, o = o)
  
  # log prior
  ll <- sum(dnorm(x = p$beta, mean = prior$beta.mn, sd = prior$beta.sd, log = TRUE))
  
  
  ll    <- ll + sum(logLikeY(d = d, p = p, c = c, o = o))
  
  return(-ll)
}

neg_log_post_a <- function(q, d, p, c, o, prior) {
  # q: log(a)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  p$a <- exp(q)
  
  # extract from the list and get calculated quantities
  alpha1m <- 1 - p$alpha
  
  # log prior
  # Remember: q = log(a) and jacobian is included
  lc <- logc(b = p$b, alpha = p$alpha)
  ll <- sum(-p$alpha / alpha1m * q - exp(lc) * p$a^(-p$alpha / alpha1m))
  
  # data and log likelihood
  c$aw  <- getAW(d = d, p = p, c = c, o = o)
  c$theta <- getTheta(d = d, p = p, c = c, o = o)
  ll <- ll + sum(logLikeY(d = d, p = p, c = c, o = o))
  
  return (-ll)
}

neg_log_post_alpha <- function(q, d, p, c, o, prior) {
  # q: logit(alpha)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  p$alpha <- transform$inv.logit(q)
  c$aw  <- getAW(d = d, p = p, c = c, o = o)
  c$theta <- getTheta(d = d, p = p, c = c, o = o)
  alpha1m <- 1 - p$alpha
  
  # extract from the list and get calculated quantities
  lc <- logc(b = p$b, alpha = p$alpha)
  
  nt     <- ncol(p$a)
  nknots <- nrow(p$a)
  
  # log density for PS random variable
  ll <- sum(log(p$alpha) - log(alpha1m) - 1 / alpha1m * log(p$a) + lc - 
              exp(lc) * p$a^(-p$alpha / alpha1m))
  
  # data and log likelihood
  ll    <- ll + sum(logLikeY(d = d, p = p, c = c, o = o))
  
  # log prior
  # ll <- ll + sum(dnorm(q, 0, 0.5, log = TRUE))
  
  return(-ll)
}

neg_log_post_a_alpha <- function(q, d, p, c, o, prior) {
  # q: c(log(a), logit(alpha))
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  nknots  <- nrow(p$a)
  nt      <- ncol(p$a)
  q.a     <- matrix(q[1:(nknots * nt)], nknots, nt)
  q.alpha <- tail(q, 1)
  p$a     <- exp(q.a)
  p$alpha <- transform$inv.logit(q.alpha)
  alpha1m <- 1 - p$alpha
  
  # log prior
  # Remember: q = log(a) and jacobian is included
  lc      <- logc(b = p$b, alpha = p$alpha)
  ll <- sum(-p$alpha / alpha1m * q.a + lc - exp(lc) * p$a^(-p$alpha / alpha1m)) +
    nknots * nt * (log(p$alpha) - log(alpha1m))
  
  # data and log likelihood
  c$aw  <- getAW(d = d, p = p, c = c, o = o)
  c$theta <- getTheta(d = d, p = p, c = c, o = o)
  ll <- ll + sum(logLikeY(d = d, p = p, c = c, o = o))
  
  return (-ll)
}

llps <- function(q, d, p, c, o, prior, addup = TRUE) {
  # q: logit(b)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  p$b <- transform$inv.logit(q)
  
  # extract from the list and get calculated quantities
  alpha1m <- 1 - p$alpha
  
  # log prior
  lc <- logc(b = p$b, alpha = p$alpha)
  ll <- lc - exp(lc) * a^(-p$alpha / alpha1m) + log(p$b) + log(1 - p$b)
  
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

neg_log_post_b <- function(q, d, p, c, o, prior, addup = TRUE) {
  # q: logit(b)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  p$b <- transform$inv.logit(q)
  
  # extract from the list and get calculated quantities
  alpha1m <- 1 - p$alpha
  
  # prior is U(0, 1)
  lc <- logc(b = p$b, alpha = p$alpha)
  ll <- lc - exp(lc) * a^(-p$alpha / alpha1m) + 
        log(p$b) + log(1 - p$b)  # jacobian
  
  # b does not have an impact on the likelihood
  
  if (addup) {
    ll <- sum(ll)
  }
  
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

neg_log_post_grad_a <- function(q, d, p, c, o, prior) {
  # q: log(a)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable

  # parameter transformation
  p$a  <- exp(q)
  
  # for the gradiant of a_alpha, we calculate this before calling this function
  c$aw <- getAW(d = d, p = p, c = c, o = o)
  c$theta <- getTheta(d = d, p = p, c = c, o = o)
  
  nt     <- ncol(p$a)
  nknots <- nrow(p$a)
  
  # extract from the list and get calculated quantities
  lc      <- logc(b = p$b, alpha = p$alpha)
  alpha1m <- 1 - p$alpha
  
  # grad wrt log(a)
  grad <- matrix(0, nknots, nt)
  
  # the likelihood component should be
  # For non-observances:
  #   a.lt * wz.star.lt
  # For observances:
  #   a.lt * wz.star.lt / expm1(theta.t)
  for (t in 1:nt) {
    wz.star.t <- exp((log(c$w) - log(c$z[, t])) / p$alpha)
    for (l in 1:nknots) {
      grad[l, t] <- p$a[l, t] * (-sum(wz.star.t[d$y == 0, l]) + 
                                 sum(wz.star.t[d$y == 1, l] / expm1(c$theta[d$y == 1, t])))
    }
  }
  grad <- grad - p$alpha / alpha1m * (1 - exp(lc) * p$a^(-p$alpha / alpha1m))
  
  # whole function is written as gradient of log likelihood
  return(-grad)
}

# just to verify that the gradient is coming back as it should
# neg_log_post_grad_a2 <- function(q, d, p, c, o, prior) {
#   p$a <- exp(q)
#   nknots <- nrow(p$a)
#   nt     <- ncol(p$a)
#   ns     <- nrow(d$y)
#   
#   c$aw <- getAW(d = d, p = p, c = c, o = o)
#   c$theta <- getTheta(d = d, p = p, c = c, o = o)
#   
#   alpha1m <- 1 - p$alpha
#   lc      <- logc(b = p$b, alpha = p$alpha)
#   
#   grad <- matrix(0, nknots, nt) 
#   for (k in 1:nknots) { 
#     for (t in 1:nt) {
#       grad[k, t] <- grad[k, t] - p$alpha / alpha1m + 
#         p$alpha / alpha1m * exp(lc[k, t]) * p$a[k, t]^(-p$alpha / alpha1m) - 
#         p$a[k, t] * sum((c$w[d$y == 0, k] / c$z[d$y == 0, t])^(1 / p$alpha)) +
#         p$a[k, t] * sum((c$w[d$y == 1, k] / c$z[d$y == 1, t])^(1 / p$alpha) / (exp(c$theta[d$y == 1, t]) - 1))
#     }
#   }
#   
#   return (-grad)
# }

neg_log_post_grad_b <- function(q, d, p, c, o, prior) {
  # q: logit(b)
  # others is a list
  #   y:     data
  #   alpha: spatial dependence
  #   wz:    kernel weights
  #   a:     positive stable random effects
  #   b:     auxiliary random variable
  
  # parameter transformation
  p$b <- transform$inv.logit(q)
  
  # extract from the list and get calculated quantities
  alpha1m <- 1 - p$alpha  # used a lot in gradiant for b
  apb     <- p$alpha * pi * p$b
  a1mpb   <- alpha1m * pi * p$b
  pb      <- pi * p$b
  c       <- exp(logc(b = p$b, alpha = p$alpha))
  
  # chain rule work
  u <- sin(apb)^(p$alpha / alpha1m)
  v <- sin(pb)^(-1 / alpha1m)
  w <- sin(a1mpb)
  
  # chain rule derivatives
  dbdq <- p$b * (1 - p$b)  # Jacobian
  dudb <- p$alpha^2 * pi / alpha1m * sin(apb)^(p$alpha / alpha1m - 1) * cos(apb)
  dvdb <- - pi / alpha1m * sin(pb)^(-1 / alpha1m - 1) * cos(pb)
  dwdb <- cos(a1mpb) * alpha1m * pi
  dcdb <- u * v * dwdb + u * dvdb * w + dudb * v * w
  dcdq <- dcdb * dbdq
  
  # dcdb
  grad <- dcdq * (-1 / c + a^(-p$alpha / alpha1m)) + 2 * p$b - 1
  
  return(grad)
}

neg_log_post_grad_alpha <- function(q, d, p, c, o, prior, eps = 0.0001) {
  # q: logit(alpha)

  # using a quick and easy numerical approximation
  # might be good to eventually get the actual value, but this is pretty reliable
  # need to recalculate aw and theta when we add eps to q, otherwise it doesn't 
  #   take into account that there was some movement.

  grad <- (neg_log_post_alpha(q = q + eps, d = d, p = p, c = c, o = o, 
                              prior = prior) - 
           neg_log_post_alpha(q = q, d = d, p = p, c = c, o = o, 
                              prior = prior)) / eps
  
  return(grad)
}

neg_log_post_grad_beta <- function(q, d, p, c, o, prior, eps = 0.0001) {
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
  
  grad <- (neg_log_post_beta(q = q + eps, d = d, p = p, c = c, o = o, 
                             prior = prior) - 
           neg_log_post_beta(q = q, d = d, p = p, c = c, o = o, prior = prior)) / eps
  
  return(grad)
}

neg_log_post_grad_beta_2 <- function(q, d, p, c, o, prior) {
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
  p$beta <- q
  
  # chain rule - need to figure out how to deal with multiple days
  if (p$xi == 0) {
    dzdb <- exp(-c$x.beta) * (-d$x)
  }
  
  # logprior
  ll <- -2 * (p$beta - prior$beta.mn) / prior$beta.sd^2
  
  # data
  
  return(grad)
}

neg_log_post_grad_a_alpha <- function(q, d, p, c, o, prior) {
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
  nknots  <- ncol(c$w)# for the gradiant of a_alpha, we calculate this before calling this function
  nt      <- ncol(d$y)
  q.a     <- matrix(q[1:nknots], nknots, nt)  # a needs to be a matrix
  q.alpha <- tail(q, 1)
  
  grad <- rep(NA, length(q))
  
  p$a <- exp(q.a)
  p$alpha <- transform$inv.logit(q.alpha)
  c$aw <- getAW(d = d, p = p, c = c, o = o)
  c$theta <- getTheta(d = d, p = p, c = c, o = o)
  
  grad[1:nknots] <- neg_log_post_grad_a(q = q.a, d = d, p = p, c = c, o = o, 
                                        prior = prior)
  grad[(nknots + 1)] <- neg_log_post_grad_alpha(q = q.alpha, d = d, p = p, 
                                                c = c, o = o, prior = prior)
  
  return(grad)
}
