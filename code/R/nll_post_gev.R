library(Rcpp)

log_ps_full <- function(data, alpha, b, a, calc, others) {
  alpha1m <- 1 - alpha$cur
  lc <- logc(b = b$cur, alpha = alpha$cur) 
  ll <- sum(log(alpha$cur) - log(alpha1m) - 1 / alpha1m * log(a$cur) + 
              lc - exp(lc) * a$cur^(-alpha$cur / alpha1m) + 
              log(a$cur) + log(b$cur) + log(1 - b$cur))  # Jacobian
  
  return (ll)
}

log_post_full <- function(data, beta, xi, a, b, alpha, rho, calc, others) {
  alpha1m <- 1 - alpha$cur
  lc      <- logc(b = b$cur, alpha = alpha$cur)
  
  # redo calculated
  x.beta <- getXBeta(y = data$y, x = data$x, beta = beta$cur)
  z      <- getZ(xi = xi$cur, x.beta = x.beta, thresh = others$thresh)
  w      <- getW(rho = rho$cur, dw2 = others$dw2, a.cutoff = others$a.cutoff)
  w.star <- getWStarIDs(alpha = alpha$cur, w = w, IDs = others$IDs)
  aw     <- getAW(a = a$cur, w.star = w.star)
  theta  <- getTheta(alpha = alpha$cur, z = z, aw = aw)
  
#   ll <- log_ps_full(d = d, p = p, c = c, o = o) + 
#         sum(logLikeY(d = d, p = p, c = c, o = o)) +
#         sum(dnorm(p$beta, prior$beta.mn, prior$beta.sd, log = TRUE))
  
  ll <- sum(log(alpha$cur) - log(alpha1m) - log(a$cur) / alpha1m + lc - 
              exp(lc) * a$cur^(-alpha$cur / alpha1m) +
              log(a$cur) + log(b$cur) + log(1 - b$cur))
  ll <- ll + sum(-theta[data$y == 0]) + 
    sum(log(1 - exp(-theta[data$y == 1])))
  ll <- ll + sum(dnorm(beta$cur, beta$mn, beta$sd, log = TRUE))
  
  return(ll)
}

neg_log_post_beta_full <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                                   others) {
  beta$cur <- q
  ll <- log_post_full(data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                      rho = rho, calc = calc, others = others)
  return (-ll)
}

neg_log_post_a_full <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                                others) {
  a$cur <- exp(q)
  ll <- log_post_full(data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                      rho = rho, calc = calc, others = others)
  return (-ll)
}

neg_log_post_b_full <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                                others) {
  b$cur <- transform$inv.logit(q)
  ll <- log_post_full(data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                      rho = rho, calc = calc, others = others)
  return (-ll)
}

neg_log_post_alpha_full <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                                    others) {
  alpha$cur <- transform$inv.logit(q)
  ll <- log_post_full(data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                      rho = rho, calc = calc, others = others)
  return(-ll)
}

neg_log_post_beta <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                              others) {
  # extract from the list and get calculated quantities
  beta$cur <- q
  
  x.beta <- getXBeta(y = data$y, x = data$x, beta = beta$cur)
  z      <- getZ(xi = xi$cur, x.beta = x.beta, thresh = others$thresh)
  theta  <- getTheta(alpha = alpha$cur, z = z, aw = calc$aw)
  
  # log prior
  ll <- sum(dnorm(x = beta$cur, mean = beta$mn, sd = beta$sd, log = TRUE))
  ll <- ll + sum(logLikeY(y = data$y, theta = theta))
  
  return(-ll)
}

neg_log_post_a <- function(q, data, beta, xi, a, b, alpha, rho, calc, others) {
  # parameter transformation
  a$cur <- exp(q)
  
  # extract from the list and get calculated quantities
  alpha1m <- 1 - alpha$cur
  
  # log prior
  # Remember: q = log(a) and jacobian is included
  lc <- logc(b = b$cur, alpha = alpha$cur)
#   if (any(is.nan(lc))) {
#     b.trouble.lc <<- b
#     alpha.trouble.lc <<- alpha
#     stop("nan in logc")
#   }
  ll <- sum(-alpha$cur / alpha1m * q - exp(lc) * a$cur^(-alpha$cur / alpha1m))
  
  # data and log likelihood
  aw <- getAW(a = a$cur, w.star = calc$w.star)
#   if (any(is.nan(aw))) {
#     a.trouble.aw <<- a
#     calc.trouble.aw <<- calc
#     stop("nan in getAW")
#   }
  
  theta <- getTheta(alpha = alpha$cur, z = calc$z, aw = aw)
#   if (any(is.nan(theta))) {
#     alpha.trouble.theta <<- alpha
#     calc.trouble.theta  <<- calc
#     aw.trouble.theta    <<- aw
#     stop("nan in theta")
#   }
  
  ll <- ll + sum(logLikeY(y = data$y, theta = theta))
  
  return (-ll)
}

neg_log_post_alpha <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                               others) {
  # parameter transformation
  alpha$cur <- transform$inv.logit(q)
  w.star    <- getWStarIDs(alpha = alpha$cur, w = calc$w, IDs = others$IDs)
  aw        <- getAW(a = a$cur, w.star = w.star)
  theta     <- getTheta(alpha = alpha$cur, z = calc$z, aw = aw)
  alpha1m   <- 1 - alpha$cur
  
  # extract from the list and get calculated quantities
  lc <- logc(b = b$cur, alpha = alpha$cur)
  
  nt     <- ncol(a$cur)
  nknots <- nrow(a$cur)
  
  # log density for PS random variable
  ll <- sum(log(alpha$cur) - log(alpha1m) - 1 / alpha1m * log(a$cur) + lc - 
              exp(lc) * a$cur^(-alpha$cur / alpha1m))
  
  # data and log likelihood
  ll <- ll + sum(logLikeY(y = data$y, theta = theta))
  
  # prior
  ll <- ll + (alpha$a - 1) * log(alpha$cur) + (alpha$b - 1) * log(alpha1m) +
    log(alpha$cur) - log(alpha1m)  # jacobian
  
  return(-ll)
}

neg_log_post_a_alpha <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                                 others) {
  # parameter transformation
  nknots    <- nrow(a$cur)
  nt        <- ncol(a$cur)
  q.a       <- matrix(q[1:(nknots * nt)], nknots, nt)
  q.alpha   <- tail(q, 1)
  a$cur     <- exp(q.a)
  alpha$cur <- transform$inv.logit(q.alpha)
  alpha1m   <- 1 - alpha$cur
  
  # log prior
  # Remember: q = log(a) and jacobian is included
  lc <- logc(b = b$cur, alpha = alpha$cur)
  ll <- sum(-alpha$cur / alpha1m * q.a + lc - exp(lc) * a$cur^(-alpha$cur / alpha1m)) +
    nknots * nt * (log(alpha$cur) - log(alpha1m))
  
  ll <- ll + (alpha$a - 1) * log(alpha$cur) + (alpha$b - 1) * log(alpha1m) + 
    log(alpha$cur) - log(alpha1m)  # jacobian
  # data and log likelihood
  w.star <- getWStarIDs(alpha = alpha$cur, w = calc$w, IDs = others$IDs)
  aw     <- getAW(a = a$cur, w.star = w.star)
  theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = aw)
  ll <- ll + sum(logLikeY(y = data$y, theta = theta))
  
  return (-ll)
}

llps <- function(q, data, beta, xi, a, b, alpha, rho, calc, others, 
                 addup = TRUE) {
  # parameter transformation
  b$cur <- transform$inv.logit(q)
  
  # extract from the list and get calculated quantities
  alpha1m <- 1 - alpha$cur
  
  # log prior
  lc <- logc(b = b$cur, alpha = alpha$cur)
  ll <- lc - exp(lc) * a$cur^(-alpha$cur / alpha1m) + log(b$cur) + log(1 - b$cur)
  
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

neg_log_post_b <- function(q, data, beta, xi, a, b, alpha, rho, calc, others, 
                           addup = TRUE) {
  # parameter transformation
  b$cur <- transform$inv.logit(q)
  
  # extract from the list and get calculated quantities
  alpha1m <- 1 - alpha$cur
  
  # prior is U(0, 1)
  lc <- logc(b = b$cur, alpha = alpha$cur)
  ll <- lc - exp(lc) * a$cur^(-alpha$cur / alpha1m) + 
        log(b$cur) + log(1 - b$cur)  # jacobian
  
  # b does not have an impact on the likelihood
  
  if (addup) {
    ll <- sum(ll)
  }
  
  return (-ll)
}

neg_log_post_grad_a <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                                others, recalc = TRUE) {
  # parameter transformation
  a$cur  <- exp(q)
  
  # for the gradiant of a_alpha, we calculate this before calling this function
  if (recalc) {
    aw    <- getAW(a = a$cur, w.star = calc$w.star)
    theta <- getTheta(alpha = alpha$cur, z = calc$z, aw = aw)
  } else {
    aw    <- calc$aw
    theta <- calc$theta
  }
  
  nt     <- ncol(a$cur)
  nknots <- nrow(a$cur)
  
  # extract from the list and get calculated quantities
  lc      <- logc(b = b$cur, alpha = alpha$cur)
  alpha1m <- 1 - alpha$cur
  
  # grad wrt log(a)
  grad <- matrix(0, nknots, nt)
  
  # the likelihood component should be
  # For non-observances:
  #   a.lt * wz.star.lt
  # For observances:
  #   a.lt * wz.star.lt / expm1(theta.t)
  # lw.star <- log(calc$w.star)
  for (t in 1:nt) {
    these <- data$y[, t] == 1
    wz.star.t <- -exp((calc$lw - calc$lz[, t]) / alpha$cur)
    wz.star.t[these, ] <- -wz.star.t[these, ] / expm1(theta[these, t])
    theta.t <- theta[these, t]
    grad[, t] <- a$cur[, t] * colSums(wz.star.t)
  }
  
  checkpoint <- 1
  if (any(is.nan(grad))) {
    # print(a$cur[which(is.nan(grad)), ])
    # print(paste("nan at a checkpoint ", checkpoint, sep = ""))
  }
  
  grad <- grad - alpha$cur / alpha1m * (1 - exp(lc) * a$cur^(-alpha$cur / alpha1m))
  
#   if (any(is.nan(grad))) {
#     beta.trouble.a.grad <<- beta
#     xi.trouble.a.grad <<- xi
#     a.trouble.a.grad <<- a
#     b.trouble.a.grad <<- b
#     alpha.trouble.a.grad <<- alpha
#     rho.trouble.a.grad <<- rho
#     calc.trouble.a.grad <<- calc
#     aw.trouble.a.grad <<- aw
#     theta.trouble.a.grad <<- theta
#     lc.trouble.a.grad <<- lc
#   }
  
  checkpoint <- checkpoint + 1
  if (any(is.nan(grad))) {
    # print(a$cur[which(is.nan(grad)), ])
    # print(paste("nan at a checkpoint ", checkpoint, sep = ""))
  }
  
  # whole function is written as gradient of log likelihood
  return(-grad)
}

# neg_log_post_grad_a_s <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
#                                 others) {
#   # parameter transformation
#   a$cur  <- exp(q)
#   
#   # for the gradiant of a_alpha, we calculate this before calling this function
#   aw    <- getAW(a = a$cur, w.star = calc$w.star)
#   theta <- getTheta(alpha = alpha$cur, z = calc$z, aw = aw)
#   
#   nt     <- ncol(a$cur)
#   nknots <- nrow(a$cur)
#   
#   # extract from the list and get calculated quantities
#   lc      <- logc(b = b$cur, alpha = alpha$cur)
#   alpha1m <- 1 - alpha$cur
#   
#   # grad wrt log(a)
#   grad <- matrix(0, nknots, nt)
#   
#   # the likelihood component should be
#   # For non-observances:
#   #   a.lt * wz.star.lt
#   # For observances:
#   #   a.lt * wz.star.lt / expm1(theta.t)
#   lw.star <- log(calc$w.star)
#   for (t in 1:nt) {
#     wz.star.t <- exp(lw.star - log(calc$z[, t]) / alpha$cur)
#     these <- data$y[, t] == 1
#     theta.t <- theta[these, t]
#     for (l in 1:nknots) {
#     grad[l, t] <- a$cur[l, t] * 
#             (-sum(wz.star.t[!these, l]) + 
#                sum(wz.star.t[these, l] / expm1(theta.t)))
#     }
#     
#   }
#   grad <- grad - alpha$cur / alpha1m * (1 - exp(lc) * a$cur^(-alpha$cur / alpha1m))
#   
#   # whole function is written as gradient of log likelihood
#   return(-grad)
# }

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

neg_log_post_grad_b <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                                others) {
  # parameter transformation
  b$cur <- transform$inv.logit(q)
  
  # extract from the list and get calculated quantities
  alpha1m <- 1 - alpha$cur  # used a lot in gradiant for b
  apb     <- alpha$cur * pi * b$cur
  a1mpb   <- alpha1m * pi * b$cur
  pb      <- pi * b$cur
  c       <- exp(logc(b = b$cur, alpha = alpha$cur))
  
  # chain rule work
  u <- sin(apb)^(alpha$cur / alpha1m)
  v <- sin(pb)^(-1 / alpha1m)
  w <- sin(a1mpb)
  
  # chain rule derivatives
  dbdq <- b$cur * (1 - b$cur)  # Jacobian
  dudb <- alpha$cur^2 * pi / alpha1m * sin(apb)^(alpha$cur / alpha1m - 1) * cos(apb)
  dvdb <- - pi / alpha1m * sin(pb)^(-1 / alpha1m - 1) * cos(pb)
  dwdb <- cos(a1mpb) * alpha1m * pi
  dcdb <- u * v * dwdb + u * dvdb * w + dudb * v * w
  dcdq <- dcdb * dbdq
  
  # dcdb
  grad <- dcdq * (-1 / c + a$cur^(-alpha$cur / alpha1m)) + 2 * b$cur - 1
  
  return(grad)
}

neg_log_post_grad_alpha_s <- function(q, data, beta, xi, a, b, alpha, rho, calc,
                                    others, eps = 0.0001) {
  # q: logit(alpha)

  # using a quick and easy numerical approximation
  # might be good to eventually get the actual value, but this is pretty reliable
  # need to recalculate aw and theta when we add eps to q, otherwise it doesn't 
  #   take into account that there was some movement.

  grad <- (neg_log_post_alpha(q = q + eps, data = data, beta = beta, xi = xi, 
                              a = a, b = b, alpha = alpha, rho = rho, 
                              calc = calc, others = others) - 
           neg_log_post_alpha(q = q, data = data, beta = beta, xi = xi, a = a, 
                              b = b, alpha = alpha, rho = rho, calc = calc, 
                              others = others)) / eps
  
  return(grad)
}

neg_log_post_grad_alpha <- function(q, data, beta, xi, a, b, alpha, rho, calc,
                                    others, recalc = TRUE) {
  # q: logit(alpha)
  nknots <- nrow(a$cur)
  nt     <- ncol(data$y)
  ns     <- nrow(data$y)
  
  alpha$cur <- transform$inv.logit(q)
  if (recalc) {
    w.star <- getWStarIDs(alpha = alpha$cur, w = calc$w, IDs = others$IDs)
    aw     <- getAW(a = a$cur, w.star = w.star)
    theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = aw)
  } else {
    w.star <- calc$w.star
    aw     <- calc$aw
    theta  <- calc$theta
  }
  
  grad <- 0
  alpha.sq <- alpha$cur * alpha$cur
  la <- log(a$cur)
  alpha1m <- 1 - alpha$cur
  alpha1m.sq <- alpha1m * alpha1m
  apb <- alpha$cur * pi * b$cur
  pb  <- pi * b$cur
  a1mpb <- alpha1m * pi * b$cur
  spb  <- sin(pb)
  capb <- cos(apb)
  sapb <- sin(apb)
  # a.aa1m <- a$cur^(-alpha / alpha1m)
  a.aa1m <- exp(-alpha$cur * la / alpha1m)
  
  # likelihood
  # very concise way to express this. Likely as fast as we can get it
  # below, we set some of these to zero for numerical stability
  #   when w = 0, log(w / z) = -Inf, so diff * exp(diff) = 0 * -Inf = NaN
  for (t in 1:nt) {
    these <- data$y[, t] == 1
    diff.t <- calc$lw - calc$lz[, t]  # ns x nknots
    diff.t <- diff.t * exp(diff.t / alpha$cur)
    diff.t[these, ] <- -diff.t[these, ] / expm1(theta[these, t])
    diff.t[is.nan(diff.t)] <- 0  # gradient should be zero when weight is 0
    grad <- grad + sum(diff.t %*% a$cur)
  }
  
  checkpoint <- 1
  if (is.nan(grad)) {
    # print(alpha$cur)
    # print(paste("nan at alpha checkpoint ", checkpoint, sep = ""))
    return(grad)
  }
  
  grad <- grad / alpha.sq

  checkpoint <- checkpoint + 1
  if (is.nan(grad)) {
    # print(alpha$cur)
    # print(paste("nan at alpha checkpoint ", checkpoint, sep = ""))
    return(grad)
  }
  
#   part1 <- sum(1 / alpha + 1 / alpha1m - la / alpha1m.sq + 
#                  pb * capb / (alpha1m * sapb) +
#                  log(sapb) / alpha1m.sq - log(spb) / alpha1m.sq - 
#                  pb * cos(a1mpb) / sin(a1mpb) - pb * capb / sapb)
#   part2 <- pb * a.aa1m * capb * sin(-a1mpb) / 
#                  ((sapb/spb)^(-1 / alpha1m) * sapb^2) - 
#                  pb * a.aa1m * cos(-a1mpb) / ((sapb / spb)^(-1 / alpha1m) * sapb) - 
#                  a.aa1m * (-1 / alpha1m - alpha / alpha1m.sq) * la * sin(-a1mpb) / 
#                  ((sapb / spb)^(-1 / alpha1m) * sapb) + a.aa1m * 
#                  (pb * capb / (-alpha1m * sapb) - log(sapb / spb) / alpha1m.sq) * 
#                  sin(-a1mpb) / ((sapb / spb)^(-1 / alpha1m) * sapb)
#   
#   part2a <- (sapb/spb)^(-1 / alpha1m)
#   part2b <- log(sapb / spb)
#   
#   if (is.nan(part1)) {
#     print("part 1")
#   }
#   if (any(is.nan(part2))) {
#     print("part 2")
#     print(which(is.nan(part2)))
#     print(paste("pb", pb[which(is.nan(part2))]))
#     print(paste("a.aa1m", a.aa1m[which(is.nan(part2))]))
#     print(paste("capb", capb[which(is.nan(part2))]))
#     print(paste("a", a$cur[which(is.nan(part2))]))
#     print(paste("b", b$cur[which(is.nan(part2))]))
#   }
#   if (any(is.nan(part2a))) {
#     print("part 2a")
#   }
#   if (any(is.nan(part2b))) {
#     print("part 2b")
#   }
  
  grad <- grad + sum(1 / alpha$cur + 1 / alpha1m - la / alpha1m.sq + 
                       pb * capb / (alpha1m * sapb) +
                       log(sapb) / alpha1m.sq - log(spb) / alpha1m.sq - 
                       pb * cos(a1mpb) / sin(a1mpb) - pb * capb / sapb) - 
    sum(pb * a.aa1m * capb * sin(-a1mpb) / 
          ((sapb/spb)^(-1 / alpha1m) * sapb^2) - 
          pb * a.aa1m * cos(-a1mpb) / ((sapb / spb)^(-1 / alpha1m) * sapb) - 
          a.aa1m * (-1 / alpha1m - alpha$cur / alpha1m.sq) * la * sin(-a1mpb) / 
          ((sapb / spb)^(-1/alpha1m) * sapb) + a.aa1m * 
          (pb * capb / (-alpha1m * sapb) - log(sapb / spb) / alpha1m.sq) * 
          sin(-a1mpb) / ((sapb / spb)^(-1 / alpha1m) * sapb))
  
  checkpoint <- checkpoint + 1
  if (is.nan(grad)) {
    # print(alpha$cur)
    # print(paste("nan at alpha checkpoint ", checkpoint, sep = ""))
    return(grad)
  }

  # prior
  grad <- grad + (alpha$a - 1) / alpha$cur + (alpha$b - 1) / alpha1m

  grad <- grad * alpha$cur * alpha1m + 1 # Jacobian
  # grad <- grad - 1 / alpha$cur + 1 / alpha1m
  
  checkpoint <- checkpoint + 1
  if (is.nan(grad)) {
    # print(alpha$cur)
    # print(paste("nan at alpha checkpoint ", checkpoint, sep = ""))
    return(grad)
  }
  
  return(-grad)
}

neg_log_post_grad_beta <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                                   others, eps = 0.0001) {
  # q: beta
  
  # using a quick and easy numerical approximation
  # might be good to eventually get the actual value, but this is pretty reliable
  
  grad <- (neg_log_post_beta(q = q + eps, data = data, beta = beta, xi = xi, 
                             a = a, b = b, alpha = alpha, rho = rho, 
                             calc = calc, others = others) - 
           neg_log_post_beta(q = q, data = data, beta = beta, xi = xi, a = a,
                             b = b, alpha = alpha, rho = rho, calc = calc, 
                             others = others)) / eps
  
  return(grad)
}

neg_log_post_grad_a_alpha <- function(q, data, beta, xi, a, b, alpha, rho, calc, 
                                      others) {
  # parameter transformation
  nknots  <- ncol(calc$w)
  nt      <- ncol(data$y)
  q.a     <- matrix(q[1:nknots], nknots, nt)  # a needs to be a matrix
  q.alpha <- tail(q, 1)
  
  grad <- rep(NA, length(q))
  
  a$cur       <- exp(q.a)
  alpha$cur   <- transform$inv.logit(q.alpha)
  
  if (any(a$cur == Inf)) {
    print("Infinity in a")
    print(a$cur[which(a$cur == Inf)])
    print(q.a[which(a$cur == Inf)])
  }
  
  if (any (a$cur == 0)) {
    print("0 in a")
    print(which(a$cur == 0))
    print(a$cur[which(a$cur == 0)])
    print(q.a[which(a$cur == 0)])
  }
  
  # a and alpha only impact a.star, aw, and theta
  calc$w.star <- getWStarIDs(alpha = alpha$cur, w = calc$w, IDs = others$IDs)
  calc$aw     <- getAW(a = a$cur, w.star = calc$w.star)
  calc$theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
  
  grad[1:nknots] <- neg_log_post_grad_a(q = q.a, data = data, beta = beta, 
                                        xi = xi, a = a, b = b, alpha = alpha, 
                                        rho = rho, calc = calc, others = others,
                                        recalc = FALSE)
  grad[(nknots + 1)] <- neg_log_post_grad_alpha(q = q.alpha, data = data, 
                                                beta = beta, xi = xi, a = a, 
                                                b = b, alpha = alpha, rho = rho, 
                                                calc = calc, others = others,
                                                recalc = FALSE)
  
  return(grad)
}

# neg_log_post_grad_a_alpha <- function(q, data, beta, xi, a, b, alpha, rho, calc,
#                                       others) {
#   
#   # parameter transformation
#   nknots  <- ncol(calc$w)
#   nt      <- ncol(data$y)
#   nkt     <- nknots * nt
#   q.a     <- matrix(q[1:nknots], nknots, nt)  # a needs to be a matrix
#   q.alpha <- tail(q, 1)
#   
#   grad <- rep(NA, length(q))
#   
#   a$cur       <- exp(q.a)
#   alpha$cur   <- transform$inv.logit(q.alpha)
#   
#   # for the gradiant of a_alpha, we calculate this before calling this function
#   w.star <- getWStar(alpha = alpha$cur, w = calc$w)
#   aw    <- getAW(a = a$cur, w.star = w.star)
#   theta <- getTheta(alpha = alpha$cur, z = calc$z, aw = aw)
# 
#   # extract from the list and get calculated quantities
#   lc      <- logc(b = b$cur, alpha = alpha$cur)
#   alpha1m <- 1 - alpha$cur
#   
#   alpha.sq <- alpha * alpha
#   la <- log(a$cur)
#   alpha1m <- 1 - alpha
#   alpha1m.sq <- alpha1m * alpha1m
#   apb <- alpha * pi * b$cur
#   pb  <- pi * b$cur
#   a1mpb <- alpha1m * pi * b$cur
#   spb  <- sin(pb)
#   capb <- cos(apb)
#   sapb <- sin(apb)
#   a.aa1m <- a$cur^(-alpha / alpha1m)
#   
#   # the likelihood component should be
#   # For non-observances:
#   #   a.lt * wz.star.lt
#   # For observances:
#   #   a.lt * wz.star.lt / expm1(theta.t)
#   # lw.star <- log(calc$w.star)
# 
#   # likelihood
#   # very concise way to express this. Likely as fast as we can get it
#   # below, we set some of these to zero for numerical stability
#   #   when w = 0, log(w / z) = -Inf, so diff * exp(diff) = 0 * -Inf = NaN
#   for (t in 1:nt) {
#     these <- data$y[, t] == 1
#     diff.t <- calc$lw - calc$lz[, t]  # ns x nknots
#     diff.t <- diff.t * exp(diff.t / alpha)
#     diff.t[these, ] <- -diff.t[these, ] / expm1(theta[these, t])
#     diff.t[is.nan(diff.t)] <- 0  # gradient should be zero when weight is 0
#     grad <- grad + sum(diff.t %*% a$cur)
#   }
#   
#   for (t in 1:nt) {
#     these <- data$y[, t] == 1
#     wz.star.t <- -exp((calc$lw - calc$lz[, t]) / alpha$cur)
#     wz.star.t[these, ] <- -wz.star.t[these, ] / expm1(theta[these, t])
#     theta.t <- theta[these, t]
#     grad[, t] <- a$cur[, t] * colSums(wz.star.t)
#   }
#   
#   grad <- grad - alpha$cur / alpha1m * (1 - exp(lc) * a$cur^(-alpha$cur / alpha1m))
#   
#   if (any(is.nan(grad))) {
#     beta.trouble.a.grad <<- beta
#     xi.trouble.a.grad <<- xi
#     a.trouble.a.grad <<- a
#     b.trouble.a.grad <<- b
#     alpha.trouble.a.grad <<- alpha
#     rho.trouble.a.grad <<- rho
#     calc.trouble.a.grad <<- calc
#     aw.trouble.a.grad <<- aw
#     theta.trouble.a.grad <<- theta
#     lc.trouble.a.grad <<- lc
#   }
#   
#   grad <- grad / alpha.sq
#   
#   grad <- grad + sum(1 / alpha + 1 / alpha1m - la / alpha1m.sq + 
#                        pb * capb / (alpha1m * sapb) +
#                        log(sapb) / alpha1m.sq - log(spb) / alpha1m.sq - 
#                        pb * cos(a1mpb) / sin(a1mpb) - pb * capb / sapb) - 
#     sum(pb * a.aa1m * capb * sin(-a1mpb) / 
#           ((sapb/spb)^(-1 / alpha1m) * sapb^2) - 
#           pb * a.aa1m * cos(-a1mpb) / ((sapb / spb)^(-1 / alpha1m) * sapb) - 
#           a.aa1m * (-1 / alpha1m - alpha / alpha1m.sq) * la * sin(-a1mpb) / 
#           ((sapb / spb)^(-1/alpha1m) * sapb) + a.aa1m * 
#           (pb * capb / (-alpha1m * sapb) - log(sapb / spb) / alpha1m.sq) * 
#           sin(-a1mpb) / ((sapb / spb)^(-1 / alpha1m) * sapb))
#   
#   grad <- grad * alpha * alpha1m  # Jacobian
#   
#   # whole function is written as gradient of log likelihood
#   return(-grad)
# }
