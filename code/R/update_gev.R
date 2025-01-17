updateBeta <- function(data, beta, xi, alpha, calc, others) {
  np  <- length(beta$cur)
  ns  <- nrow(data$y)
  nt  <- ncol(data$y)
  
  cur.lly <- logLikeY(y = data$y, theta = calc$theta)
  cur.beta <- beta$cur
  
  can.beta <- rnorm(np, beta$cur, beta$eps)
  x.beta   <- getXBeta(y = data$y, x = data$x, beta = can.beta) 
  
  if (any(xi$cur * (x.beta - others$thresh) > 1)) {  # numerical stability
    can.lly <- -Inf
  } else {
    z       <- getZ(xi = xi$cur, x.beta = x.beta, thresh = others$thres)
    theta   <- getTheta(alpha = alpha$cur, z = z, aw = calc$aw)
    can.lly <- logLikeY(y = data$y, theta = theta)
  }
  
  R <- sum(can.lly - cur.lly) +
    sum(dnorm(can.beta, beta$mn, beta$sd, log = TRUE) -
          dnorm(cur.beta, beta$mn, beta$sd, log = TRUE))
  
  if (!is.na(R)) { if (log(runif(1)) < R) {
    results <- list(q = can.beta, accept = TRUE)
  } else {
    results <- list(q = cur.beta, accept = FALSE)
  }} else {
    results <- list(q = cur.beta, accept = FALSE)
  }
  
  return(results)
}

updateXi <- function(data, xi, alpha, calc, others) {
  nt  <- ncol(data$y)
  
  cur.lly <- logLikeY(y = data$y, theta = calc$theta)
  cur.xi <- xi$cur
  
  can.xi <- rnorm(1, xi$cur, xi$eps)
  
  if (any(can.xi * (calc$x.beta - others$thresh) > 1)) {
    can.lly <- -Inf
  } else {
    z       <- getZ(xi = can.xi, x.beta = calc$x.beta, thresh = others$thresh)
    theta   <- getTheta(alpha = alpha$cur, z = z, aw = calc$aw)
    can.lly <- logLikeY(y = data$y, theta = theta)
  }
  
  R <- sum(can.lly - cur.lly) +
    dnorm(can.xi, xi$mn, xi$sd, log = TRUE) - 
    dnorm(cur.xi, xi$mn, xi$sd, log = TRUE)
  
  if (!is.na(R)) { if (log(runif(1)) < R) {
    results <- list(q = can.xi, accept = TRUE)
  } else {
    results <- list(q = cur.xi, accept = FALSE)
  }} else {
    results <- list(q = cur.xi, accept = FALSE)
  }
  
  return(results)
}

updateRho <- function(data, a, alpha, rho, calc, others) {
  nt     <- ncol(data$y)
  nknots <- nrow(a$cur)
  
  cur.lly <- logLikeY(y = data$y, theta = calc$theta)
  cur.rho <- rho$cur
  
  # can.rho <- exp(rnorm(1, log(rho$cur), rho$eps))
  cur.rho.star <- transform$logit(cur.rho, rho$lower, rho$upper)
  can.rho.star <- rnorm(1, cur.rho.star, rho$eps)
  can.rho      <- transform$inv.logit(can.rho.star, rho$lower, rho$upper)
  
  if (can.rho < max(sqrt(others$dw2)) & can.rho > 1e-6) {
    w       <- getW(rho = can.rho, dw2 = others$dw2, a.cutoff = others$a.cutoff)
    w.star  <- getWStarIDs(alpha = alpha$cur, w = w, IDs = others$IDs)
    aw      <- getAW(a = a$cur, w.star = w.star)
    theta   <- getTheta(alpha = alpha$cur, z = calc$z, aw = aw)
    can.lly <- logLikeY(y = data$y, theta = theta)
    
    R <- sum(can.lly - cur.lly) +
      log(can.rho - rho$lower) + log(rho$upper - can.rho) - # Jacobian
      log(cur.rho - rho$lower) - log(rho$upper - cur.rho)
      # dnorm(log(can.rho), rho$mn, rho$sd, log=TRUE) -
      # dnorm(log(cur.rho), rho$mn, rho$sd, log=TRUE)
  } else {
    R <- -Inf
  }
  
  if (!is.na(R)) { if (log(runif(1)) < R) {
    results <- list(q = can.rho, accept = TRUE)
  } else {
    results <- list(q = cur.rho, accept = FALSE)
  }} else {
    results <- list(q = cur.rho, accept = FALSE)
  }

  return(results)
}

# update the alpha term for theta.star
updateAlpha <- function(data, a, b, alpha, calc, others) {
  nt     <- ncol(data$y)
  nknots <- nrow(a$cur)
  
  cur.alpha   <- alpha$cur
  cur.alpha1m <- 1 - alpha$cur
  lc          <- logc(b = b$cur, alpha = cur.alpha)
  
  cur.lly <- logLikeY(y = data$y, theta = calc$theta)
  cur.llps <- sum(log(cur.alpha) - log(cur.alpha1m) - 
                    1 / cur.alpha1m * log(a$cur) + lc - 
                    exp(lc) * a$cur^(-cur.alpha / cur.alpha1m))
  
  alpha.star     <- transform$logit(cur.alpha)
  can.alpha.star <- rnorm(1, alpha.star, alpha$eps)
  can.alpha      <- transform$inv.logit(can.alpha.star)
  # can.alpha    <- pnorm(rnorm(1, qnorm(cur.alpha), alpha$eps))
  can.alpha1m    <- 1 - can.alpha
  can.w.star     <- getWStarIDs(alpha = can.alpha, w = calc$w, IDs =others$IDs)
  can.aw         <- getAW(a = a$cur, w.star = can.w.star)
  can.theta      <- getTheta(alpha = can.alpha, z = calc$z, aw = can.aw)
  lc             <- logc(b = b$cur, alpha = can.alpha)
  
  can.lly <- logLikeY(y = data$y, theta = can.theta)
  can.llps <- sum(log(can.alpha) - log(can.alpha1m) - 
                    1 / can.alpha1m * log(a$cur) + lc - 
                    exp(lc) * a$cur^(-can.alpha / can.alpha1m))
  
  R <- sum(can.lly - cur.lly) + sum(can.llps - cur.llps) +
    dbeta(can.alpha, alpha$a, alpha$b, log = TRUE) -
    dbeta(cur.alpha, alpha$a, alpha$b, log = TRUE) + 
    log(can.alpha - 0) + log(1 - can.alpha) - # Jacobian of the prior
    log(cur.alpha - 0) - log(1 - cur.alpha)
  
  if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
    results <- list(q = can.alpha, accept = TRUE)
  } else {
    results <- list(q = cur.alpha, accept = FALSE)
  }} else {
    results <- list(q = cur.alpha, accept = FALSE)
  }
  
  return(results)
}

# rewrite to use theta function not theta.star
pred.spgev <- function(mcmcoutput, s.pred, x.pred, knots, start = 1, end = NULL,
                       thin = 1, thresh = 0, update = NULL) {
  if (is.null(end)) {
    end <- length(mcmcoutput$xi)
  }
  
  np     <- nrow(s.pred)
  niters <- length(start:end)
  nknots <- dim(mcmcoutput$a)[2]
  
  if (is.null(dim(mcmcoutput$beta))) {
    p <- 1
  } else {
    p <- dim(mcmcoutput$beta)[2]
  }
  
  if (is.null(dim(x.pred))) {
    nt <- 1
  } else {
    if (p == 1) {
      nt <- ncol(x.pred)
    } else {
      nt <- dim(x.pred)[2]
    }
  }
  
  niters <- length(start:end)
  # b is only an auxiliary variable for model fitting
  beta  <- matrix(mcmcoutput$beta[start:end, , drop = F], niters, p)
  xi    <- mcmcoutput$xi[start:end]
  a     <- mcmcoutput$a[start:end, , , drop = F]
  alpha <- mcmcoutput$alpha[start:end]
  rho   <- mcmcoutput$rho[start:end]
  
  dw2p  <- as.matrix(rdist(s.pred, knots))^2
  a.cutoff <- mcmcoutput$a.cutoff
  IDs.p <- getIDs(dw2 = dw2p, a.cutoff = a.cutoff)
  
  # y.pred <- matrix(NA, nrow=(niters / thin), ncol = np)
  prob.success <- matrix(NA, nrow=(niters / thin), ncol=np)
  x.beta <- matrix(NA, np, nt)
  
  for (i in 1:length(start:end)) {
    if (i %% thin == 0) {
      alpha.i <- alpha[i]
      x.beta  <- getXBeta(x = x.pred, beta = beta[i, ], ns = np, nt = nt)
      z       <- getZ(xi = xi[i], x.beta = x.beta, thresh = thresh)
      w       <- getW(rho = rho[i], dw2 = dw2p, a.cutoff = a.cutoff)
      # there are certain circumstances when a place where we want to predict
      # is too far away from all the knot locations, so the kernel functions 
      # are 0 at all sites. when this happens w = 0 / sum(0) = nan, so setting 
      # w = 0 when we get NaN from getW.
      w[is.nan(w)] <- 0
      
      w.star  <- getWStarIDs(alpha = alpha.i, w = w, IDs = IDs.p)
      aw      <- getAW(a = a[i, , ], w.star = w.star)
      theta   <- getTheta(alpha = alpha.i, z = z, aw = aw)
      prob.success[(i / thin), ] <- 1 - exp(-theta)
      
      # if z is nan, it means that x.beta is such a large number there is
      # basically 0 probability that z < 0
      if (any(is.nan(z))) {
        these <- which(is.nan(z))
        prob.success[(i / thin), these] <- 1
      }
      
      # y.pred[(i / thin), ] <- rbinom(n = np, size = 1, prob = prob.success)
      
    }
    if (!is.null(update)) {
      if (i %% update == 0) {
        cat("\t Iter", i, "\n")
      }
    }
  }
  
  return(prob.success)
  
}

# updateBeta <- function(y, theta, alpha, a.star, z, w, wz, beta, beta.m, beta.s,
#                        x.beta, xi, x, cur.lly, acc, att, mh, thresh = 0) {
#   np  <- length(beta)
#   ns  <- nrow(y)
#   nt  <- ncol(y)
# 
#   for (p in 1:np) {
#     att[p]      <- att[p] + 1
#     can.beta    <- rnorm(1, beta[p], mh[p])
#     # trying to save a little time by just updating for the pth column of x
#     if (nt == 1) {
#       can.x.beta <- x.beta + x[, p] * (can.beta - beta[p])
#     } else {
#       can.x.beta <- matrix(NA, ns, nt)
#       for (t in 1:nt) {
#         start <- (t - 1) * ns + 1
#         end   <- t * ns
#         can.x.beta[, t]  <- x.beta[, t] + x[start:end, p] * (can.beta - beta[p])
#       }
#     }
# 
#     if (any(xi * (can.x.beta - thresh) > 1)) {  # numerical stability
#       can.lly <- -Inf
#     } else {
#       can.z     <- getZ(xi = xi, x.beta = can.x.beta, thresh = thresh)
#       can.wz    <- getwzCPP(z = can.z, w = w)
#       can.theta <- getThetaCPP(wz = can.wz, a_star = a.star, alpha = alpha)
#       can.lly   <- logLikeY(y = y, theta = can.theta)
#     }
# 
#     R <- sum(can.lly - cur.lly) +
#          dnorm(can.beta, beta.m, beta.s, log = TRUE) -
#          dnorm(beta[p], beta.m, beta.s, log = TRUE)
# 
#     if (!is.na(R)) { if (log(runif(1)) < R) {
#       beta[p] <- can.beta
#       x.beta  <- can.x.beta
#       z       <- can.z
#       wz      <- can.wz
#       theta   <- can.theta
#       cur.lly <- can.lly
#       acc[p]  <- acc[p] + 1
#     }}
#   }
# 
#   results <- list(beta = beta, x.beta = x.beta, z = z, wz = wz,
#                   theta = theta, cur.lly = cur.lly, att = att, acc = acc)
#   return(results)
# }

# updateXi <- function(y, theta, alpha, a.star, z, w, wz, x.beta, xi, xi.m, xi.s,
#                      cur.lly, acc, att, mh, thresh = 0) {
#   nt  <- ncol(y)
# 
#   att <- att + 1
#   can.xi     <- rnorm(1, xi, mh)
# 
#   if (sum(can.xi * (x.beta - thresh) > 1) > 0) {
#     can.lly <- -Inf
#   } else {
#     can.z     <- getZ(xi = can.xi, x.beta = x.beta, thresh = thresh)
#     can.wz    <- getwzStarCPP(z = can.z, w = w)
#     can.theta <- getThetaCPP(wz = can.wz, a_star = a.star, alpha = alpha)
#     can.lly   <- logLikeY(y = y, theta = can.theta)
#   }
# 
#   R <- sum(can.lly - cur.lly) +
#        dnorm(can.xi, xi.m, xi.s, log = TRUE) - dnorm(xi, xi.m, xi.s, log = TRUE)
# 
#   if (!is.na(R)) { if (log(runif(1)) < R) {
#     xi      <- can.xi
#     z       <- can.z
#     wz      <- can.wz
#     theta   <- can.theta
#     cur.lly <- can.lly
#     acc     <- acc + 1
#   }}
# 
#   results <- list(xi = xi, z = z, wz = wz, theta = theta,
#                   cur.lly = cur.lly, att = att, acc = acc)
#   return(results)
# }

# updateRho <- function(y, theta, a.star, alpha, cur.lly, w, z, wz, dw2,
#                       rho, logrho.m, logrho.s, rho.upper = Inf, A.cutoff,
#                       acc, att, mh) {
#   nt     <- ncol(y)
#   nknots <- nrow(a.star)
#   
#   att <- att + 1
#   
#   rho.star     <- log(rho)
#   can.rho.star <- rnorm(1, rho.star, mh)
#   can.rho      <- exp(can.rho.star)
#   can.w        <- stdW(makeW(dw2 = dw2, rho = can.rho, A.cutoff = A.cutoff))
#   can.wz       <- getwzCPP(z = z, w = can.w)
#   can.theta    <- getThetaCPP(wz = can.wz, a_star = a.star, alpha = alpha)
#   can.lly      <- logLikeY(y = y, theta = can.theta)
#   
#   R <- sum(can.lly - cur.lly) +
#     dnorm(can.rho.star, logrho.m, logrho.s, log=TRUE) -
#     dnorm(rho.star, logrho.m, logrho.s, log=TRUE)
#   
#   if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
#     rho     <- can.rho
#     w       <- can.w
#     wz      <- can.wz
#     theta   <- can.theta
#     cur.lly <- can.lly
#     acc     <- acc + 1
#   }}
#   
#   results <- list(rho = rho, w = w, wz = wz, theta = theta,
#                   cur.lly = cur.lly, att = att, acc = acc)
#   return(results)
# }

updateXiBeta <- function(y, alpha, z, w, wz, beta, theta, a.star,
                         x.beta, xi, x, xt, xtx.inv, xi.m, xi.s, cur.lly,
                         xi.fix, beta.fix, acc.p, att.p, mh.p, acc.xi, att.xi,
                         mh.xi, thresh = 0) {
  ns <- nrow(y)
  nt <- ncol(y)
  np <- length(beta)
  alpha.inv <- 1 / alpha

  att.p  <- att.p + 1
  if (!xi.fix) {
    att.xi <- att.xi + 1
  }

  # trying to get xi and beta using the marginal P(Y = 0)
  cur.p     <- exp(-1 / z)  # will be ns x nt matrix
  if (np == 1) {
    # intercept only model all elements in cur.p should be the same
    can.pstar <- rnorm(1, log(cur.p[1]/(1 - cur.p[1])), mh.p)
    can.p     <- 1 / (1 + exp(-can.pstar))
    can.p     <- matrix(can.p, ns, nt)
  } else {  # covariates in the model
    can.pstar <- matrix(rnorm(ns * nt, log(cur.p / (1 - cur.p)), mh.p), ns, nt)
    can.p     <- 1 / (1 + exp(-can.pstar))
  }

  # extract beta for the new xi term
  if (xi.fix) {
    can.xi <- xi
  } else {
    can.xi <- rnorm(1, xi, mh.xi)
  }

  if (xi < 1e-6 & xi > -1e-6) {
    can.x.beta <- log(-log(can.p))

  } else {
    can.x.beta <- (1 - (-1 / log(can.p))^can.xi) / can.xi
  }

  if (any(can.xi * (can.x.beta - thresh) > 1)) {
    can.lly <- -Inf
  } else {
    can.z     <- getZ(xi = can.xi, x.beta = can.x.beta, thresh = thresh)
    can.wz    <- getwzCPP(z = can.z, w = w)
    can.theta <- getThetaCPP(wz = can.wz, a_star = a.star, alpha = alpha)
    can.lly   <- logLikeY(y = y, theta = can.theta)
  }

  R <- sum(can.lly - cur.lly) +
       dbeta(can.p[1], 1, 1, log = TRUE) - dbeta(cur.p[1], 1, 1, log = TRUE)

  if (!xi.fix) {
    R <- R + dnorm(can.xi, xi.m, xi.s, log = TRUE) -
             dnorm(xi, xi.m, xi.s, log = TRUE)
  }

  if (!is.na(R)) { if (log(runif(1)) < R) {
    beta    <- xtx.inv %*% (xt %*% can.x.beta)
    x.beta  <- can.x.beta
    z       <- can.z
    wz      <- can.wz
    theta   <- can.theta
    cur.lly <- can.lly
    if (!xi.fix) {
      xi      <- can.xi
      acc.xi  <- acc.xi + 1
    }
    acc.p   <- acc.p + 1
  }}

  results <- list(beta = beta, x.beta = x.beta, xi = xi, z = z,
                  wz = wz, theta = theta, cur.lly = cur.lly,
                  att.p = att.p, acc.p = acc.p,
                  att.xi = att.xi, acc.xi = acc.xi)
}

# update the random effects for theta
updateA <- function(y, theta, a, a.star, alpha, wz, cur.lly, cur.llps,
                    mid.points, bin.width, mh, cuts, IDs, threads = 1) {
  nt     <- ncol(y)
  nknots <- nrow(a)
  
  alpha.inv <- 1 / alpha

  for (t in 1:nt) {
    cur.lly.t <- cur.lly[, t]
    for (k in 1:nknots) {
      cur.a <- a[k, t]
      l1    <- get.level(cur.a, cuts)  # keeps sd reasonable
      can.a <- exp(rnorm(1, log(cur.a), mh[l1]))
      can.a.star <- can.a^alpha
      l2    <- get.level(can.a, cuts)

      # can.theta only changes at a site when it's near the knot
      # just a vector for the day's theta values
      these <- IDs[[k]]  # get sites that are impacted by knot location
      wz.these <- wz[these, k, t]
      can.theta <- theta[these, t] + 
                   (wz.these * can.a.star)^alpha.inv - 
                   (wz.these * cur.a^alpha)^alpha.inv
      
      # Note: There are potential numerical stability issues here. If alpha
      # is really small, then we could end up with wz.star = 0 depending on 
      # how small w was to begin with. This really only becomes an issue when
      # Inf is drawn for the candidate value as 0 * Inf = NaN.
      # This is why we do the multiplication as 
      # (A_l^alpha * w_li / z_li)^(1 / alpha)
      
      if (can.a == Inf) { # numerical stability checks
        can.llps  <- -Inf
        can.lly.t <- -Inf
      } else if (any(can.theta <= 0)) {  
        can.llps  <- -Inf
        can.lly.t <- -Inf
      } else {
        can.llps  <- dPS.Rcpp(a = can.a, alpha = alpha,
                              mid.points = mid.points, bin.width = bin.width,
                              threads = threads)
        can.lly.t <- logLikeY(y = y[these, t, drop = F], theat = can.theta)
      }

      R <- sum(can.lly.t - cur.lly.t[these]) +
           can.llps - cur.llps[k, t] +
           dlognormal(cur.a, can.a, mh[l2]) - # candidate sd changes
           dlognormal(can.a, cur.a, mh[l1])

      if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
        a[k, t]          <- can.a
        a.star[k, t]     <- can.a.star
        theta[these, t]  <- can.theta
        cur.lly.t[these] <- can.lly.t
        cur.llps[k, t]   <- can.llps
      }}
    }
  }

  cur.lly <- logLikeY(y = y, theta = theta)

  results <- list(a = a, a.star = a.star, theta = theta,
                  cur.lly = cur.lly, cur.llps = cur.llps)
  return(results)
}

# # rewrite to use theta function not theta.star
# pred.spgev <- function(mcmcoutput, s.pred, x.pred, knots, start = 1, end = NULL,
#                        thin = 1, thresh = 0, update = NULL) {
#   if (is.null(end)) {
#     end <- length(mcmcoutput$xi)
#   }
# 
#   np     <- nrow(s.pred)
#   niters <- length(start:end)
#   nknots <- dim(mcmcoutput$a)[2]
# 
#   if (is.null(dim(mcmcoutput$beta))) {
#     p <- 1
#   } else {
#     p <- dim(mcmcoutput$beta)[2]
#   }
# 
#   if (is.null(dim(x.pred))) {
#     nt <- 1
#   } else {
#     if (p == 1) {
#       nt <- ncol(x.pred)
#     } else {
#       nt <- dim(x.pred)[2]
#     }
#   }
# 
#   niters <- length(start:end)
#   beta  <- matrix(mcmcoutput$beta[start:end, , drop = F], niters, p)
#   xi    <- mcmcoutput$xi[start:end]
#   a     <- mcmcoutput$a[start:end, , , drop = F]
#   alpha <- mcmcoutput$alpha[start:end]
#   rho   <- mcmcoutput$rho[start:end]
# 
#   dw2p  <- as.matrix(rdist(s.pred, knots))^2
#   a.cutoff <- mcmcoutput$a.cutoff
# 
#   prob.success <- matrix(NA, nrow=niters, ncol=np)
#   x.beta <- matrix(NA, np, nt)
#   
#   for (i in 1:length(start:end)) {
#     alpha.i <- alpha[i]
#     x.beta  <- getXBeta(x.pred, ns = np, nt = nt, beta = beta[i, ])
#     z       <- getZ(xi = xi[i], x.beta=x.beta, thresh=thresh)
#     w       <- stdW(makeW(dw2 = dw2p, rho = rho[i], a.cutoff = a.cutoff))
#     wz      <- getwzCPP(z = z, w = w)
#     theta   <- getThetaCPP(wz = wz, 
#                             a_star = matrix(a[i, , ], nknots, nt)^alpha.i, 
#                             alpha = alpha.i)
#     prob.success[i, ] <- 1 - exp(-theta)
# 
#     # if z is nan, it means that x.beta is such a large number there is
#     # basically 0 probability that z < 0
#     if (any(is.nan(z))) {
#       these <- which(is.nan(z))
#       prob.success[i, these] <- 1
#     }
# 
#     if (!is.null(update)) {
#       if (i %% update == 0) {
#         cat("\t Iter", i, "\n")
#       }
#     }
#   }
# 
#   return(prob.success)
# 
# }

predictY <- function(mcmcoutput, s.pred, x.pred, knots, start = 1, end = NULL,
                     thin = 1, thresh = 0, update = NULL) {

  if (is.null(end)) {
    end <- nrow(mcmcoutput$beta)
  }

  np   <- nrow(s.pred)
  iters <- end - start + 1
  dw2p <- as.matrix(rdist(s.pred, knots))^2
  yp <- matrix(NA, nrow = iters, ncol = np)
  for (i in start:end) {
    yp[i, ] <- rRareBinarySpat(x = x.pred, s = s.pred, knots = knots,
                               beta = mcmcoutput$beta[i, ],
                               xi = mcmcoutput$xi[i],
                               alpha = mcmcoutput$alpha[i],
                               rho = mcmcoutput$rho[i],
                               thresh = thresh, dw2 = dw2p,
                               a = mcmcoutput$a[i, ])$y
    if (sum(is.nan(yp[i, ])) > 0) {
      cat("beta", mcmcoutput$beta[i, ], "\n")
      cat("xi", mcmcoutput$xi[i], "\n")
      cat("alpha", mcmcoutput$alpha[i], "\n")
      cat("rho", mcmcoutput$rho[i], "\n")
      stop()
    }
    if (!is.null(update)) {
      if (i %% update == 0) {
        cat("\t Iter", i, "\n")
      }
    }
  }

  return(yp)

}

# # We might try to do a block update for everything
# # Yeah, that's what I was thinking.  I realize this would be more coding, 
# # but it runs 20 times faster and converges better it would be worth it.  
# # And maybe you could update mu, xi, and A in a batch?
# updateXiBetaA <- function(y, alpha, z, w, wz, beta, kernel, a, a.star,
#                           x.beta, xi, x, xi.m, xi.s,xi.fix, beta.fix, 
#                           acc.beta, att.beta, mh.beta, 
#                           acc.xi, att.xi, mh.xi, 
#                           cur.lly, cur.llps, thresh = 0,
#                           mh.a, cuts, IDs, threads = 1) {
#   
#   np  <- length(beta)
#   ns  <- nrow(y)
#   nt  <- ncol(y)
#   
#   # get candidate for beta
#   att.beta   <- att.beta + 1
#   can.beta   <- rnorm(np, beta, mh.beta)
#   can.x.beta <- getXBeta(x = x, ns = ns, nt = nt, beta = can.beta)
#     
#   # get candidate for xi
#   if (!xi.fix) {
#     att.xi <- att.xi + 1
#     can.xi <- rnorm(1, xi, mh.xi)
#   } else {
#     can.xi <- xi
#   }
#   
#   # get candidate for A terms
#   # TODO: Make adjustment to candidate for Langevin updates
#   can.a <- matrix(NA, nknots, nt)
#   for (t in 1: nt) {
#     for (k in 1:nknots) {
#       cur.a       <- a[k, t]
#       l1          <- get.level(cur.a, cuts)  # keeps sd reasonable
#       can.a.kt    <- exp(rnorm(1, log(cur.a), mh[l1]))
#       l2          <- get.level(can.a.kt, cuts)
#       can.a[k, t] <- can.a.kt
#     }
#   }
#     
#   can.a.star <- exp(alpha * can.a)
#   can.kernel <- getKernelCPP(wz = can.wz, a_star = a.star, alpha = alpha)
#   
#   if (any(can.xi * (can.x.beta - thresh) > 1)) {  # numerical stability
#     can.lly  <- -Inf
#     can.llps <- -Inf
#   } else if (any(can.a) == Inf) {
#     can.lly  <- -Inf
#     can.llps <- -Inf
#   } else if (any(can.kernel <= 0)) {
#     can.lly  <- -Inf
#     can.llps <- -Inf
#   } else {
#     can.z      <- getZ(xi = can.xi, x.beta = can.x.beta, thresh = thresh)
#     can.wz     <- getwzCPP(z = can.z, w = w)
#     can.kernel <- getKernelCPP(wz = can.wz, a_star = can.a.star, alpha = alpha)
#     can.lly    <- logLikeY(y = y, kernel = can.kernel)
#     can.llps   <- ld(u = u, a = can.a, alpha = alpha)
#   }
#   
#   R <- can.lly - cur.lly + can.llps - cur.llps + 
#        dnorm(can.beta, beta.m, beta.s, log = TRUE) -
#        dnorm(beta, beta.m, beta.s, log = TRUE)
#        # TODO: Add in langevin update here 
#        
#   
#   if (!xi.fix) {
#     R <- R + dnorm(can.xi, xi.m, xi.s, log = TRUE) - 
#              dnorm(xi, xi.m, xi.s, log = TRUE)
#   }
#   
#   if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
#     beta     <- can.beta
#     x.beta   <- can.x.beta
#     acc.beta <- acc.beta + 1
#     if (!xi.fix) {
#       xi     <- can.xi
#       acc.xi <- acc.xi + 1
#     }
#     z        <- can.z
#     wz       <- can.wz
#     a        <- can.a
#     a.star   <- can.a.star
#     kernel   <- can.kernel
#     cur.lly  <- can.lly
#     cur.llps <- can.llps
#   }}
#        
#   results <- list(beta = beta, x.beta = x.beta, xi = xi, z = z,
#                   wz = wz, kernel = kernel,
#                   att.beta = att.beta, acc.beta = acc.beta,
#                   att.xi = att.xi, acc.xi = acc.xi,
#                   a = a, a.star = a.star,
#                   cur.lly = cur.lly, cur.llps = cur.llps)
#   return(results)
# }