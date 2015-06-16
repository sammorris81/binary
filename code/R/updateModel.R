# update the beta term for the linear model
# TODO: still need to update this for using the kernel of logLikeY2
updateBeta <- function(y, theta.star, alpha, z, z.star, beta, beta.m, beta.s,
                       x.beta, xi, x, cur.lly, acc, att, mh, thresh=0) {
  # tried a block update for the beta update, but it doesn't do as well
  # as individual updates for each beta term separately.
  np  <- length(beta)
  ns  <- nrow(y)
  nt  <- ncol(y)
  alpha.inv <- 1 / alpha

  for (p in 1:np) {
    att[p]      <- att[p] + 1
    can.beta    <- rnorm(1, beta[p], mh[p])
    # trying to save a little time
    if (nt == 1) {
      can.x.beta <- x.beta + x[, p] * (can.beta - beta[p])
    } else {
      can.x.beta <- matrix(NA, ns, nt)
      for (t in 1:nt) {
        start <- (t - 1) * ns + 1
        end   <- t * ns
        can.x.beta[, t]  <- x.beta[, t] + x[start:end, p] * (can.beta - beta[p])
      }
    }

    if (sum(xi * (can.x.beta - thresh) > 1) > 0) {  # numerical stability
      can.lly <- -Inf
    } else {
      can.z       <- getZ(xi=xi, x.beta=can.x.beta, thresh=thresh)
      can.z.star  <- can.z^(alpha.inv)
      can.lly     <- logLikeY(y=y, theta.star=theta.star, z.star=can.z.star)
    }

    # if (sum(is.nan(can.lly)) > 0) {
    #   print("z in beta update")
    #   print(mean(can.z))
    #   print(mean(can.z.star))
    # }

    R <- sum(can.lly - cur.lly) +
         dnorm(can.beta, beta.m, beta.s, log=TRUE) -
         dnorm(beta[p], beta.m, beta.s, log=TRUE)

    if (!is.na(R)) { if (log(runif(1)) < R) {
      beta[p] <- can.beta
      x.beta  <- can.x.beta
      z       <- can.z
      z       <- can.z.star
      cur.lly <- can.lly
      acc[p]  <- acc[p] + 1
    }}
  }

  results <- list(beta=beta, x.beta=x.beta, z=z, z.star=z.star, cur.lly=cur.lly,
                  att=att, acc=acc)
  return(results)
}

# update the xi term for the linear model
updateXi <- function(y, theta.star, alpha, z, z.star, x.beta, xi, xi.m, xi.s,
                     cur.lly, acc, att, mh, thresh=0) {
  nt  <- ncol(y)
  att <- att + 1
  alpha.inv <- 1 / alpha

  can.xi     <- rnorm(1, xi, mh)

  if (sum(can.xi * (x.beta - thresh) > 1) > 0) {  # numerical stability
    # print("should generate NaNs")
    can.lly <- -Inf
  } else {
    can.z      <- getZ(xi=can.xi, x.beta=x.beta, thresh=thresh)
    can.z.star <- can.z^(alpha.inv)

    can.lly <- logLikeY(y=y, theta.star=theta.star, z.star=can.z.star)

    # if (sum(is.nan(can.lly)) > 0) {
    #   print("z in xi update")
    #   print(mean(can.z))
    #   print(mean(can.z.star))
    # }
  }

  R <- sum(can.lly - cur.lly) +
       dnorm(can.xi, xi.m, xi.s, log=TRUE) -
       dnorm(xi, xi.m, xi.s, log=TRUE)
       # dTNorm(xi, can.xi, mh, log=TRUE) -  # truncated so not symmetric
       # dTNorm(can.xi, xi, mh, log=TRUE)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    xi      <- can.xi
    z       <- can.z
    z.star  <- can.z.star
    cur.lly <- can.lly
    acc     <- acc + 1
  }}

  results <- list(xi=xi, z=z, z.star=z.star, cur.lly=cur.lly, att=att, acc=acc)
  return(results)
}

updateXiBeta2 <- function(y, theta.star, alpha, z, z.star, beta, beta.m, beta.s,
                         x.beta, xi, x, xi.m, xi.s, cur.lly,
                         can.mean, can.var,
                         acc.beta, att.beta, mh.beta,
                         acc.xi, att.xi, mh.xi, thresh=0) {
  np  <- length(beta) + 1  # using later to get the end of the vector
  ns  <- nrow(y)
  nt  <- ncol(y)
  att.beta <- att.beta + 1
  att.xi   <- att.xi + 1
  alpha.inv <- 1 / alpha
  
  cor <- diag(rep(1, np))
  cor[1, 2] <- cor[2, 1] <- -0.5
  cov <- diag(c(mh.xi, mh.beta)) %*% cor %*% diag(c(mh.xi, mh.beta))
  can.xibeta <- c(xi, beta) + t(chol(cov)) %*% rnorm(np)
  
  can.xi   <- can.xibeta[1]
  can.beta <- can.xibeta[2:np]
  
  can.x.beta <- matrix(NA, ns, nt)
  if (nt == 1) {
    if (np == 1) {
      can.x.beta <- x * can.beta
    } else {
      can.x.beta <- x %*% can.beta
    }
  } else {
    for (t in 1:nt) {
      start <- (t - 1) * ns + 1
      end   <- t * ns
      if (np == 1) {
        can.x.beta[, t] <- x[start:end, t] * can.beta
      } else {
        can.x.beta[, t] <- x[start:end, t] %*% can.beta
      }
    }
  }
  
  if (any(can.xi * (can.x.beta - thresh) > 1)) {
    can.lly <- -Inf
  } else {
    can.z <- getZ(xi=can.xi, x.beta=can.x.beta, thresh=thresh)
    can.z.star <- can.z^alpha.inv
    can.lly <- logLikeY(y=y, theta.star=theta.star, z.star=can.z.star)
  }
  
  R <- sum(can.lly - cur.lly) +
       sum(dnorm(can.beta, beta.m, beta.s, log=TRUE)) -
       sum(dnorm(beta, beta.m, beta.s, log=TRUE)) + 
       dnorm(can.xi, xi.m, xi.s, log=TRUE) -
       dnorm(xi, xi.m, xi.s, log=TRUE)
  
  if (!is.na(R)) { if (log(runif(1)) < R) {
    beta     <- can.beta
    x.beta   <- can.x.beta
    xi       <- can.xi
    z        <- can.z
    z.star   <- can.z.star
    cur.lly  <- can.lly
    acc.beta <- acc.beta + 1
    acc.xi   <- acc.xi + 1
  }}
  
  results <- list(beta=beta, x.beta=x.beta,
                  xi=xi, z=z, z.star=z.star, cur.lly=cur.lly,
                  att.beta=att.beta, acc.beta=acc.beta,
                  att.xi=att.xi, acc.xi=acc.xi)
  return(results)
}

updateXiBeta <- function(y, theta.star, alpha, z, z.star, beta, 
                         x.beta, xi, x, xt, xtx.inv, xi.m, xi.s, cur.lly,
                         acc.p, att.p, mh.p, acc.xi, att.xi, mh.xi, thresh=0) {
  # xtx.inv = solve(t(x) %*% x)
  # xt = t(x)
  
  # z should be ns x nt
  
  att.p  <- att.p + 1
  att.xi <- att.xi + 1
  
  ns <- nrow(y)
  nt <- ncol(y)
  np <- length(beta)
  alpha.inv <- 1 / alpha

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
  can.xi <- rnorm(1, xi, mh.xi)
  can.x.beta <- (1 - (-1 / log(can.p))^can.xi) / can.xi
  
  if (any(can.xi * (can.x.beta - thresh) > 1)) {
    can.lly <- -Inf
  } else {
    can.z <- getZ(xi=can.xi, x.beta=can.x.beta, thresh=thresh)
    can.z.star <- can.z^alpha.inv
    can.lly <- logLikeY(y=y, theta.star=theta.star, z.star=can.z.star)
  }
  
  R <- sum(can.lly - cur.lly) +
       dbeta(can.p[1], 1, 1, log=TRUE) - dbeta(cur.p[1], 1, 1, log=TRUE) +
       dnorm(can.xi, xi.m, xi.s, log=TRUE) - dnorm(xi, xi.m, xi.s, log=TRUE)
  
  if (!is.na(R)) { if (log(runif(1)) < R) {
    beta    <- xtx.inv %*% (xt %*% can.x.beta)
    x.beta  <- can.x.beta
    xi      <- can.xi
    z       <- can.z
    z.star  <- can.z.star
    cur.lly <- can.lly
    acc.xi  <- acc.xi + 1
    acc.p   <- acc.p + 1
  }}
  
  results <- list(beta=beta, x.beta=x.beta, xi=xi, z=z, z.star=z.star,
                  cur.lly=cur.lly, att.p=att.p, acc.p=acc.p,
                  att.xi=att.xi, acc.xi=acc.xi)
}

# updateXiBeta <- function(y, theta.star, alpha, z, z.star, beta, beta.m, beta.s,
#                          x.beta, xi, x, xi.m, xi.s, cur.lly,
#                          can.mean, can.var,
#                          acc.beta, att.beta, mh.beta,
#                          acc.xi, att.xi, mh.xi, thresh=0) {
#   np  <- length(beta) + 1  # using later to get the end of the vector
#   ns  <- nrow(y)
#   nt  <- ncol(y)
#   att.beta <- att.beta + 1
#   att.xi   <- att.xi + 1
#   alpha.inv <- 1 / alpha
# 
#   can.xibeta <- rmvt(n = 1, sigma = can.var, df = 15, delta = can.mean)
# 
#   can.xi   <- can.xibeta[1]
#   can.beta <- can.xibeta[2:np]
# 
#   can.x.beta <- matrix(NA, ns, nt)
#   if (nt == 1) {
#     if (np == 1) {
#       can.x.beta <- x * can.beta
#     } else {
#       can.x.beta <- x %*% can.beta
#     }
#   } else {
#     for (t in 1:nt) {
#       start <- (t - 1) * ns + 1
#       end   <- t * ns
#       if (np == 1) {
#         can.x.beta[, t] <- x[start:end, t] * can.beta
#       } else {
#         can.x.beta[, t] <- x[start:end, t] %*% can.beta
#       }
#     }
#   }
# 
#   if (any(can.xi * (can.x.beta - thresh) > 1)) {
#     can.lly <- -Inf
#   } else {
#     can.z <- getZ(xi=can.xi, x.beta=can.x.beta, thresh=thresh)
#     can.z.star <- can.z^alpha.inv
#     can.lly <- logLikeY(y=y, theta.star=theta.star, z.star=can.z.star)
#   }
# 
#   R <- sum(can.lly - cur.lly) +
#        sum(dnorm(can.beta, beta.m, beta.s, log=TRUE)) -
#        sum(dnorm(beta, beta.m, beta.s, log=TRUE)) +
#        dnorm(can.xi, xi.m, xi.s, log=TRUE) -
#        dnorm(xi, xi.m, xi.s, log=TRUE)
# 
#   if (!is.na(R)) { if (log(runif(1)) < R) {
#     beta     <- can.beta
#     x.beta   <- can.x.beta
#     xi       <- can.xi
#     z        <- can.z
#     z.star   <- can.z.star
#     cur.lly  <- can.lly
#     acc.beta <- acc.beta + 1
#     acc.xi   <- acc.xi + 1
#   }}
# 
#   results <- list(beta=beta, x.beta=x.beta,
#                   xi=xi, z=z, z.star=z.star, cur.lly=cur.lly,
#                   att.beta=att.beta, acc.beta=acc.beta,
#                   att.xi=att.xi, acc.xi=acc.xi)
#   return(results)
# }

# update the random effects for theta.star
updateA <- function(y, kernel, a, alpha, cur.lly, cur.llps, wz.star,
                    mid.points, bin.width, mh, cuts) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  for (t in 1:nt) {
    cur.lly.t <- cur.lly[, t]
    for (k in 1:nknots) {
      cur.a          <- a[k, t]
      l1             <- get.level(cur.a, cuts)  # keeps sd reasonable
      can.a          <- exp(rnorm(1, log(cur.a), mh[l1]))
      l2             <- get.level(can.a, cuts)
      # can.theta.star only changes at a site when it's near the knot
      can.kernel <- kernel[, t] + wz.star[, k, t] * (can.a - cur.a)
      if (sum(can.theta.star <= 0) > 0) {  # numerical stability
        can.llps  <- -Inf
        can.lly.t <- -Inf
      } else {
        can.llps  <- dPS.Rcpp(a=can.a, alpha=alpha,
                              mid.points=mid.points, bin.width=bin.width)
        can.lly.t <- logLikeY2(y=y[, t], kernel = can.kernel)
      }

      R <- sum(can.lly.t - cur.lly.t) +
           can.llps - cur.llps[k, t] +
           dlognormal(cur.a, can.a, mh[l2]) - # candidate sd changes
           dlognormal(can.a, cur.a, mh[l1])

      if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
        a[k, t]         <- can.a
        kernel[, t]     <- can.kernel
        cur.lly.t       <- can.lly.t
        cur.llps[k, t]  <- can.llps
      }}
    }
  }

  cur.lly <- logLikeY(y=y, theta.star=theta.star, z.star=z.star)

  if (sum(is.nan(cur.lly)) > 0) {
    cur.lly <<- cur.lly
    theta.star <<- theta.star
    z.star <<- z.star
    a <<- a
    print("NaN error in lly for a update")
  }

  results <- list(a = a, kernel = kernel, 
                  cur.lly = cur.lly, cur.llps = cur.llps)
  return(results)
}

# update the alpha term for theta.star
updateAlpha <- function(y, theta.star, a, alpha, cur.lly, cur.llps, z, z.star,
                        alpha.a, alpha.b, w, w.star, mid.points, bin.width,
                        acc, att, mh, iter) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  att <- att + 1

  # not exactly U(0, 1) - using truncation for numerical stability
  # can.alpha.star <- rnorm(1, log(alpha / (1 - alpha)), mh)
  # can.alpha      <- exp(can.alpha.star) / (1 + exp(can.alpha.star))
  alpha.star     <- qnorm(alpha)
  can.alpha.star <- rnorm(1, alpha.star, mh)
  can.alpha      <- pnorm(can.alpha.star)
  can.alpha.inv  <- 1 / can.alpha
  can.w.star     <- w^can.alpha.inv
  can.z.star     <- z^can.alpha.inv
  can.theta.star <- getThetaStar(can.w.star, a)
  can.lly        <- logLikeY(y=y, theta.star=can.theta.star, z.star=can.z.star)
  can.llps       <- matrix(NA, nknots, nt)
  can.llps       <- dPS.Rcpp(a, can.alpha, mid.points, bin.width)
  
  can.kernel <- kernelCalc(z = z, w = w, a = a, alpha = can.alpha)
  cur.kernel <- kernelCalc(z = z, w = w, a = a, alpha = alpha)
  cur.lly    <- logLikeY2(y = y, kernel = cur.kernel)
  can.lly    <- logLikeY2(y = y, kernel = can.kernel)

  if (sum(is.nan(can.w.star)) > 0) {
    print(can.alpha)
    print(w)
    stop("NaN error in can.w.star for alpha update")
  }

  if (sum(is.nan(can.z.star)) > 0) {
    print(can.alpha)
    print(z)
    stop("NaN error in can.z.star for alpha update")
  }

  if (sum(is.nan(can.lly)) > 0) {
    can.lly <<- can.lly
    can.theta.star <<- can.theta.star
    can.w.star <<- can.w.star
    can.z.star <<- can.z.star
    can.theta <<- can.theta.star^can.alpha
    can.alpha <<- can.alpha
    z <<- z
    a <<- a
    stop("NaN error in can.lly for alpha update")
  }

  # for(t in 1:nt) {
  #   for (k in 1:nknots) {
  #     can.llps[k, t] <- dPS(a[k, t], can.alpha, mid.points, bin.width)
  #   }
  # }

  R <- sum(can.lly - cur.lly) + sum(can.llps - cur.llps) +
       # dnorm(can.alpha.star, log=TRUE) - dnorm(alpha.star, log=TRUE)
       dbeta(can.alpha, alpha.a, alpha.b, log=TRUE) -
       dbeta(alpha, alpha.a, alpha.b, log=TRUE)

#   if ((iter > 3000) & (iter %% 50 == 0)) {
#     can.lly <- logLikeY(y = y, theta.star = can.theta.star, z.star = can.z.star, 
#                         print = TRUE)
#     cur.lly <- logLikeY(y = y, theta.star = theta.star, z.star = z.star, 
#                         print = TRUE)
#     cur.llps <- dPS.Rcpp(a, alpha, mid.points, bin.width)
#     print(paste("can.alpha = ", can.alpha, sep=""))
#     print(paste("alpha = ", alpha, sep=""))
#     print(paste("lldiff = ", sum(can.lly - cur.lly), sep=""))
#     print(paste("alldiff = ", sum(can.llps - cur.llps), sep=""))
#     print(paste("R = ", R, sep=""))
#   }

  if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
    alpha      <- can.alpha
    w.star     <- w^can.alpha.inv
    z.star     <- z^can.alpha.inv
    theta.star <- can.theta.star
    cur.lly    <- can.lly
    cur.llps   <- can.llps
    acc        <- acc + 1
  }}

  results <- list(alpha=alpha, w.star=w.star, z.star=z.star,
                  theta.star=theta.star, cur.lly=cur.lly, cur.llps=cur.llps,
                  att=att, acc=acc)
  return(results)
}

updateRho <- function(y, theta.star, a, alpha, cur.lly, z.star, w, w.star, dw2,
                      rho, rho.upper=Inf, acc, att, mh, iter) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  att <- att + 1
  rho.star       <- log(rho)
  can.rho.star   <- rnorm(1, rho.star, mh)
  can.rho        <- exp(can.rho.star)
  # rho.star       <- transform$probit(rho, lower=1e-6, upper=rho.upper)
  # can.rho.star   <- rnorm(1, rho.star, mh)
  # can.rho        <- transform$inv.probit(can.rho.star, lower=1e-6, upper=rho.upper)
  can.w          <- stdW(makeW(dw2=dw2, rho=can.rho))
  can.w.star     <- can.w^(1 / alpha)
  can.theta.star <- getThetaStar(w.star=can.w.star, a=a)
  can.lly <- logLikeY(y=y, theta.star=can.theta.star, z.star=z.star)

  logrho.m <- -1
  logrho.s <- 2

  R <- sum(can.lly - cur.lly) +
       # dnorm(can.rho.star, log=TRUE) - dnorm(rho.star, log=TRUE)
       dnorm(can.rho.star, logrho.m, logrho.s, log=TRUE) -
       dnorm(rho.star, logrho.m, logrho.s, log=TRUE)
       # dgamma(can.rho, 1, 1, log=T) -
       # dgamma(rho, 1, 1, log=T)

  # if ((iter > 3000) & (iter %% 50 == 0)) {
  #   cur.lly <- logLikeY(y=y, theta.star=theta.star, z.star=z.star)
  #   print(paste("can.rho = ", can.rho, sep=""))
  #   print(paste("rho = ", rho, sep=""))
  #   print(paste("lldiff = ", sum(can.lly - cur.lly), sep=""))
  #   print(paste("R = ", R, sep=""))
  # }

  if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
    rho        <- can.rho
    w          <- can.w
    w.star     <- can.w.star
    theta.star <- can.theta.star
    cur.lly    <- can.lly
    acc        <- acc + 1
  }}

  results <- list(rho=rho, w=w, w.star=w.star, theta.star=theta.star,
                  cur.lly=cur.lly, att=att, acc=acc)
  return(results)
}

pred.spgev <- function(mcmcoutput, s.pred, x.pred, knots, start=1, end=NULL,
                        thin=1, thresh=0, update=NULL) {
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
  beta  <- matrix(mcmcoutput$beta[start:end, , drop = F], niters, p)
  xi    <- mcmcoutput$xi[start:end]
  a     <- mcmcoutput$a[start:end, , , drop = F]
  alpha <- mcmcoutput$alpha[start:end]
  rho   <- mcmcoutput$rho[start:end]
  
  dw2p  <- as.matrix(rdist(s.pred, knots))^2
  prob.success <- matrix(NA, nrow=niters, ncol=np)
  x.beta <- matrix(NA, np, nt)

  for (i in 1:length(start:end)) {
    x.beta       <- getXBeta(x.pred, ns = np, nt = nt, beta = beta[i, ])
    z            <- getZ(xi = xi[i], x.beta=x.beta, thresh=thresh)
    z.star       <- z^(1 / alpha[i])
    w            <- stdW(makeW(dw2 = dw2p, rho = rho[i]))
    w.star       <- w^(1 / alpha[i])
    theta.star   <- getThetaStar(w.star = w.star, a = a[i, , ])
    theta.z.star <- -theta.star / z.star
    prob.success[i, ] <- 1 - exp(theta.z.star)

    # if z is nan, it means that x.beta is such a large number there is
    # basically 0 probability that z < 0
    if (sum(is.nan(z)) > 0) {
      these <- which(is.nan(z))
      prob.success[i, these] <- 1
    }

    if (!is.null(update)) {
      if (i %% update == 0) {
        cat("\t Iter", i, "\n")
      }
    }
  }

  return(prob.success)

}

predictY <- function(mcmcoutput, s.pred, x.pred, knots, start=1, end=NULL,
                     thin=1, thresh=0, update=NULL) {

  if (is.null(end)) {
    end <- nrow(mcmcoutput$beta)
  }

  np   <- nrow(s.pred)
  iters <- end - start + 1
  dw2p <- as.matrix(rdist(s.pred, knots))^2
  yp <- matrix(NA, nrow=iters, ncol=np)
  for (i in start:end) {
    yp[i, ] <- rRareBinarySpat(x=x.pred, s=s.pred, knots=knots,
                             beta=mcmcoutput$beta[i, ], xi=mcmcoutput$xi[i],
                             alpha=mcmcoutput$alpha[i], rho=mcmcoutput$rho[i],
                             thresh=thresh, dw2=dw2p, a=mcmcoutput$a[i, ])$y
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