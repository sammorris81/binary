pairwise.rarebinaryR <- function(par, y, dw2, cov) {
  # par: parameter vector (xi, alpha, rho, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  xi    <- par[1]
  alpha <- par[2]
  rho   <- par[3]
  beta  <- par[4:npars]

  ns     <- length(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  if (length(beta) == 1) {
    x.beta <- cov * beta
  } else {
    x.beta <- cov %*% beta           # should be ns x np
  }
  z      <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (alpha < 0 | alpha > 1) {
    return(1e99)
  }
  if (rho < 0) {
    return(1e99)
  }
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)
  print(kernel[1:5, 1:5])

  ll <- 0
  for (i in 1:(ns - 1)) {
    for (j in (i + 1):ns) {
      # in all cases, we need the joint
      # first take the kernel bits and add i and j for each knot
      # joint <- -sum(apply(kernel[c(i, j), ], 2, sum)^alpha)
      joint <- -getJoint2(kernel = kernel[c(i, j), ], alpha = alpha)
      if (y[i] == 0 & y[j] == 0) {
        ll <- ll + joint
      } else if (y[i] == 1 & y[j] == 0) {  # add in marginal for i
        ll <- ll + log(exp(-1 / z[j]) - exp(joint))
      } else if (y[i] == 0 & y[j] == 1) {
        ll <- ll + log(exp(-1 / z[i]) - exp(joint))
      } else if (y[i] == 1 & y[j] == 1) {  # do 1,1 likelihood
        ll <- ll + log(1 - exp(-1 / z[i]) - exp(-1 / z[j]) + exp(joint))
      }
    }
  }

  return(-ll)
}

pairwise.rarebinaryCPP.1 <- function(par, y,
                                     dw2, d, cov, max.dist,
                                     alpha.min = 0.2, alpha.max = 0.9,
                                     threads = 1) {
  # par: parameter vector (xi, alpha, rho, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  xi         <- par[1]
  alpha.star <- par[2]
  rho.star   <- par[3]
  beta       <- par[4:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  # transform to actual space
  # alpha in (alpha.min, alpha.max)
  alpha <- (exp(alpha.star) / (1 + exp(alpha.star))) *
    (alpha.max - alpha.min) + alpha.min
  rho <- exp(rho.star)

  ns     <- nrow(y)
  nt     <- ncol(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  if (length(beta) == 1) {
    x.beta <- matrix(cov * beta, ns, nt)
  } else {
    x.beta <- cov %*% beta           # should be ns long
  }
  z <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }
  if (rho > 0.2 * max(d)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    rho = rho, d = d, max_dist = max(d), threads = threads)

  return(-ll)
}

pairwise.rarebinaryCPP.2 <- function(par, y, xi,
                                     dw2, d, cov, max.dist,
                                     alpha.min = 0.2, alpha.max = 0.9,
                                     threads = 1) {
  # par: parameter vector (alpha, rho, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  alpha.star <- par[1]
  rho.star   <- par[2]
  beta       <- par[3:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  # transform to actual space
  # alpha in (alpha.min, alpha.max)
  alpha <- (exp(alpha.star) / (1 + exp(alpha.star))) *
    (alpha.max - alpha.min) + alpha.min
  rho <- exp(rho.star)

  ns     <- nrow(y)
  nt     <- ncol(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  if (length(beta) == 1) {
    x.beta <- matrix(cov * beta, ns, nt)
  } else {
    x.beta <- cov %*% beta           # should be ns long
  }
  z <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }
  if (rho > 0.2 * max(d)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    rho = rho, d = d, max_dist = max(d), threads = threads)

  return(-ll)
}

pairwise.rarebinaryCPP.3 <- function(par, y, alpha.star,
                                     dw2, d, cov, max.dist,
                                     alpha.min = 0.2, alpha.max = 0.9,
                                     threads = 1) {
  # par: parameter vector (xi, rho, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  xi         <- par[1]
  rho.star   <- par[2]
  beta       <- par[3:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  # transform to actual space
  # alpha in (alpha.min, alpha.max)
  alpha <- (exp(alpha.star) / (1 + exp(alpha.star))) *
    (alpha.max - alpha.min) + alpha.min
  rho <- exp(rho.star)

  ns     <- nrow(y)
  nt     <- ncol(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  if (length(beta) == 1) {
    x.beta <- matrix(cov * beta, ns, nt)
  } else {
    x.beta <- cov %*% beta           # should be ns long
  }
  z <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }
  if (rho > 0.2 * max(d)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    rho = rho, d = d, max_dist = max(d), threads = threads)

  return(-ll)
}

pairwise.rarebinaryCPP.4 <- function(par, y, rho.star,
                                     dw2, d, cov, max.dist,
                                     alpha.min = 0.2, alpha.max = 0.9,
                                     threads = 1) {
  # par: parameter vector (xi, alpha, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  xi         <- par[1]
  alpha.star <- par[2]
  beta       <- par[3:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  # transform to actual space
  # alpha in (alpha.min, alpha.max)
  alpha <- (exp(alpha.star) / (1 + exp(alpha.star))) *
    (alpha.max - alpha.min) + alpha.min
  rho <- exp(rho.star)

  ns     <- nrow(y)
  nt     <- ncol(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  if (length(beta) == 1) {
    x.beta <- matrix(cov * beta, ns, nt)
  } else {
    x.beta <- cov %*% beta           # should be ns long
  }
  z <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }
  if (rho > 0.2 * max(d)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    rho = rho, d = d, max_dist = max(d), threads = threads)

  return(-ll)
}

pairwise.rarebinaryCPP.5 <- function(par, y, xi, alpha.star,
                                     dw2, d, cov, max.dist,
                                     alpha.min = 0.2, alpha.max = 0.9,
                                     threads = 1) {
  # par: parameter vector (rho, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  rho.star   <- par[1]
  beta       <- par[2:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  # transform to actual space
  # alpha in (alpha.min, alpha.max)
  alpha <- (exp(alpha.star) / (1 + exp(alpha.star))) *
    (alpha.max - alpha.min) + alpha.min
  rho <- exp(rho.star)

  ns     <- nrow(y)
  nt     <- ncol(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  if (length(beta) == 1) {
    x.beta <- matrix(cov * beta, ns, nt)
  } else {
    x.beta <- cov %*% beta           # should be ns long
  }
  z <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }
  if (rho > 0.2 * max(d)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    rho = rho, d = d, max_dist = max(d), threads = threads)

  return(-ll)
}

pairwise.rarebinaryCPP.6 <- function(par, y, xi, rho.star,
                                     dw2, d, cov, max.dist,
                                     alpha.min = 0.2, alpha.max = 0.9,
                                     threads = 1) {
  # par: parameter vector (alpha, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  alpha.star <- par[1]
  beta       <- par[2:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  # transform to actual space
  # alpha in (alpha.min, alpha.max)
  alpha <- (exp(alpha.star) / (1 + exp(alpha.star))) *
    (alpha.max - alpha.min) + alpha.min
  rho <- exp(rho.star)

  ns     <- nrow(y)
  nt     <- ncol(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  if (length(beta) == 1) {
    x.beta <- matrix(cov * beta, ns, nt)
  } else {
    x.beta <- cov %*% beta           # should be ns long
  }
  z <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }
  if (rho > 0.2 * max(d)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    rho = rho, d = d, max_dist = max(d), threads = threads)

  return(-ll)
}

pairwise.rarebinaryCPP.7 <- function(par, y, alpha.star, rho.star,
                                     dw2, d, cov, max.dist,
                                     alpha.min = 0.2, alpha.max = 0.9,
                                     threads = 1) {
  # par: parameter vector (xi, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  xi         <- par[1]
  beta       <- par[2:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  # transform to actual space
  # alpha in (alpha.min, alpha.max)
  alpha <- (exp(alpha.star) / (1 + exp(alpha.star))) *
    (alpha.max - alpha.min) + alpha.min
  rho <- exp(rho.star)

  ns     <- nrow(y)
  nt     <- ncol(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  if (length(beta) == 1) {
    x.beta <- matrix(cov * beta, ns, nt)
  } else {
    x.beta <- cov %*% beta           # should be ns long
  }
  z <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }
  if (rho > 0.2 * max(d)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    rho = rho, d = d, max_dist = max(d), threads = threads)

  return(-ll)
}

pairwise.rarebinaryCPP.8 <- function(par, y, xi, alpha.star, rho.star,
                                     dw2, d, cov, max.dist,
                                     alpha.min = 0.2, alpha.max = 0.9,
                                     threads = 1) {
  # par: parameter vector (beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  beta  <- par[1:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  # transform to actual space
  # alpha in (alpha.min, alpha.max)
  alpha <- (exp(alpha.star) / (1 + exp(alpha.star))) *
    (alpha.max - alpha.min) + alpha.min
  rho <- exp(rho.star)

  ns     <- nrow(y)
  nt     <- ncol(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  if (length(beta) == 1) {
    x.beta <- matrix(cov * beta, ns, nt)
  } else {
    x.beta <- cov %*% beta           # should be ns long
  }
  z <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }
  if (rho > 0.2 * max(d)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    rho = rho, d = d, max_dist = max(d), threads = threads)

  return(-ll)
}