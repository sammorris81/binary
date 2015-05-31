# Rcpp functions
if (!exists("dPS.Rcpp")) {
  source('llps_cpp.R')
}

################################################################################
# Common data transformations
################################################################################
transform <- list(
  logit = function(x, lower=0, upper=1) {
    x <- (x - lower) / (upper - lower)
    return(log(x / (1 - x)))
  },
  inv.logit = function(x, lower=0, upper=1) {
    p <- exp(x) / (1 + exp(x))
    p <- p * (upper - lower) + lower
    return(p)
  },
  probit = function(x, lower=0, upper=1) {
    x <- (x - lower) / (upper - lower)
    return(qnorm(x))
  },
  inv.probit = function(x, lower=0, upper=1) {
    p <- pnorm(x)
    p <- p * (upper - lower) + lower
    return(p)
  },
  log = function(x) log(x),
  exp = function(x) exp(x),
  copula = function(dens) {
    this.dens <- paste("p", dens, sep="")
    function(x, ...) qnorm(do.call(this.dens, args=list(x, ...)))
  },
  inv.copula = function(dens) {
    this.dens <- paste("q", dens, sep="")
    function(x, ...) do.call(this.dens, args=list(pnorm(x), ...))
  }
)

# TODO: There needs to be more checks here.
# (1 + xi * (thresh - x.beta)) > 0 for all x and x.pred
getZ <- function(xi, x.beta, thresh=0) {
  if (xi != 0) {
    z <- (1 + xi * (thresh - x.beta))^(1 / xi)
  } else {
    z <- exp(thresh - x.beta)
  }

  return(z)
}

# theta.star = theta^(1 / alpha) = sum_l=1^L a_l * w_l^(1 / alpha)
getThetaStar <- function(w.star, a) {
  # theta.star is ns x nt
  # w.star = w^(1 / alpha) is ns x nknots
  # a is nknots x nt
  # alpha in (0,1)
  if (length(a) == 1) {theta.star <- w.star * a}
  if (length(a) > 1) {theta.star <- w.star %*% a}
  return(theta.star)
}

# get the kernel weighting
makeW <- function(dw2, rho) {
  w <- exp(-0.5 * dw2 / (rho^2))
  return(w)
}

# standardize the kernel weights
stdW <- function(x, single=FALSE) {
  if (single) {x <- x / sum(x)}
  if (!single) {x <- sweep(x, 1, rowSums(x), "/")}
  return(x)
}


# get find the ll for y - returns ns x nt matrix for each site/day
logLikeY <- function(y, theta.star, z.star) {
  theta.z.star <- -theta.star / z.star
  if (!is.null(dim(y))) {
    ll.y <- matrix(-Inf, dim(y))
  } else {
    ll.y <- rep(-Inf, length(y))
  }

  # numerical stability issue. originally was using
  # (1 - y) * P(Y = 0) + y * P(Y = 1)
  # would return NaN because 0 * -Inf is not a number
  these <- which(y == 1)
  ll.y[-these] <- theta.z.star[-these]
  ll.y[these]  <- log(1 - exp(theta.z.star[these]))
  # if (sum(is.nan(ll.y)) > 0) {
  #   these <- which(is.nan(ll.y))
  #   print(theta.z.star[these])
  #   # print(1 - exp(theta.z.star[these]))
  # }
  return(ll.y)
}

################################################################################
#### positive stable density functions
################################################################################
dPS <- function(a, alpha, mid.points, bin.width) {
  l <- -Inf

  if (a > 0) {
    l <- log(sum(bin.width * ld(mid.points, a, alpha)))
  }

  return(l)
}

# used when evaluating the postive stable density
ld <- function(u, a, alpha) {
  psi <- pi * u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)
  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * log(a) +
          log(c) - c * (1 / a^(alpha / (1 - alpha)))

  return(exp(logd))
}

rPS <- function(n, alpha) {
  #### PS(alpha) generation as given by Stephenson(2003)
  unif <- runif(n) * pi
  stdexp.ps <- rexp(n, 1)
  logs <- (1 - alpha) / alpha * log(sin((1 - alpha) * unif)) +
          log(sin(alpha * unif)) - (1 - alpha) / alpha * log(stdexp.ps) -
          1 / alpha * log(sin(unif))
  return(exp(logs))
}

# generate rare binary data
# for now working with independent to test functions for beta and xi
rRareBinaryInd <- function(x, beta, xi, prob.success = 0.05) {
  nt <- dim(x)[2]
  ns <- dim(x)[1]

  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  z <- matrix(rgev(n=ns * nt, 1, 1, 1), ns, nt)
  h <- x.beta + (z^xi - 1) / xi

  # set the threshold for success at whatever value will give us
  # our desired percentage of 1s.
  thresh <- quantile(h, probs = prob.success)
  y <- ifelse(h > thresh, 1, 0)

  results <- list(y = y)
  return(results)
}

# generate dependent rare binary data
rRareBinarySpat <- function(x, s, knots, beta, xi, alpha, rho,
                            prob.success = 0.05, dw2 = NULL, a = NULL) {
  ns     <- dim(x)[1]
  nt     <- dim(x)[2]
  p      <- length(beta)
  y      <- matrix(NA, ns, nt)
  nknots <- nrow(knots)

  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    if (p > 1) {
      x.beta[, t] <- x[, t, ] %*% beta
    } else {
      x.beta[, t] <- x[, t] * beta
    }
  }

  # get weights
  if (is.null(dw2)) {  # for predictions, already have dw2
    dw2    <- as.matrix(rdist(s, knots))^2  # dw2 is ns x nknots
  }
  w      <- stdW(makeW(dw2 = dw2, rho = rho))  # w is ns x nknots
  w.star <- w^(1 / alpha)

  # get random effects and theta.star
  if (is.null(a)) {
    a <- matrix(rPS(n = nknots * nt, alpha = alpha), nknots, nt)
  }
  theta.star <- getThetaStar(w.star = w.star, a = a)

  # get underlying latent variable
  u <- matrix(rgev(n = ns * nt, 1, alpha, alpha), ns, nt)
  z <- u * theta.star^alpha

  h <- x.beta + (z^xi - 1) / xi

  # set the threshold for success at whatever value will give us
  # our desired percentage of 1s.
  thresh <- quantile(h, probs = (1 - prob.success))

  y <- ifelse(h > thresh, 1, 0)

  results <- list(y = y, a = a, thresh = thresh)
  return(results)
}

# sigma.sq is partial sill
rLogitSpat <- function(x, s, knots, beta, rho, sigma.sq, nu = 0.5) {
  ns     <- nrow(s)
  nt     <- dim(x)[2]
  p      <- length(beta)
  y      <- matrix(NA, ns, nt)
  nknots <- nrow(knots)

  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    if (p > 1) {
      x.beta[, t] <- x[, t, ] %*% beta
    } else {
      x.beta[, t] <- x[, t] * beta
    }
  }

  # we use 22 for covariance between knots
  # we use 12 for covariance between all sites and knots
  d.22 <- as.matrix(rdist(knots))  # d.knots is nknots x nknots
  d.12 <- as.matrix(rdist(s, knots))  #d.sknots is ns x nknots
  diag(d.22) <- 0

  if (nu == 0.5) {
    Sigma.12 <- simple.cov.sp(D = d.12, sp.type = "exponential",
                              sp.par = c(sigma.sq, rho), error.var = 0,
                              finescale.var = 0)
    Sigma.22 <- simple.cov.sp(D = d.22, sp.type = "exponential",
                              sp.par = c(sigma.sq, rho), error.var = 0,
                              finescale.var = 0)
  } else {
    Sigma.12 <- simple.cov.sp(D = d.12, sp.type = "matern",
                              sp.par = c(sigma.sq, rho), error.var = 0,
                              smoothness = nu, finescale.var = 0)
    Sigma.22 <- simple.cov.sp(D = d.22, sp.type = "matern",
                              sp.par = c(sigma.sq, rho), error.var = 0,
                              smoothness = nu, finescale.var = 0)
  }

  Sigma.22.inv <- chol2inv(chol(Sigma.22))
  Sigma.12.22.inv <- Sigma.12 %*% Sigma.22.inv
  sd.mtx <- t(chol(Sigma.22))

  # generate the random effects
  w <- matrix(NA, nknots, nt)
  w.tilde <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    w[, t] <- sd.mtx %*% rnorm(nknots, 0, 1)
    w.tilde[, t] <- Sigma.12.22.inv %*% w[, t]
  }

  p.success <- 1 / (1 + exp(-(x.beta + w.tilde)))
  for (t in 1:nt) {
    y[, t] <- rbinom(n = ns, size = 1, prob = p.success[, t])
  }

  results <- list(y = y, w.tilde = w.tilde)

  return(results)
}

# update mh settings
mhUpdate <- function(acc, att, mh, nattempts = 50, lower = 0.8, higher = 1.2) {
  acc.rate     <- acc / att
  these.update <- att > nattempts
  these.low    <- (acc.rate < 0.25) & these.update
  these.high   <- (acc.rate > 0.50) & these.update

  mh[these.low]  <- mh[these.low] * lower
  mh[these.high] <- mh[these.high] * higher

  acc[these.update] <- 0
  att[these.update] <- 0

  results <- list(acc=acc, att=att, mh=mh)
  return(results)
}


# ECkern <- function(h, alpha, gamma, Lmax=50) {
#   dw2 <- rdist(c(0, h), seq(-Lmax, Lmax, 1))
#   W <- fac2FAC(make.fac(dw2, gamma))^(1 / alpha)
#   for (j in 1:length(h)) {
#     h[j]<-sum((W[1, ] + W[j + 1, ])^alpha)
#   }

#   return(h)
# }



# ld2 is used in generating positive stable random variables
ld2 <- function(u, logs, alpha, shift=0, log=TRUE) {

  logs <- logs - shift / alpha
  s <- exp(logs)
  psi <- pi * u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)

  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * logs +
          log(c) - c * (1 / s^(alpha / (1 - alpha))) +
          logs

  return(logd)
}

dlognormal <- function(x, mu, sig) {
  dnorm(log(x), log(mu), sig, log=T) - log(x)
}

# logd <- function(theta, v) {
#   sum(log(theta) - theta * v)
# }

#### Used to set the standard deviation for the candidate distribution
#### for the A terms in the random effect. When log(A) is large means
#### sd is smaller, and log(A) small means sd is larger.
get.level.1 <- function(a, cuts) {
  lev <- a * 0 + 1
  for (j in 1:length(cuts)) {
    lev <- ifelse(a > cuts[j], j + 1, lev)
  }
  return(lev)
}

get.level <- function(a, cuts) {
  if (length(a) > 1) {
    warning("get.level should only be used when length(a) = 1")
  }
  lev <- sum(a > cuts) + 1
  return(lev)
}

logdet <- function(X) {
  determinant(X)$modulus
}

trunc <- function(x, eps=0.1) {
  x <- ifelse(x < eps, eps, x)
  x <- ifelse(x > 1-eps, 1-eps, x)
  return(x)
}

#########################################################################
# Arguments:
#   mn(nt): mean
#   sd(nt, nt): standard deviation
#   lower(1): lower truncation point (default=-Inf)
#   upper(1): upper truncation point (default=Inf)
# fudge(1): small number for numerical stability (used for lower bound)
#
# Returns:
#   y(nt): truncated normal data
#########################################################################
rTNorm <- function(mn, sd, lower=-Inf, upper=Inf, fudge=0) {
  lower.u <- pnorm(lower, mn, sd)
  upper.u <- pnorm(upper, mn, sd)

  # replace <- ((mn / sd) > 5) & (lower == 0)
  # lower.u[replace] <- 0
  lower.u <- ifelse( mn / sd > 5 & lower == 0, 0, lower.u )
  U <- tryCatch(runif(length(mn), lower.u, upper.u),
                warning=function(e) {
                  cat("mn =", mn, "\n")
                  cat("sd =", sd, "\n")
                  cat("lower.u =", lower.u, "\n")
                  cat("upper.u =", upper.u, "\n")
                })
  y <- qnorm(U, mn, sd)

  return(y)
}

dTNorm <- function(y, mn, sd, lower=-Inf, upper=Inf, fudge=0, log=TRUE) {
  lower.u <- pnorm(lower, mn, sd)
  upper.u <- pnorm(upper, mn, sd)

  ld <- dnorm(y, mn, sd, log=TRUE) - log(upper.u - lower.u)
  if (!log) {
    ld <- exp(ld)
  }

  return(ld)
}

################################################################
# Arguments:
#   post.prob(iters, yp): posterior probability of exceeding
#   validate(np): validation data
#
# Returns:
#   score(1): Brier score for dataset
################################################################
BrierScore <- function(post.prob, validate) {
  # iters <- nrow(post.prob)
  # np    <- ncol(post.prob)

  # scores <- rep(NA, iters)
  # for (i in 1:iters) {
  #   scores[i] <- mean((validate - post.prob[i, ])^2)
  # }
  probs <- apply(post.prob, 2, median)
  score <- mean((validate - probs)^2)

  return(score)
}


pairwise.rarebinaryR <- function(par, y, dw2, cov) {
  # par: parameter vector (alpha, rho, xi, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  alpha <- par[1]
  rho   <- par[2]
  xi    <- par[3]
  beta  <- par[4:npars]

  ns     <- length(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  x.beta <- cov %*% beta           # should be ns x np
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

pairwise.rarebinaryCPP <- function(par, y, dw2, cov, threads = 1) {
  # par: parameter vector (alpha, rho, xi, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # x: covariates (ns x np)

  npars <- length(par)
  alpha <- par[1]
  rho   <- par[2]
  xi    <- par[3]
  beta  <- par[4:npars]

  if (xi < 1e-10 & xi > -1e-10) {
    xi <- 0
  }

  ns     <- length(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots
  x.beta <- cov %*% beta           # should be ns long
  z      <- getZ(xi, x.beta)       # should be ns long

  # keep estimates in their actual space
  if (alpha < 1e-10 | alpha > (1 - 1e-10)) {
    return(1e99)
  }
  if (rho < 1e-10) {
    return(1e99)
  }
  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    threads = threads)

  return(-ll)
}

pairwise.rarebinary2CPP <- function(par, alpha, rho, y, W, cov, threads = 1) {
  # par: parameter vector (xi, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # cov: covariates (ns x np)

  npars <- length(par)
  xi    <- par[1]
  beta  <- par[2:npars]

  x.beta <- cov %*% beta           # should be ns long
  z      <- getZ(xi, x.beta)       # should be ns long

  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    threads = threads)

  return(-ll)  # optim minimizes, so return neg loglike
}

pairwise.rarebinary3CPP <- function(par, rho, y, W, cov, threads = 1) {
  # par: parameter vector (xi, alpha, beta)
  # y: observed data (ns)
  # dw2: squared distance between sites and knots (ns x nknots)
  # cov: covariates (ns x np)

  npars <- length(par)
  xi    <- par[1]
  alpha <- par[2]
  beta  <- par[3:npars]

  x.beta <- cov %*% beta           # should be ns long
  z      <- getZ(xi, x.beta)       # should be ns long

  if (any((xi * x.beta) > 1)) {
    return(1e99)
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)

  ll <- pairwiseCPP(kernel = kernel, alpha = alpha, z = z, y = y,
                    threads = threads)

  return(-ll)  # optim minimizes, so return neg loglike
}

fit.rarebinaryCPP <- function(init.par, rho, y, dw2, cov, threads = 1) {

  ns     <- length(y)
  np     <- ncol(cov)
  nknots <- ncol(dw2)

  W      <- stdW(makeW(dw2, rho))  # should be ns x nknots

  results <- optim(init.par, pairwise.rarebinary3CPP, rho = rho,
                   y = y, W = W,
                   cov = x, threads = threads,
                   # lower = c(1e-6, 1e-6, -1, -Inf),
                   # upper = c(1 - 1e-6, Inf, 3, Inf),
                   hessian = TRUE)
  # results$alpha <- alpha
  results$rho   <- rho

  return(results)
}


getJoint2 <- function(kernel, alpha) {
  knot.contribute <- colSums(kernel)^alpha
  # print(knot.contribute)
  joint <- sum(knot.contribute)
  return(joint)
}