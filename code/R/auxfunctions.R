if (file.exists("../../../usefulR/usefulfunctions.R")) {
  source("../../../usefulR/usefulfunctions.R")
}

# Rcpp functions
if (!exists("dPSCPP")) {
  sourceCpp(file = "./llps.cpp")
}

if (!exists("ifelsematCPP")) {
  sourceCpp(file = "./ifelse.cpp")
}

if (!exists("getThetaCPP")) {
  sourceCpp(file = "./getTheta.cpp")
}

if (!exists("pairwiseCPP")) {
  sourceCpp(file = "./pairwise.cpp")
}

source('hmc_aux.R')

source('pairwise.R')

# calculated values for the MCMC
getAW <- function(d, p, c, o) {
  w <- c$w
  alpha <- p$alpha
  
  return(getawCPP(a_star = p$a^p$alpha, w = c$w, alpha = p$alpha))
} 

getTheta <- function(d, p, c, o) {
  return(c$z^(-1 / p$alpha) * c$aw)
}

getXBeta <- function(d, p, c, o) {
  ns <- nrow(d$y)
  nt <- ncol(d$y)
  if (nt == 1) {
    return(d$x %*% p$beta)
  } else {
    x.beta <- matrix(NA, ns, nt)
    for (t in 1:nt) {
      start <- (t - 1) * ns + 1
      end   <- t * ns
      x.beta[, t] <- d$x[start:end, ] %*% p$beta
    }
  }
}

# (1 + xi * (thresh - x.beta)) > 0 for all x and x.pred
getZ <- function(d, p, c, o) {
  if (p$xi != 0) {
    z <- (1 + p$xi * (thresh - c$x.beta))^(1 / p$xi)
  } else {
    z <- exp(thresh - c$x.beta)
  }
  
  return(z)
}

getW <- function(d, p, c, o) {
  w <- stdW(makeW(dw2 = d$dw2, rho = p$rho, A.cutoff = o$A.cutoff))
}

# useful functions
adjustX <- function(x, y) {
  # always want a matrix that has ns * nt rows and np cols
  # we do this to simplify the multiplication of x %*% beta
  ns <- nrow(y)
  nt <- ncol(y)

  if (length(x) %% length(y) != 0) {
    stop("the number of sites in x is probably not the same as in y")
  }

  np <- length(x) / length(y)
  if (is.null(dim(x))) {
    if (length(x) != length(y)) {
      stop("x cannot be coerced to proper dimensions")
    } else {
      x <- matrix(x, nrow = ns * nt, ncol = np)
      return(x)
    }
  } else {
    if (is.na(dim(x)[3])) {  # i.e. x is a matrix
      if ((nrow(x) == (ns * nt)) & (ncol(x) == np)) {
        return(x)
      } else {
        stop("if x is a matrix, it should have ns * nt rows and np cols")
      }
    } else {
      if ((dim(x)[3] == np) & (dim(x)[2] == nt)) {
        x.temp <- matrix(NA, nrow = ns * nt, ncol = np)
        for (p in 1:np) { for (t in 1:nt) {
          start <- (t - 1) * ns + 1
          end   <- t * ns
          x.temp[start:end, p] <- x[, t, p]
        } }
        return(x.temp)
      }
    }
  }
}



# storing each day as an element of a list
getwzStar <- function(z, w, alpha) {
  nknots  <- ncol(w)
  ns      <- nrow(w)
  nt      <- ncol(z)

  wz.star <- array(NA, dim=c(ns, nknots, nt))
  for (t in 1:nt) {
    z.t            <- matrix(rep(z[, t], nknots), ns, nknots)
    wz.star[, , t] <- exp((log(w) - log(z.t)) / alpha)
    wz.star[, , t] <- ifelsematCPP(wz.star[, , t], 1e-7)
  }

  # wz.star <- getwzstarCPP(z = z, w = w, alpha = alpha)

  return(wz.star)
}

# getwz <- function(z, w) {
#   nknots  <- ncol(w)
#   ns      <- nrow(w)
#   nt      <- ncol(z)
#   
#   wz <- array(NA, dim=c(ns, nknots, nt))
#   for (t in 1:nt) {
#     wz[, , t] <- w / rep(z[, t], nknots)
#   }
#   
#   # wz.star <- getwzstarCPP(z = z, w = w, alpha = alpha)
#   
#   return(wz)
# }

# # trying a slightly different calculation
# getTheta <- function(wz.star, a, IDs = NULL) {
#   nt <- dim(wz.star)[3]
#   nknots <- dim(wz.star)[2]
#   ns <- dim(wz.star)[1]
# 
#   theta <- matrix(NA, nrow = ns, ncol = nt)
# 
#   if (is.null(IDs)) {
#     for (t in 1:nt) {
#       if (nknots == 1) {
#         theta[, t] <- wz.star[, , t] * a[ , t]
#       } else if (nknots > 1) {
#         theta[, t] <- wz.star[, , t] %*% a[, t]
#       }
#     }
#   } else {  # when using IDs, there needs to be more than one knot
#     for (t in 1:nt) {
#       a.t <- a[, t]
#       for (i in 1:ns) {
#         these <- IDs[[i]]
#         theta[i, t] <- wz.star[i, these, t] %*% a.t[these]
#       }
#     }
#   }
# 
#   return(theta)
# }

getIDs <- function(dw2, A.cutoff) {
  nknots <- ncol(dw2)
  IDs <- vector(mode = "list", length = nknots)
  for (k in 1:nknots) {
    IDs[[k]] <- which(sqrt(dw2[, k]) <= A.cutoff)
  }
  return(IDs)
}

# get the kernel weighting
makeW <- function(dw2, rho, A.cutoff = NULL) {
  if (is.null(A.cutoff)) {
    A.cutoff <- max(sqrt(dw2))
  }
  w <- exp(-0.5 * dw2 / (rho^2))
  
  # only include sites that are close to the knot
  # we need to do this in addition to the IDs because the weights over the 
  # active knots needs to sum to 1. If we don't set the knots beyond the 
  # cutoff to 0, then when we don't include them later on, the weights will 
  # sum to something slightly smaller than 1
  w[sqrt(dw2) > A.cutoff] <- 0  
  return(w)
}

# standardize the kernel weights
stdW <- function(x, single = FALSE) {
  if (single) {x <- x / sum(x)}
  if (!single) {x <- sweep(x, 1, rowSums(x), "/")}
  return(x)
}


logLikeY <- function(y, theta) {
  ll.y <- matrix(-Inf, nrow(y), ncol(y))

  # numerical stability issue. originally was using
  # (1 - y) * P(Y = 0) + y * P(Y = 1)
  # would return NaN because 0 * -Inf is not a number
  ll.y[y == 0] <- -theta[y == 0]
  ll.y[y == 1]  <- log(1 - exp(-theta[y == 1]))

  return(ll.y)
}

################################################################################
#### positive stable density functions
################################################################################
dPS <- function(a, alpha, mid.points, bin.width) {
  l <- -Inf

  if (a > 0) {
    l <- log(sum(bin.width * ld(mid.points, a, alpha, log = FALSE)))
  }

  return(l)
}

# used when evaluating the postive stable density
ld <- function(u, a, alpha, log = TRUE) {
  psi  <- pi * u
  c    <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c    <- c * sin((1 - alpha) * psi) / sin(alpha * psi)
  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * log(a) +
          log(c) - c * (1 / a^(alpha / (1 - alpha)))
  
  if (!log) {
    logd <- exp(logd)
  }
  
  return(logd)
}

# ld2 is used in generating positive stable random variables
# also used when updating A and we include auxiliary variable from 
# Stephenson
ld2 <- function(u, logs, alpha, shift = 0, log = TRUE) {
  
  logs <- logs - shift / alpha
  s <- exp(logs)
  psi <- pi * u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)
  
  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * logs +
    log(c) - c * (1 / s^(alpha / (1 - alpha))) +
    logs
  
  if (!log) {
    logd <- exp(logd)
  }
  return(logd)
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

# theta.star = theta^(1 / alpha) = sum_l=1^L a_l * w_l^(1 / alpha)
getThetaStar <- function(w.star, a) {
  # theta.star is ns x nt
  # w.star = w^(1 / alpha) is ns x nknots
  # a is nknots x nt
  # alpha in (0,1)
  w.star <- ifelse(w.star < 1e-7, 0, w.star)
  if (length(a) == 1) {theta.star <- w.star * a}
  if (length(a) > 1) {theta.star <- w.star %*% a}
  return(theta.star)
}

# generate dependent rare binary data
rRareBinarySpat <- function(x, s, knots, beta, xi, alpha, rho, nt = 1,
                            prob.success = 0.05, dw2 = NULL, a = NULL) {
  p <- length(beta)
  if (nt == 1) {
    ns <- nrow(x)
  } else {
    if (nrow(x) %% nt != 0) {
      stop("The number of rows of x must be a multiple of the number of days")
    }
    ns <- nrow(x) / nt
  }

  y      <- matrix(NA, ns, nt)
  nknots <- nrow(knots)

  x.beta <- getXBeta(x = x, ns = ns, nt = nt, beta = beta)

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

  if (xi == 0) {
    h <- x.beta + log(z)
  } else {
    h <- x.beta + (z^xi - 1) / xi
  }

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

# sigma.sq is partial sill
rProbitSpat <- function(x, s, knots, beta, rho, sigma.sq, nu = 0.5,
                        prob.success) {
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

  h <- x.beta + w.tilde

  thresh <- quantile(h, probs = (1 - prob.success))
  y <- ifelse(h > thresh, 1, 0)

  results <- list(y = y, w.tilde = w.tilde, thresh = thresh)

  return(results)
}

rHotSpotSpat <- function(x, s, xlim = NULL, ylim = NULL, bw, 
                         nhotspots = 1, 
                         prob.in = 0.99, prob.out = 0.0001) {
  ## Generate data at hotspots
  ## x: covariates
  ## s: sites
  ## xlim: limits in the x direction for the hotspot
  ## ylim: limits in the y direction for the hotspot 
  ## bw: size of the hotspots
  ## nhotspots: how many hotspots
  ## prob.in: probability of success near hotspot
  ## prob.out: probability of success away from hotspot
  if (is.null(xlim)) {
    xlim = range(s[, 1])
  }
  if (is.null(ylim)) {
    ylim = range(s[, 2])
  }
  
  if ((length(xlim) != 2) | (length(ylim) != 2)) {
    stop("xlim and ylim must be length 2")
  }
  
  ns <- dim(s)[1]
  y  <- matrix(0, ns, 1)
  
  # generate the hotspots
  hotspots <- cbind(runif(nhotspots, xlim[1], xlim[2]), 
                    runif(nhotspots, ylim[1], ylim[2]))

  # figure out how many knots to select
  d  <- rdist(s, hotspots)
  
  for (i in 1:ns) {
    if (any(d[i, ] < bw)) {
      p <- prob.in
    } else {
      p <- prob.out
    }
    y[i, 1] <- rbinom(1, size = 1, prob = p)
  }
  
  results <- list(y = y, hotspots = hotspots)
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


dPS.Rcpp <- function(a, alpha, mid.points, bin.width, threads = 1) {
  if (is.null(dim(a))) {
    ns <- length(a)
    nt <- 1
    a <- matrix(a, ns, nt)  # turn it into a matrix
  } else {
    ns <- nrow(a)
    nt <- ncol(a)
  }

  results <- dPSCPP(a = a, alpha = alpha, mid_points = mid.points,
                    bin_width = bin.width, threads = 1)
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

dlognormal <- function(x, mu, sig) {
  dnorm(log(x), log(mu), sig, log = TRUE) - log(x)
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
    if (!is.matrix(a)) {
      a <- matrix(a, length(a), 1)
      lev <- as.vector(getLevelCPP(a = a, cuts = cuts))
    } else {
      lev <- getLevelCPP(a = a, cuts = cuts)
    }
  } else {
    lev <- sum(a > cuts) + 1
  }

  return(lev)
}

logdet <- function(X) {
  determinant(X)$modulus
}

trunc <- function(x, eps = 0.1) {
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
rTNorm <- function(mn, sd, lower = -Inf, upper = Inf, fudge = 0) {
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

dTNorm <- function(y, mn, sd, lower = -Inf, upper = Inf, fudge = 0,
                   log = TRUE) {
  lower.u <- pnorm(lower, mn, sd)
  upper.u <- pnorm(upper, mn, sd)

  ld <- dnorm(y, mn, sd, log = TRUE) - log(upper.u - lower.u)
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

fit.rarebinaryCPP <- function(xi.init, alpha.init, rho.init, beta.init,
                              xi.fix = FALSE, alpha.fix = FALSE,
                              rho.fix = FALSE, beta.fix = FALSE,
                              y, dw2, d, max.dist = NULL, cov,
                              method = "BFGS", hessian = TRUE,
                              alpha.min = 0, alpha.max = 1, threads = 1) {
  if (is.null(max.dist)) {
    max.dist <- max(d)
  }
  
  if (method == "Nelder-Mead") {
    hessian <- FALSE
  }

#   alpha.star.init <- log(alpha.init / (1 - alpha.init))
#   rho.star.init <- log(rho.init)

  # different combinations of parameters for fixing
  if (!xi.fix & !alpha.fix & !rho.fix & !beta.fix) {  # xi, alpha, rho, beta

    init.par <- c(xi.init, alpha.init, rho.init, beta.init)
    results  <- optim(init.par, pairwise.rarebinaryCPP.1, y = y,
                      dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                      alpha.min = alpha.min, alpha.max = alpha.max,
                      threads = threads, method = method, hessian = hessian)
    results$fixed <- "none"
    results$param.names <- c("xi", "alpha", "rho", "beta")

  } else if (xi.fix & !alpha.fix & !rho.fix & !beta.fix) {  # alpha, rho, and beta

    init.par <- c(alpha.init, rho.init, beta.init)
    results  <- optim(init.par, pairwise.rarebinaryCPP.2, y = y,
                      xi = xi.init,
                      dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                      alpha.min = alpha.min, alpha.max = alpha.max,
                      threads = threads, method = method, hessian = hessian)
    results$fixed <- "xi"
    results$param.names <- c("alpha", "rho", "beta")

  } else if (!xi.fix & alpha.fix & !rho.fix & !beta.fix) {  # xi, rho, and beta

    init.par <- c(xi.init, rho.init, beta.init)
    results  <- optim(init.par, pairwise.rarebinaryCPP.3, y = y,
                      alpha = alpha.init,
                      dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                      alpha.min = alpha.min, alpha.max = alpha.max,
                      threads = threads, method = method, hessian = hessian)
    results$fixed <- "alpha"
    results$param.names <- c("xi", "rho", "beta")

  } else if (!xi.fix & !alpha.fix & rho.fix & !beta.fix) {  # xi, alpha, and beta

    init.par <- c(xi.init, alpha.init, beta.init)
    results  <- optim(init.par, pairwise.rarebinaryCPP.4, y = y,
                      rho = rho.init,
                      dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                      alpha.min = alpha.min, alpha.max = alpha.max,
                      threads = threads, method = method, hessian = hessian)
    results$fixed <- "rho"
    results$param.names <- c("xi", "alpha", "beta")

  } else if (xi.fix & alpha.fix & !rho.fix & !beta.fix) {  # rho and beta

    init.par <- c(rho.init, beta.init)
    results  <- optim(init.par, pairwise.rarebinaryCPP.5, y = y,
                      xi = xi.init, alpha = alpha.init,
                      dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                      alpha.min = alpha.min, alpha.max = alpha.max,
                      threads = threads, method = method, hessian = hessian)
    results$fixed <- "xi and alpha"
    results$param.names <- c("rho", "beta")

  } else if (xi.fix & !alpha.fix & rho.fix & !beta.fix) {  # alpha and beta

    init.par <- c(alpha.init, beta.init)
    results  <- optim(init.par, pairwise.rarebinaryCPP.6, y = y,
                      xi = xi.init, rho = rho.init,
                      dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                      alpha.min = alpha.min, alpha.max = alpha.max,
                      threads = threads, method = method, hessian = hessian)
    results$fixed <- "xi and rho"
    results$param.names <- c("alpha", "beta")

  } else if (!xi.fix & alpha.fix & rho.fix & !beta.fix) {  # xi and beta

    init.par <- c(xi.init, beta.init)
    results  <- optim(init.par, pairwise.rarebinaryCPP.7, y = y,
                      alpha = alpha.init, rho = rho.init,
                      dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                      alpha.min = alpha.min, alpha.max = alpha.max,
                      threads = threads, method = method, hessian = hessian)
    results$fixed <- "alpha and rho"
    results$param.names <- c("xi", "beta")

  } else if (xi.fix & alpha.fix & rho.fix & !beta.fix) {  # beta only

    init.par <- c(beta.init)
    results  <- optim(init.par, pairwise.rarebinaryCPP.8, y = y,
                      xi = xi.init, alpha = alpha.init,
                      rho = rho.init,
                      dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                      alpha.min = alpha.min, alpha.max = alpha.max,
                      threads = threads, method = method, hessian = hessian)
    results$fixed <- "xi, alpha, and rho"
    results$param.names <- "beta"

  } else if (xi.fix & !alpha.fix & !rho.fix & beta.fix) {

    results.beta <- optim(par = beta.init, beta.hat, y = y,
                          cov = cov, xi = xi.init, method = "BFGS",
                          hessian = TRUE)
    beta.hat <- results.beta$par
    init.par <- c(alpha.init, rho.init)
    results <- optim(init.par, pairwise.rarebinaryCPP.9, y = y,
                     xi = xi.init, beta = beta.hat,
                     dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                     alpha.min = alpha.min, alpha.max = alpha.max,
                     threads = threads, method = method, hessian = hessian)
    results$fixed <- "xi"
    results$param.names <- c("alpha", "rho")
    results$beta <- beta.hat
    results$beta.cov <- solve(results.beta$hessian)

  } else if (xi.fix & !alpha.fix & rho.fix & beta.fix) {

    results.beta <- optim(par = beta.init, beta.hat, y = y,
                          cov = cov, xi = xi.init, method = "BFGS",
                          hessian = TRUE)
    beta.hat <- results.beta$par
    init.par <- alpha.init
    results <- optim(init.par, pairwise.rarebinaryCPP.10, y = y,
                     xi = xi.init, rho = rho.init, beta = beta.hat,
                     dw2 = dw2, d = d, max.dist = max.dist, cov = cov,
                     alpha.min = alpha.min, alpha.max = alpha.max,
                     threads = threads, method = method, hessian = hessian)
    results$fixed <- "xi, rho"
    results$param.names <- c("alpha")
    results$beta <- beta.hat
    results$beta.cov <- solve(results.beta$hessian)

  } else if (xi.fix & alpha.fix & !rho.fix & beta.fix) {

  }

  results$method <- method

  return(results)
}

getJoint2 <- function(kernel, alpha) {
  knot.contribute <- colSums(kernel)^alpha
  # print(knot.contribute)
  joint <- sum(knot.contribute)
  return(joint)
}

dic.spgev <- function(mcmcoutput, y, x, dw2, start = 1, end = NULL, thin = 1,
                      thresh = 0, update = NULL) {

  ns <- nrow(y)
  nt <- ncol(y)

  if (is.null(end)) {
    end <- length(mcmcoutput$xi)
  }

  niters <- length(start:end)
  nknots <- dim(mcmcoutput$a)[2]

  dbar <- -2 * mean(mcmcoutput$lly[start:end])

  if (is.null(dim(mcmcoutput$beta))) {
    p    <- 1
  } else {
    p    <- ncol(mcmcoutput$beta)
  }

  beta  <- matrix(mcmcoutput$beta[start:end, , drop = F], niters, p)
  xi    <- mcmcoutput$xi[start:end]
  a     <- mcmcoutput$a[start:end, , , drop = F]
  alpha <- mcmcoutput$alpha[start:end]
  rho   <- mcmcoutput$rho[start:end]

  betabar  <- apply(beta, 2, mean)
  xibar    <- mean(xi)
  abar     <- apply(a, c(2, 3), mean)
  alphabar <- mean(alpha)
  rhobar   <- mean(rho)

  x.beta     <- x %*% betabar
  z          <- getZ(xi = xibar, x.beta=x.beta, thresh=thresh)
  z.star     <- z^(1 / alphabar)
  w          <- stdW(makeW(dw2 = dw2, rho = rhobar))
  w.star     <- w^(1 / alphabar)
  theta.star <- getThetaStar(w.star = w.star, a = abar)
  prob       <- 1 - exp(-theta.star / z.star)
  dthetabar  <- -2 * sum(logLikeY(y, theta.star = theta.star,
                                    z.star = z.star))

  pd  <- dbar - dthetabar
  dic <- dbar + pd

  results <- list(dic = dic, pd = pd, dbar = dbar, dthetabar = dthetabar)
  return(results)
}


fit.rarebinaryInd <- function(init.par, y, cov) {
  results <- optim(init.par, ind.rarebinary, y = y, cov = cov,
                   method = "BFGS", hessian = TRUE)
}

# needs to loop over all pairs of sites and calculate the hessian
rarebinary.ij <- function(par, y, rho, dw2, cov, threads = 1) {
  npars <- length(par)
  alpha <- par[1]
  xi    <- par[2]
  beta  <- par[3:npars]

  W <- stdW(makeW(dw2, rho))
  ns <- length(y)

  H <- matrix(0, npars, npars)

  if (length(beta) == 1) {
    x.beta <- matrix(cov * beta, 2, 1)
  } else {
    x.beta <- cov %*% beta
  }

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)  # should be a matrix
  z      <- getZ(xi, x.beta)

  ll <- pairwiseCPPij(kernel = kernel[c(i, j), ], alpha = alpha,
                      z1 = z[1, 1], z2 = z[2, 1], y1 = y[1], y2 = y[2])

  return(ll)
}

# needs to loop over all pairs of sites and calculate the jacobian
jacobian.rarebinaryCPP <- function(par, y, rho, d, dw2, cov, threads = 1) {
  alpha <- par[1]
  xi    <- par[2]
  beta  <- par[3:npars]

  W <- stdW(makeW(dw2, rho))
  ns <- length(y)

  D <- matrix(0, 1, length(par))

  kernel <- exp((log(W) - log(as.vector(z))) / alpha)
  z      <- getZ(xi, x.beta)

  for (i in 1:(ns - 1)) {
    for (j in 1:ns) {
      if (d[i, j] < 0.4) {
        D <- D + jacobian(pairwiseCPPij(kernel = kernel[c(i, j), ],
                                        alpha = alpha,
                                        z1 = z[i], z2 = z[j],
                                        y1 = y[i], y2 = y[j]))
      }
    }
  }
}


##########
#### Old
##########
# get find the ll for y - returns ns x nt matrix for each site/day
# logLikeY2 <- function(y, theta.star, z.star, print = F) {
#   theta.z.star <- -theta.star / z.star
#
#   if (!is.null(dim(y))) {
#     ll.y <- matrix(-Inf, dim(y))
#   } else {
#     ll.y <- rep(-Inf, length(y))
#   }
#
#   # numerical stability issue. originally was using
#   # (1 - y) * P(Y = 0) + y * P(Y = 1)
#   # would return NaN because 0 * -Inf is not a number
#   these <- which(y == 1)
#   if (print) {
#     print(as.vector(theta.star))
#     print(unique(z.star))
#   }
#   ll.y[-these] <- theta.z.star[-these]
#   ll.y[these]  <- log(1 - exp(theta.z.star[these]))
#   # if (sum(is.nan(ll.y)) > 0) {
#   #   these <- which(is.nan(ll.y))
#   #   print(theta.z.star[these])
#   #   # print(1 - exp(theta.z.star[these]))
#   # }
#   return(ll.y)
# }