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
  ll.y <- (1 - y) * theta.z.star + y * log(1 - exp(theta.z.star))
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
rRareBinaryInd <- function(x, beta, xi, thresh=0) {
  nt <- dim(x)[2]
  ns <- dim(x)[1]

  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  z <- matrix(rgev(n=ns * nt, 1, 1, 1), ns, nt)
  h <- x.beta + (z^xi - 1) / xi
  y <- ifelse(h > thresh, 1, 0)
  return(y)
}

# generate dependent rare binary data
rRareBinarySpat <- function(x, s, knots, beta, xi, alpha, rho, thresh=0,
                            dw2=NULL, a=NULL) {
  ns     <- nrow(s)
  nt     <- dim(x)[2]
  p      <- dim(x)[3]
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
  w      <- stdW(makeW(dw2=dw2, rho=rho))      # w is ns x nknots
  w.star <- w^(1 / alpha)

  # get random effects and theta.star
  if (is.null(a)) {
    a <- matrix(rPS(n=nknots * nt, alpha=alpha), nknots, nt)
  }
  theta.star <- getThetaStar(w.star=w.star, a=a)

  # get underlying latent variable
  u <- matrix(rgev(n=ns * nt, 1, alpha, alpha), ns, nt)
  z <- u * theta.star^alpha

  h <- x.beta + (z^xi - 1) / xi

  y <- ifelse(h > thresh, 1, 0)

  results <- list(y=y, a=a)
  return(results)
}


# update mh settings
mhUpdate <- function(acc, att, mh, nattempts=50, lower=0.8, higher=1.2) {
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
ld2 <- function(u, logs, alpha, shift=0, log=T) {

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

################################################################
# Arguments:
#   post.prob(iters, yp): posterior probability of exceeding
#   validate(np): validation data
#
# Returns:
#   scores(iters): posterior density of the brier scores
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
