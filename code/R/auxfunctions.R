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

getZ <- function(xi, x.beta) {
  if (xi != 0) {
    z <- (1 + xi * x.beta)^(1 / xi)
  } else {
    z <- exp(x.beta)
  }
}

# TODO: Is this theta or theta^alpha. I think it's theta^alpha
# theta sum_l=1^L a_l * w_l^(1 / alpha)
getTheta <- function(w, a, alpha) {
  # theta is nxnF
  # s is nFxnt
  # alpha in (0,1)
  w.star <- w^(1 / alpha)
  if (length(a) == 1) {xxx <- w.star * a}
  if (length(a) > 1) {xxx <- w.star %*% a}
  return(xxx)
}

# get the kernel weighting
makeW <- function(dw2, logrho) {
  # rho.star = log(rho)
  rho2 <- exp(logrho)^2
  w <- exp(-0.5 * dw2 / rho2)
  return(w)
}

# standardize the kernel weights
stdW <- function(x, single=FALSE) {
  if (single) {x <- x / sum(x)}
  if (!single) {x <- sweep(x, 1, rowSums(x), "/")}
  return(x)
}


# get find the ll for y - returns ns x nt matrix for each site/day
logLikeY <- function(y, theta, alpha, z) {
  z.star <- -(theta / z)^(1 / alpha)
  ll.y <- (1 - y) * z.star + y * log(1 - exp(z.star))
  return(ll.y)
}

################################################################################
#### positive stable density functions
################################################################################
dPS <- function(a, alpha, npts=100) {
  Ubeta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
  MidPoints <- (Ubeta[-1] + Ubeta[-(npts + 1)]) / 2
  BinWidth <- Ubeta[-1] - Ubeta[-(npts + 1)]
  l <- -Inf

  if (a > 0) {
    l <- log(sum(BinWidth * ld(MidPoints, a, alpha)))
  }

  return(l)
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
rRareBinaryInd <- function(x, beta, xi) {
  nt <- dim(x)[2]
  ns <- dim(x)[1]

  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  z    <- (1 + xi * x.beta)^(1 / xi)
  prob <- 1 - exp(-1 / z)  # we need P(Y = 1) for rbinom
  y    <- matrix(rbinom(n=ns * nt, size=1, prob=prob), nrow=ns, ncol=nt)
  return(y)
}

# generate dependent rare binary data
rRareBinarySpat <- function(x, s, knots, beta, xi, alpha, rho) {
  y      <- matrix(NA, ns, nt)
  p      <- dim(x)[3]
  nknots <- nrow(knots)

  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    if (p > 1) {
      x.beta[, t] <- x[, t, ] %*% beta
    } else {
      x.beta[, t] <- x[, t] * beta
    }
  }
  z <- getZ(xi=xi, x.beta=x.beta)

  # get weights
  dw2 <- as.matrix(rdist(s, knots))^2  # dw2 is ns x nknots
  w <- stdW(makeW(dw2, log(rho)))      # w is ns x nknots

  # get random effects and theta
  a     <- matrix(rPS(n=nknots * nt, alpha=alpha), nknots, nt)
  theta <- matrix(1, ns, nt)
  for (t in 1:nt) {
    theta[, t] <- getTheta(w, a[, t], alpha)
  }

  prob <- 1 - exp(-(theta / z)^(1 / alpha))  # we need P(Y = 1) for rbinom

  for (t in 1:nt) {
    y[, t] <- rbinom(n=ns, size=1, prob=prob)
  }

  results <- list(y=y, a=a)
  return(results)
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


ECkern <- function(h, alpha, gamma, Lmax=50) {
  dw2 <- rdist(c(0, h), seq(-Lmax, Lmax, 1))
  W <- fac2FAC(make.fac(dw2, gamma))^(1 / alpha)
  for (j in 1:length(h)) {
    h[j]<-sum((W[1, ] + W[j + 1, ])^alpha)
  }

  return(h)
}



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



logd <- function(theta, v) {
  sum(log(theta) - theta * v)
}

#### Used to set the standard deviation for the candidate distribution
#### for the A terms in the random effect. When log(A) is large means
#### sd is smaller, and log(A) small means sd is larger.
get.level <- function(a, cuts) {
  lev <- a * 0 + 1
  for (j in 1:length(cuts)) {
    lev <- ifelse(a > cuts[j], j + 1, lev)
  }
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

