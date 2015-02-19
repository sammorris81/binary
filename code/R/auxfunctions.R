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
makeW <- function(dw2, rho.star) {
  # rho.star is in (-Inf, Inf)
  rho2 <- exp(rho.star)^2
  fac <- exp(-0.5 * dw2 / rho2)
  return(fac)
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
  ll.y <- (1 - y) * z.star + y * log(1 - z.star)
  return(ll.y)
}

npts <- 50
Ubeta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
MidPoints <- (Ubeta[-1] + Ubeta[-(npts + 1)]) / 2
BinWidth <- Ubeta[-1] - Ubeta[-(npts + 1)]

dPS <- function(a, alpha, npts=100) {
  l <- -Inf
  if (a > 0) {
    l <- log(sum(BinWidth * ld(MidPoints, a, alpha)))
  }

  return(l)
}




ECkern <- function(h, alpha, gamma, Lmax=50) {
  dw2 <- rdist(c(0, h), seq(-Lmax, Lmax, 1))
  W <- fac2FAC(make.fac(dw2, gamma))^(1 / alpha)
  for (j in 1:length(h)) {
    h[j]<-sum((W[1, ] + W[j + 1, ])^alpha)
  }

  return(h)
}



ld <- function(u, A, alpha) {
  psi <- pi * u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)
  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * log(A) +
          log(c) - c * (1 / A^(alpha / (1 - alpha)))

  return(exp(logd))
}


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

# generate rare binary data
# for now working with independent to test functions for beta and xi
rBinaryRareInd <- function(x, beta, xi) {
  nt <- dim(x)[2]
  ns <- dim(x)[1]

  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    x.beta[, t] <- x[, t, ] %*% beta
  }

  z <- (1 + xi * x.beta)^(1 / xi)
  p <- 1 - exp(-1 / z)  # we need P(Y = 1) for rbinom
  y <- rbinom(n=ns * nt, size=1, prob=p)
  return(y)
}

# generate dependent rare binary data
rRareBinary <- function(x, beta, xi, alpha, rho) {

}
