# update the beta term for the linear model
updateBeta <- function(y, theta, alpha, z, beta, beta.m, beta.s, xi, x, cur.lly,
                       acc, att, mh) {
  # at the moment coding a block update for all parameters
  np  <- dim(x)[3]
  nt  <- dim(x)[2]
  att <- att + 1

  can.beta <- rnorm(np, beta, mh)
  can.x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    if (np > 1) {
      can.x.beta[, t] <- x[, t, ] %*% can.beta
    } else {
      can.x.beta[, t] <- x[, t] * can.beta
    }
  }
  can.z <- getZ(xi=xi, x.beta=can.x.beta)

  can.lly <- logLikeY(y=y, theta=theta, alpha=alpha, z=can.z)

  R <- sum(can.lly - cur.lly) +
       sum(dnorm(can.beta, beta.m, beta.s, log=TRUE)) -
       sum(dnorm(beta, beta.m, beta.s, log=TRUE))

  if (!is.na(R)) { if (log(runif(1)) < R) {
    beta    <- can.beta
    x.beta  <- can.x.beta
    z       <- can.z
    cur.lly <- can.lly
    acc     <- acc + 1
  }}

  # for (p in 1:np) {
  #   att[p] <- att[p] + 1
  #   can.beta[p] <- rnorm(1, beta[p], mh[p])
  #   can.x.beta  <- matrix(NA, ns, nt)
  #   for (t in 1:nt) {
  #     if (np > 1) {
  #       can.x.beta[, t] <- x[, t, ] %*% can.beta
  #     } else {
  #       can.x.beta[, t] <- x[, t] * can.beta
  #     }
  #   }
  #   can.z <- (1 + xi * can.x.beta)^(1 / xi)

  #   # treat as independent at the moment
  #   cur.lly <- y * (-1 / z) + (1 - y) * log(1 - exp(-1 / z))
  #   can.lly <- y * (-1 / can.z) + (1 - y) * log(1 - exp(-1 / can.z))
  #   # can.lly <- logLikeY(y=y, theta=theta, alpha=alpha, z=can.z)

  #   R <- sum(can.lly - cur.lly) +
  #        dnorm(can.beta[p], log=TRUE) - dnorm(beta[p], log=TRUE)

  #   if (!is.na(R)) { if (log(runif(1)) < R) {
  #     beta[p] <- can.beta[p]
  #     x.beta  <- can.x.beta
  #     z       <- can.z
  #     cur.lly <- can.lly
  #     acc[p]  <- acc[p] + 1
  #   }}
  # }

  results <- list(beta=beta, x.beta=x.beta, z=z, cur.lly=cur.lly,
                  att=att, acc=acc)
  return(results)
}

# update the xi term for the linear model
updateXi <- function(y, theta, alpha, z, x.beta, xi, xi.m, xi.s, cur.lly,
                     acc, att, mh) {
  nt <- dim(x)[2]
  att <- att + 1

  can.xi <- rnorm(1, xi, mh)
  can.z <- getZ(xi=can.xi, x.beta=x.beta)

  can.lly <- logLikeY(y=y, theta=theta, alpha=alpha, z=can.z)

  R <- sum(can.lly - cur.lly) +
       dnorm(can.xi, xi.m, xi.s, log=TRUE) -
       dnorm(xi, xi.m, xi.s, log=TRUE)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    xi      <- can.xi
    z       <- can.z
    cur.lly <- can.lly
    acc     <- acc + 1
  }}

  results <- list(xi=xi, z=z, cur.lly=cur.lly, att=att, acc=acc)
}

# update the random effects for theta
updateA <- function(y, theta, a, alpha, cur.lly, cur.llps, z, w,
                    mh, cuts) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  for (t in 1:nt) {
    for (k in 1:nknots) {
      l1        <- get.level(a[k, t], cuts)  # keeps sd reasonable
      can.a     <- exp(rnorm(1, log(a[k, t]), mh[l1]))
      l2        <- get.level(can.a, cuts)
      www       <- w[, k]^(1 / alpha)
      can.theta <- theta[, t] + www * (can.a - a[k, t])
      can.llps  <- dPS(can.a, alpha)
      can.lly   <- logLikeY(y=y[, t], theta=theta[, t], alpha=alpha, z=z[, t])

      R <- sum(can.lly - cur.lly[, t]) +
           can.llps - cur.llps[k, t] +
           dlognormal(a[k, t], cana, mh[l2]) - # candidate sd changes
           dlognormal(cana, a[k, t], mh[l1])

      if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
        a[t, k]        <- can.a
        theta[, t]     <- can.theta
        cur.lly[, t]   <- can.lly[, t]
        cur.llps[k, t] <- can.llps
      }}
    }
  }

  can.lly <- logLikeY(y=y, theta=theta, alpha=alpha, z=z)

  results <- list(a=a, theta=theta, cur.lly=cur.lly, cur.llps,
                  att=att, acc=acc)
}

# update the alpha term for theta
updateAlpha <- function(y, theta, a, alpha, cur.lly, cur.llps, z, w,
                        acc, att, mh) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  att <- att + 1
  cur.alpha.star <- transform$probit(alpha, 0, 1)
  can.alpha.star <- rnorm(1, cur.alpha.star, mh)
  can.alpha      <- transform$inv.probit(can.alpha.star, 0, 1)
  can.theta      <- getTheta(w, a, can.alpha)

  can.lly <- logLikeY(y=y, theta=can.theta, alpha=can.alpha, z=z)
  for(t in 1:nt) {
    for (k in 1:knots) {
      can.llps[k, t] <- dPS(a, can.alpha)
    }
  }

  R <- sum(can.lly - cur.lly) + sum(can.llps - cur.llps) +
       dnorm(can.alpha.star, log=TRUE) - dnorm(cur.alpha.star, log=TRUE)

  if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
    alpha    <- can.alpha
    theta    <- can.theta
    cur.lly  <- can.lly
    cur.llps <- can.llps
    acc      <- acc + 1
  }}


  results <- list(alpha=alpha, theta=theta, cur.lly=cur.lly, att=att, acc=acc)
}

updateRho <- function() {

  results <- list(rho=rho, theta=theta, cur.lly=cur.lly, att=att, acc=acc)
}