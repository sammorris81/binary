# update the beta term for the linear model
updateBeta <- function(y, theta.star, alpha, z, beta, beta.m, beta.s, xi, x,
                       cur.lly, acc, att, mh) {
  # tried a block update for the beta update, but it doesn't do as well
  # as individual updates for each beta term separately.
  np  <- dim(x)[3]
  ns  <- nrow(y)
  nt  <- dim(x)[2]
  # att <- att + 1

  # can.beta <- rnorm(np, beta, mh)
  # can.x.beta <- matrix(NA, ns, nt)
  # for (t in 1:nt) {
  #   if (np > 1) {
  #     can.x.beta[, t] <- x[, t, ] %*% can.beta
  #   } else {
  #     can.x.beta[, t] <- x[, t] * can.beta
  #   }
  # }
  # can.z <- getZ(xi=xi, x.beta=can.x.beta)

  # can.lly <- logLikeY(y=y, theta.star=theta.star, alpha=alpha, z=can.z)

  # R <- sum(can.lly - cur.lly) +
  #      sum(dnorm(can.beta, beta.m, beta.s, log=TRUE)) -
  #      sum(dnorm(beta, beta.m, beta.s, log=TRUE))

  # if (!is.na(R)) { if (log(runif(1)) < R) {
  #   beta    <- can.beta
  #   x.beta  <- can.x.beta
  #   z       <- can.z
  #   cur.lly <- can.lly
  #   acc     <- acc + 1
  # }}

  for (p in 1:np) {
    att[p] <- att[p] + 1
    can.beta <- beta
    can.beta[p] <- rnorm(1, beta[p], mh[p])
    can.x.beta  <- x.beta
    # trying to save a little time
    can.x.beta <- can.x.beta + can.beta[p] * x[, , p] - beta[p] * x[, , p]
    # for (t in 1:nt) {
    #   if (np > 1) {
    #     can.x.beta[, t] <- x[, t, ] %*% can.beta
    #   } else {
    #     can.x.beta[, t] <- x[, t] * can.beta
    #   }
    # }
    can.z <- getZ(xi=xi, x.beta=can.x.beta)

    # treat as independent at the moment
    can.lly <- logLikeY(y=y, theta.star=theta.star, alpha=alpha, z=can.z)

    R <- sum(can.lly - cur.lly) +
         dnorm(can.beta[p], beta.m, beta.s, log=TRUE) -
         dnorm(beta[p], beta.m, beta.s, log=TRUE)

    if (!is.na(R)) { if (log(runif(1)) < R) {
      beta[p] <- can.beta[p]
      x.beta  <- can.x.beta
      z       <- can.z
      cur.lly <- can.lly
      acc[p]  <- acc[p] + 1
    }}
  }

  results <- list(beta=beta, x.beta=x.beta, z=z, cur.lly=cur.lly,
                  att=att, acc=acc)
  return(results)
}

# update the xi term for the linear model
updateXi <- function(y, theta.star, alpha, z, x.beta, xi, xi.m, xi.s, cur.lly,
                     acc, att, mh) {
  nt <- ncol(y)
  att <- att + 1

  can.xi <- rnorm(1, xi, mh)
  can.z <- getZ(xi=can.xi, x.beta=x.beta)

  can.lly <- logLikeY(y=y, theta.star=theta.star, alpha=alpha, z=can.z)

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
  return(results)
}

updateBetaXi <- function(y, theta.star, alpha, z, beta, beta.m, beta.s, x,
                         xi, xi.m, xi.s, cur.lly, acc.beta, att.beta, mh.beta,
                         acc.xi, att.xi, mh.xi) {
  np  <- dim(x)[3]
  ns  <- nrow(y)
  nt  <- dim(x)[2]
  att.beta <- att.beta + 1
  att.xi   <- att.xi + 1

  can.beta <- rnorm(np, beta, mh.beta)
  can.x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    if (np > 1) {
      can.x.beta[, t] <- x[, t, ] %*% can.beta
    } else {
      can.x.beta[, t] <- x[, t] * can.beta
    }
  }

  can.xi <- rnorm(1, xi, mh.xi)
  can.z <- getZ(xi=can.xi, x.beta=can.x.beta)

  can.lly <- logLikeY(y=y, theta.star=theta.star, alpha=alpha, z=can.z)

  R <- sum(can.lly - cur.lly) +
       sum(dnorm(can.beta, beta.m, beta.s, log=TRUE)) -
       sum(dnorm(beta, beta.m, beta.s, log=TRUE))
       dnorm(can.xi, xi.m, xi.s, log=TRUE) -
       dnorm(xi, xi.m, xi.s, log=TRUE)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    beta     <- can.beta
    x.beta   <- can.x.beta
    xi       <- can.xi
    z        <- can.z
    cur.lly  <- can.lly
    acc.beta <- acc.beta + 1
    acc.xi   <- acc.xi + 1
  }}

  results <- list(beta=beta, x.beta=x.beta,
                  xi=xi, z=z, cur.lly=cur.lly,
                  att.beta=att.beta, acc.beta=acc.beta,
                  att.xi=att.xi, acc.xi=acc.xi)
  return(results)
}

# update the random effects for theta.star
updateA <- function(y, theta.star, a, alpha, cur.lly, cur.llps, z, w,
                    mid.points, bin.width, mh, cuts) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  for (t in 1:nt) {
    cur.lly.t <- cur.lly[, t]
    # cur.lly.t <- logLikeY(y=y[, t], theta.star=theta.star[, t], alpha=alpha, z=z[, t])
    for (k in 1:nknots) {
      l1             <- get.level(a[k, t], cuts)  # keeps sd reasonable
      can.a          <- exp(rnorm(1, log(a[k, t]), mh[l1]))
      l2             <- get.level(can.a, cuts)
      www            <- w[, k]^(1 / alpha)  # w is ns x nknots
      # can.theta.star only changes at a site when it's near the knot
      can.theta.star <- theta.star[, t] + www * (can.a - a[k, t])
      can.llps       <- dPS(a=can.a, alpha=alpha,
                            mid.points=mid.points, bin.width=bin.width)
      can.lly.t      <- logLikeY(y=y[, t], theta.star=can.theta.star,
                                 alpha=alpha, z=z[, t])

      R <- sum(can.lly.t - cur.lly.t) +
           can.llps - cur.llps[k, t] +
           dlognormal(a[k, t], can.a, mh[l2]) - # candidate sd changes
           dlognormal(can.a, a[k, t], mh[l1])

      if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
        a[k, t]         <- can.a
        theta.star[, t] <- can.theta.star
        cur.lly.t       <- can.lly.t
        cur.llps[k, t]  <- can.llps
      }}
    }
  }

  cur.lly <- logLikeY(y=y, theta.star=theta.star, alpha=alpha, z=z)

  results <- list(a=a, theta.star=theta.star,
                  cur.lly=cur.lly, cur.llps=cur.llps)
  return(results)
}

# update the alpha term for theta.star
updateAlpha <- function(y, theta.star, a, alpha, cur.lly, cur.llps, z, w,
                        mid.points, bin.width, acc, att, mh) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  att <- att + 1

  # not exactly U(0, 1) - using truncation for numerical stability
  cur.alpha.star <- transform$probit(alpha, 0.000001, 0.999999)
  can.alpha.star <- rnorm(1, cur.alpha.star, mh)
  can.alpha      <- transform$inv.probit(can.alpha.star, 0.000001, 0.999999)
  can.theta.star <- getThetaStar(w, a, can.alpha)
  can.lly        <- logLikeY(y=y, theta.star=can.theta.star, alpha=can.alpha,
                             z=z)
  can.llps       <- matrix(NA, nknots, nt)
  for(t in 1:nt) {
    for (k in 1:nknots) {
      can.llps[k, t] <- dPS(a[k, t], can.alpha, mid.points, bin.width)
    }
  }

  R <- sum(can.lly - cur.lly) + sum(can.llps - cur.llps) +
       dnorm(can.alpha.star, 0, 5, log=TRUE) -
       dnorm(cur.alpha.star, 0, 5, log=TRUE)

  if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
    alpha      <- can.alpha
    theta.star <- can.theta.star
    cur.lly    <- can.lly
    cur.llps   <- can.llps
    acc        <- acc + 1
  }}

  results <- list(alpha=alpha, theta.star=theta.star, cur.lly=cur.lly,
                  cur.llps=cur.llps, att=att, acc=acc)
  return(results)
}

updateRho <- function(y, theta.star, a, alpha, cur.lly, z, w, dw2, rho,
                      rho.upper=Inf, acc, att, mh) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  att <- att + 1
  # rho.star       <- log(rho)
  # can.rho.star   <- rnorm(1, rho.star, mh)
  # can.rho        <- exp(can.rho.star)
  rho.star <- transform$probit(rho, lower=0, upper=rho.upper)
  can.rho.star <- rnorm(1, rho.star, mh)
  can.rho  <- transform$inv.probit(can.rho.star, lower=0, upper=rho.upper)
  can.w          <- stdW(makeW(dw2=dw2, rho=can.rho))
  can.theta.star <- getThetaStar(w=can.w, a=a, alpha=alpha)
  can.lly <- logLikeY(y=y, theta.star=can.theta.star, alpha=alpha, z=z)

  logrho.m <- -1
  logrho.s <- 0.1

  R <- sum(can.lly - cur.lly) +
       dnorm(can.rho.star, log=TRUE) - dnorm(rho.star, log=TRUE)

       # dnorm(can.rho.star, logrho.m, logrho.s, log=TRUE) -
       # dnorm(rho.star, logrho.m, logrho.s, log=TRUE)
       # dgamma(can.rho, 1, 1, log=T) -
       # dgamma(rho, 1, 1, log=T)

  if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
    rho        <- can.rho
    w          <- can.w
    theta.star <- can.theta.star
    cur.lly    <- can.lly
    acc        <- acc + 1
  }}

  results <- list(rho=rho, w=w, theta.star=theta.star, cur.lly=cur.lly,
                  att=att, acc=acc)
  return(results)
}