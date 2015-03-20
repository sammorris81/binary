# update the beta term for the linear model
updateBeta <- function(y, theta.star, alpha, z, z.star, beta, beta.m, beta.s,
                       x.beta, xi, x, cur.lly, acc, att, mh, thresh=0) {
  # tried a block update for the beta update, but it doesn't do as well
  # as individual updates for each beta term separately.
  np  <- dim(x)[3]
  ns  <- nrow(y)
  nt  <- dim(x)[2]
  alpha.inv <- 1 / alpha

  for (p in 1:np) {
    att[p]      <- att[p] + 1
    can.beta    <- rnorm(1, beta[p], mh[p])
    # trying to save a little time
    can.x.beta  <- x.beta + x[, , p] * (can.beta - beta[p])
    can.z       <- getZ(xi=xi, x.beta=can.x.beta, thresh=thresh)
    can.z.star  <- can.z^(alpha.inv)

    # treat as independent at the moment
    can.lly <- logLikeY(y=y, theta.star=theta.star, z.star=can.z.star)

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
  can.z      <- getZ(xi=can.xi, x.beta=x.beta, thresh=thresh)
  can.z.star <- can.z^(alpha.inv)

  can.lly <- logLikeY(y=y, theta.star=theta.star, z.star=can.z.star)

  R <- sum(can.lly - cur.lly) +
       dnorm(can.xi, xi.m, xi.s, log=TRUE) -
       dnorm(xi, xi.m, xi.s, log=TRUE)

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

# updateBetaXi <- function(y, theta.star, alpha, z, beta, beta.m, beta.s, x,
#                          xi, xi.m, xi.s, cur.lly, acc.beta, att.beta, mh.beta,
#                          acc.xi, att.xi, mh.xi, thresh=0) {
#   np  <- dim(x)[3]
#   ns  <- nrow(y)
#   nt  <- dim(x)[2]
#   att.beta <- att.beta + 1
#   att.xi   <- att.xi + 1
#   alpha.inv <- 1 / alpha
#
#   can.beta <- rnorm(np, beta, mh.beta)
#   can.x.beta <- matrix(NA, ns, nt)
#   for (t in 1:nt) {
#     if (np > 1) {
#       can.x.beta[, t] <- x[, t, ] %*% can.beta
#     } else {
#       can.x.beta[, t] <- x[, t] * can.beta
#     }
#   }
#
#   can.xi <- rnorm(1, xi, mh.xi)
#   can.z <- getZ(xi=can.xi, x.beta=can.x.beta, thresh=thresh)
#   can.z.star <- can.z^alpha.inv
#
#   can.lly <- logLikeY(y=y, theta.star=theta.star, z.star=can.z.star)
#
#   R <- sum(can.lly - cur.lly) +
#        sum(dnorm(can.beta, beta.m, beta.s, log=TRUE)) -
#        sum(dnorm(beta, beta.m, beta.s, log=TRUE))
#        dnorm(can.xi, xi.m, xi.s, log=TRUE) -
#        dnorm(xi, xi.m, xi.s, log=TRUE)
#
#   if (!is.na(R)) { if (log(runif(1)) < R) {
#     beta     <- can.beta
#     x.beta   <- can.x.beta
#     xi       <- can.xi
#     z        <- can.z
#     cur.lly  <- can.lly
#     acc.beta <- acc.beta + 1
#     acc.xi   <- acc.xi + 1
#   }}
#
#   results <- list(beta=beta, x.beta=x.beta,
#                   xi=xi, z=z, cur.lly=cur.lly,
#                   att.beta=att.beta, acc.beta=acc.beta,
#                   att.xi=att.xi, acc.xi=acc.xi)
#   return(results)
# }

# update the random effects for theta.star
updateA <- function(y, theta.star, a, alpha, cur.lly, cur.llps, z.star, w.star,
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
      can.theta.star <- theta.star[, t] + w.star[, k] * (can.a - cur.a)
      if (sum(can.theta.star <= 0) > 0) {  # numerical stability
        can.llps <- -Inf
        can.lly.t <- -Inf
      } else {
        can.llps       <- dPS.Rcpp(a=can.a, alpha=alpha,
                                   mid.points=mid.points, bin.width=bin.width)
        can.lly.t      <- logLikeY(y=y[, t], theta.star=can.theta.star,
                                   z.star=z.star[, t])
      }

      R <- sum(can.lly.t - cur.lly.t) +
           can.llps - cur.llps[k, t] +
           dlognormal(cur.a, can.a, mh[l2]) - # candidate sd changes
           dlognormal(can.a, cur.a, mh[l1])

      if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
        a[k, t]         <- can.a
        theta.star[, t] <- can.theta.star
        cur.lly.t       <- can.lly.t
        cur.llps[k, t]  <- can.llps
      }}
    }
  }

  cur.lly <- logLikeY(y=y, theta.star=theta.star, z.star=z.star)

  results <- list(a=a, theta.star=theta.star,
                  cur.lly=cur.lly, cur.llps=cur.llps)
  return(results)
}

# update the alpha term for theta.star
updateAlpha <- function(y, theta.star, a, alpha, cur.lly, cur.llps, z, z.star,
                        w, w.star, mid.points, bin.width, acc, att, mh) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  att <- att + 1

  # not exactly U(0, 1) - using truncation for numerical stability
  cur.alpha.star <- transform$probit(alpha, 0.000001, 0.999999)
  can.alpha.star <- rnorm(1, cur.alpha.star, mh)
  can.alpha      <- transform$inv.probit(can.alpha.star, 0.000001, 0.999999)
  can.alpha.inv  <- 1 / can.alpha
  can.w.star     <- w^can.alpha.inv
  can.z.star     <- z^can.alpha.inv
  can.theta.star <- getThetaStar(can.w.star, a)
  can.lly        <- logLikeY(y=y, theta.star=can.theta.star, z.star=can.z.star)
  can.llps       <- matrix(NA, nknots, nt)
  can.llps       <- dPS.Rcpp(a, can.alpha, mid.points, bin.width)
  # for(t in 1:nt) {
  #   for (k in 1:nknots) {
  #     can.llps[k, t] <- dPS(a[k, t], can.alpha, mid.points, bin.width)
  #   }
  # }

  R <- sum(can.lly - cur.lly) + sum(can.llps - cur.llps) +
       dnorm(can.alpha.star, log=TRUE) -
       dnorm(cur.alpha.star, log=TRUE)

  if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
    alpha      <- can.alpha
    w.star     <- can.w.star
    z.star     <- can.z.star
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
                      rho, rho.upper=Inf, acc, att, mh) {
  nt     <- ncol(y)
  nknots <- nrow(a)

  att <- att + 1
  # rho.star       <- log(rho)
  # can.rho.star   <- rnorm(1, rho.star, mh)
  # can.rho        <- exp(can.rho.star)
  rho.star       <- transform$probit(rho, lower=0, upper=rho.upper)
  can.rho.star   <- rnorm(1, rho.star, mh)
  can.rho        <- transform$inv.probit(can.rho.star, lower=0, upper=rho.upper)
  can.w          <- stdW(makeW(dw2=dw2, rho=can.rho))
  can.w.star     <- can.w^(1 / alpha)
  can.theta.star <- getThetaStar(w.star=can.w.star, a=a)
  can.lly <- logLikeY(y=y, theta.star=can.theta.star, z.star=z.star)

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
    w.star     <- can.w.star
    theta.star <- can.theta.star
    cur.lly    <- can.lly
    acc        <- acc + 1
  }}

  results <- list(rho=rho, w=w, w.star=w.star, theta.star=theta.star,
                  cur.lly=cur.lly, att=att, acc=acc)
  return(results)
}