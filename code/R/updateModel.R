# update the beta term for the linear model
updateBeta <- function(y, theta, alpha, z, beta, xi, x, cur.lly,
                       acc, att, mh) {
  # at the moment coding a block update for all parameters
  p <- dim(x)[3]
  nt <- dim(x)[2]
  att <- att + 1

  can.beta <- rnorm(p, beta, mh)
  for (t in 1:nt) {
    if (p > 1) {
      can.x.beta[, t] <- x[, t, ] %*% can.beta
    } else {
      can.x.beta[, t] <- x[, t] * can.beta
    }
  }
  can.z <- -(1 - xi * can.x.beta)^(-1 / xi)

  can.lly <- logLikeY(y=y, theta=theta, alpha=alpha, z=can.z)

  R <- sum(can.lly - cur.lly) +
       sum(dnorm(can.beta, log=TRUE)) - sum(dnorm(beta, log=TRUE))

  if (!is.na(R)) { if (log(runif(1)) < R) {
    beta    <- can.beta
    x.beta  <- can.x.beta
    z       <- can.z
    cur.lly <- can.lly
    acc     <- acc + 1
  }}

  results <- list(beta=beta, x.beta=x.beta, z=z, cur.lly=cur.lly,
                  att=att, acc=acc)
}

# update the xi term for the linear model
updateXi <- function(y, theta, alpha, z, x.beta, xi, cur.lly,
                     acc, att, mh) {
  nt <- dim(x)[2]
  att <- att + 1

  can.xi <- rnorm(1, xi, mh)
  can.z <- -(1 - can.xi * x.beta)^(-1 / can.xi)

  can.lly <- logLikeY(y=y, theta=theta, alpha=alpha, z=can.z)

  R <- sum(can.lly - cur.lly) +
       dnorm(can.xi, log=TRUE) - dnorm(xi, log=TRUE)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    xi      <- can.xi
    z       <- can.z
    cur.lly <- can.lly
    acc     <- acc + 1
  }}

  results <- list(xi=xi, z=z, cur.lly=cur.lly, att=att, acc=acc)
}

# update the random effects for theta
updateA <- function(y, theta, alpha, z, w, acc, att, mh) {

  for (t in 1:nt) {
    for (k in 1:nknots) {
      l1        <- get.level(a[k, t], cuts)  # keeps sd reasonable
      can.a     <- exp(rnorm(1, log(a[k, t]), mh[l1]))
      l2        <- get.level(can.a, cuts)
      www       <- w[, k]^(1 / alpha)
      can.theta <- theta[, t] + www * (can.a - a[k, t])
      can.lp    <- dPS(can.a, alpha)
    }
  }



  results <- list(a=a, theta=theta, cur.lly=cur.lly, att=att, acc=acc)
}

# update the alpha term for theta
updateAlpha <- function() {

  results <- list(alpha=alpha, theta=theta, cur.lly=cur.lly, att=att, acc=acc)
}

updateRho <- function() {

  results <- list(rho=rho, theta=theta, cur.lly=cur.lly, att=att, acc=acc)
}