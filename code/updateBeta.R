updateBeta <- function(y, A, w, alpha, z, beta, x, cur.lly, acc, att, mh) {
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

  can.lly <- logLikeY(y=y, A=A, w=w, alpha=alpha, z=can.z)
  if (!is.na(R)) { if (log(runif(1)) < R) {
    beta    <- can.beta
    x.beta  <- can.x.beta
    cur.lly <- can.lly
    acc     <- acc + 1
  }}

  results <- list(beta=beta, x.beta=x.beta, cur.lly=cur.lly, att, acc)
}