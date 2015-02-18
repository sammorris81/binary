mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL,
                 knots=NULL, thresh=NULL,
                 init.beta, init.alpha, init.range, init.bw, init.logs
    ) {

  # initial setup
  p <- dim(x)[3]
  ns <- nrow(y)
  nt <- ncol(y)

  # A <-  # TODO
  # w <-  # TODO
  # alpha <-  # TODO
  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    if (p > 1) {
      x.beta[, t] <- x[, t, ] %*% beta
    } else {
      x.beta[, t] <- x[, t] * beta
    }
  }
  # xi <-  # TODO
  z <- -(1 - xi * x.beta)^(-1 / xi)

  cur.lly <- logLikeY(y=y, A=A, w=w, alpha=alpha, z=z)  # try to keep current
  # MH tuning parameters

  # storage

  for (iter in 1:iters) { for (ttt in 1:nthin) {
    # update beta
    beta.update <- updateBeta(y=y, A=A, w=w, alpha=alpha, z=z,
                              beta=beta, x=x, cur.lly=cur.lly,
                              acc=acc.beta, att=att.beta, mh=mh.beta)  # TODO: Make function
    beta     <- beta.update$beta
    x.beta   <- beta.update$x.beta
    cur.lly  <- beta.update$cur.lly
    att.beta <- beta.update$att
    acc.beta <- beta.update$acc

    for (t in 1:nt) {
      if (p > 1) {
        x.beta[, t] <- x[, t, ] %*% beta
      } else {
        x.beta[, t] <- x[, t] * beta
      }
    }
    z <- -(1 - xi * x.beta)^(-1 / xi)

    # update xi

    z <- -(1 - xi * x.beta)^(-1 / xi)

    # update A

    # update alpha

    # update w
  }}




}