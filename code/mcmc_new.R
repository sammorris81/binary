mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL,
                 xi.init=0.1, beta.init=0,
                 knots=NULL, thresh=NULL, logrho.init=1, alpha.init=0.5,
                 init.beta, init.alpha, init.range, init.bw, init.logs,
                 iterplot=FALSE, iters=50000, burn=10000, update=100, thin=1
    ) {
  library(fields)

  # initial setup
  p <- dim(x)[3]
  ns <- nrow(y)
  nt <- ncol(y)
  nknots <- dim(knots)[1]

  # get the initial set of weights for the sites and knots
  logrho <- logrho.init
  dw2    <- as.matrix(rdist(s, knots))^2  # dw2 is ns x nknots
  w      <- stdW(makeW(dw2, logrho))      # w is ns x nknots

  # get the initial set of random effects
  alpha <- alpha.init
  a     <- matrix(1, nknots, nt)  # random intensities
  theta <- matrix(1, ns, nt)      # sum_l a_l * w_l^(1/alpha)
  for (t in 1:nt) {
    theta[, t] <- getTheta(w, a[, t], alpha)
  }

  # get initial z
  xi     <- xi.init
  if (length(beta.init) == 1) {
    beta <- rep(beta.init, p)
  } else {
    beta   <- beta.init
  }
  x.beta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
    if (p > 1) {
      x.beta[, t] <- x[, t, ] %*% beta
    } else {
      x.beta[, t] <- x[, t] * beta
    }
  }
  z <- -(1 - xi * x.beta)^(-1 / xi)

  for (t in 1:nt) {
    for (k in 1:nknots) {
      cur.lp[k, t] <- dPS(a, alpha)
  }}

  cur.lly <- logLikeY(y=y, theta=theta, alpha=alpha, z=z)  # try to keep current

  # MH tuning parameters
  acc.beta  <- att.beta  <- mh.beta  <- 1
  acc.xi    <- att.xi    <- mh.xi    <- 1
  acc.a     <- att.a     <- mh.a     <- rep(1, 100)
  acc.alpha <- att.alpha <- mh.alpha <- 1
  acc.rho   <- att.rho   <- mh.rho   <- 1

  # storage
  beta.keepers  <- matrix(NA, nrow=iters, ncol=p)
  xi.keepers    <- rep(NA, iters)
  a.keepers     <- array(NA, dim=c(iters, nknots, nt))
  alpha.keepers <- rep(NA, iters)
  rho.keepers   <- rep(NA, iters)

  for (iter in 1:iters) { for (ttt in 1:nthin) {
    # update beta
    beta.update <- updateBeta(y=y, theta=theta, alpha=alpha, z=z,
                              beta=beta, xi=xi, x=x, cur.lly=cur.lly,
                              acc=acc.beta, att=att.beta, mh=mh.beta)  # TODO: Make function
    beta     <- beta.update$beta
    x.beta   <- beta.update$x.beta
    z        <- beta.update$z
    cur.lly  <- beta.update$cur.lly
    att.beta <- beta.update$att
    acc.beta <- beta.update$acc

    # update xi
    xi.update <- updateXi(y=y, A=A, w=w, alpha=alpha, z=z, x.beta=x.beta,
                          xi=xi, cur.lly, acc=acc.xi, att=att.xi, mh=mh.xi)
    xi      <- xi.update$xi
    z       <- xi.update$z
    cur.lly <- xi.update$cur.lly
    att.xi  <- xi.update$att
    acc.xi  <- xi.update$acc

    # update a
    a.update <- updateA(  # TODO: write function
                        )
    a       <- a.update$a
    theta   <- a.update$theta
    cur.lly <- a.update$cur.lly
    att.a   <- a.update$att
    acc.a   <- a.update$ac

    # update alpha
    alpha.update <- updateAlpha(  # TODO: write function
                                )

    alpha     <- alpha.update$alpha
    theta     <- alpha.update$theta
    cur.lly   <- alpha.update$cur.lly
    att.alpha <- alpha.update$att
    acc.alpha <- alpha.update$acc

    # update rho
    rho.update <- updateRho(  # TODO: write function
                            )
    rho     <- rho.update$rho
    theta   <- rho.update$theta
    cur.lly <- rho.update$cur.lly
    att.rho <- rho.update$att
    acc.rho <- rho.update$acc

  }}




}