mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL,
                 beta.init=0, beta.m=0, beta.s=10,
                 xi.init=0.1, xi.m=0, xi.s=1, npts=100, knots=NULL,
                 thresh=NULL, rho.init=1, rho.upper=Inf, alpha.init=0.5,
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
  rho    <- rho.init
  dw2    <- as.matrix(rdist(s, knots))^2  # dw2 is ns x nknots
  w      <- stdW(makeW(dw2, rho))         # w is ns x nknots
  w.star <- w^(1 / alpha)

  # get the initial set of random effects
  alpha      <- alpha.init
  a          <- matrix(1, nknots, nt)  # random intensities
  theta.star <- getThetaStar(w.star, a)  # sum_l a_l * w_l^(1/alpha)

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
  z <- getZ(xi=xi, x.beta=x.beta)

  u.beta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
  mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
  bin.width <- u.beta[-1] - u.beta[-(npts + 1)]

  for (t in 1:nt) {
    for (k in 1:nknots) {
      cur.llps[k, t] <- dPS(a=a, alpha=alpha,
                            mid.points=mid.points, bin.width=bin.width)
  }}

  # keep current in mcmc for time savings
  cur.lly <- logLikeY(y=y, theta.star=theta.star, alpha=alpha, z=z)

  # MH tuning parameters
  acc.beta  <- att.beta  <- mh.beta  <- 1
  acc.xi    <- att.xi    <- mh.xi    <- 1
  cuts      <- exp(c(-1, 0, 1, 2, 5, 10))
  mh.a      <- rep(1, 100)
  acc.a     <- att.a     <- 0 * mh.a
  acc.alpha <- att.alpha <- mh.alpha <- 1
  acc.rho   <- att.rho   <- mh.rho   <- 1

  # storage
  keepers.beta  <- matrix(NA, nrow=iters, ncol=p)
  keepers.xi    <- rep(NA, iters)
  keepers.a     <- array(NA, dim=c(iters, nknots, nt))
  keepers.alpha <- rep(NA, iters)
  keepers.rho   <- rep(NA, iters)

  for (iter in 1:iters) { for (ttt in 1:nthin) {
    # update beta and xi
    # update beta
    beta.update <- updateBeta(y=y, theta.star=theta.star, alpha=alpha, 
                              z=z, z.star=z.star, 
                              beta=beta, beta.m=beta.m, beta.s=beta.s,
                              xi=xi, x=x, cur.lly=cur.lly,
                              acc=acc.beta, att=att.beta, mh=mh.beta)
    beta     <- beta.update$beta
    x.beta   <- beta.update$x.beta
    z        <- beta.update$z
    z.star   <- beta.update$z.star
    cur.lly  <- beta.update$cur.lly
    att.beta <- beta.update$att
    acc.beta <- beta.update$acc

    if (i < burn / 2) {
      mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta)
      acc.beta  <- mh.update$acc
      att.beta  <- mh.update$att
      mh.beta   <- mh.update$mh
    }

    # update xi
    xi.update <- updateXi(y=y, theta.star=theta.star, alpha=alpha, z=z,
                          z.star=z.star, x.beta=x.beta, 
                          xi=xi, xi.m=xi.m, xi.s=xi.s,
                          cur.lly, acc=acc.xi, att=att.xi, mh=mh.xi)
    xi      <- xi.update$xi
    z       <- xi.update$z
    z.star  <- xi.update$z.star
    cur.lly <- xi.update$cur.lly
    att.xi  <- xi.update$att
    acc.xi  <- xi.update$acc

    if (i < burn / 2) {
      mh.update <- mhUpdate(acc=acc.xi, att=att.xi, mh=mh.xi)
      acc.xi  <- mh.update$acc
      att.xi  <- mh.update$att
      mh.xi   <- mh.update$mh
    }

    # update a - NOTE: does not use acc, att, and mh like usual
    old.a    <- a
    a.update <- updateA(y=y, theta.star=theta.star, a=a, alpha=alpha,
                        cur.lly=cur.lly, cur.llps=cur.llps, z.star=z.star, 
                        w.star=w.star, mid.points=mid.points, 
                        bin.width=bin.width, mh=mh.a, cuts=cuts)
    a          <- a.update$a
    theta.star <- a.update$theta.star
    cur.lly    <- a.update$cur.lly
    cur.llps   <- a.update$cur.llps

    # not going to worry about tuning the a terms at the moment
    # level <- get.level(old.a, cuts)
    # for (j in 1:length(mh.a)) {
    #   acc.a[j] <- acc.a[j] + sum(old.a[level == j] != a[level == j])
    #   att.a[j] <- att.a[j] + sum(level == j)
    #   if ((i < burn / 2) & (att.a[j] > 100)) {
    #     if (acc.a[j] / att.a[j] < 0.3) { mh.a[j] <- mh.a[j] * 0.9 }
    #     if (acc.a[j] / att.a[j] > 0.6) { mh.a[j] <- mh.a[j] * 1.1 }
    #     acc.a[j] <- att.a[j] <- 0
    #   }
    # }

    # update alpha
    alpha.update <- updateAlpha(y=y, theta.star=theta.star, a=a, alpha=alpha,
                                cur.lly=cur.lly, cur.llps=cur.llps, 
                                z=z, z.star=z.star, w=w, w.star=w.star, 
                                mid.points=mid.points, bin.width=bin.width,
                                acc=acc.alpha, att=att.alpha, mh=mh.alpha)

    alpha     <- alpha.update$alpha
    w.star    <- alpha.update$w.star
    z.star    <- alpha.update$z.star
    theta     <- alpha.update$theta.star
    cur.lly   <- alpha.update$cur.lly
    cur.llps  <- alpha.update$cur.llps
    att.alpha <- alpha.update$att
    acc.alpha <- alpha.update$acc

    if (i < burn / 2) {
      mh.update <- mhUpdate(acc=acc.alpha, att=att.alpha, mh=mh.alpha)
      acc.alpha <- mh.update$acc
      att.alpha <- mh.update$att
      mh.alpha  <- mh.update$mh
    }

    # update rho
    rho.update <- updateRho(y=y, theta.star=theta.star, a=a, alpha=alpha,
                            cur.lly=cur.lly, z.star=z.star, w=w, w.star=w.star, 
                            dw2=dw2, rho=rho, rho.upper=rho.upper,
                            acc=acc.rho, att=att.rho, mh=mh.rho)
    rho        <- rho.update$rho
    w          <- rho.update$w
    w.star     <- rho.update$w.star
    theta.star <- rho.update$theta.star
    cur.lly    <- rho.update$cur.lly
    att.rho    <- rho.update$att
    acc.rho    <- rho.update$acc

    if (i < burn / 2) {
      mh.update <- mhUpdate(acc=acc.rho, att=att.rho, mh=mh.rho)
      acc.rho   <- mh.update$acc
      att.rho   <- mh.update$att
      mh.rho    <- mh.update$mh
    }

  }  # end thin

  # storage
  keepers.beta[iter, , ] <- beta
  keepers.xi[iter]       <- xi
  keepers.a[iter, , ]    <- a
  keepers.alpha[iter]    <- alpha
  keepers.rho[iter]      <- rho

  }

}