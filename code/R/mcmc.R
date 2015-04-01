mcmc <- function(y, s, x, s.pred=NULL, x.pred=NULL,
                 beta.init=0, beta.m=0, beta.s=10, xi.init=0.1, xi.m=0, xi.s=1,
                 npts=100, knots=NULL, thresh=0,
                 beta.tune=0.01, xi.tune=0.1,
                 alpha.tune=0.1, rho.tune=0.1, A.tune=1,
                 beta.attempts=50, xi.attempts=50,
                 alpha.attempts=200, rho.attempts=200,
                 rho.init=1, rho.upper=Inf, alpha.init=0.5, a.init=1,
                 iterplot=FALSE, iters=50000, burn=10000, update=100, thin=1
    ) {
  library(fields)

  # initial setup
  p <- dim(x)[3]
  ns <- nrow(y)
  nt <- ncol(y)
  nknots <- dim(knots)[1]

  # predictions
  predictions <- !is.null(s.pred) & !is.null(x.pred)
  if (predictions) {
    np   <- nrow(s.pred)
    dw2p <- as.matrix(rdist(s.pred, knots))^2
    keepers.y.pred <- array(0, c((iters - burn), np, nt))
  } else {
    np <- 0
    keepers.y.pred <- NULL
  }

  # get the initial set of weights for the sites and knots
  alpha  <- alpha.init
  rho    <- rho.init
  dw2    <- as.matrix(rdist(s, knots))^2  # dw2 is ns x nknots
  w      <- stdW(makeW(dw2, rho))         # w is ns x nknots
  w.star <- w^(1 / alpha)

  # get the initial set of random effects
  if (length(a.init) > 1) {
    a <- a.init
  } else {
    a <- matrix(a.init, nknots, nt)
  }
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
  z.star <- z^(1 / alpha)

  u.beta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
  mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
  bin.width <- u.beta[-1] - u.beta[-(npts + 1)]

  # keep current likelihood values in mcmc for time savings
  cur.lly  <- logLikeY(y=y, theta.star=theta.star, z.star=z.star)
  cur.llps <- dPS.Rcpp(a=a, alpha=alpha, mid.points=mid.points,
                       bin.width=bin.width)

  # MH tuning parameters
  if (length(beta.tune) == 1) {
    acc.beta <- att.beta <- mh.beta <- rep(beta.tune, p)
  } else {
    acc.beta <- att.beta <- mh.beta <- beta.tune
  }
  acc.xi    <- att.xi    <- mh.xi    <- xi.tune
  cuts      <- exp(c(-1, 0, 1, 2, 5, 10))
  mh.a      <- rep(A.tune, 100)
  acc.a     <- att.a     <- 0 * mh.a
  acc.alpha <- att.alpha <- mh.alpha <- alpha.tune
  acc.rho   <- att.rho   <- mh.rho   <- alpha.tune

  # storage
  keepers.beta  <- matrix(NA, nrow=iters, ncol=p)
  keepers.xi    <- rep(NA, iters)
  keepers.a     <- array(NA, dim=c(iters, nknots, nt))
  keepers.alpha <- rep(NA, iters)
  keepers.rho   <- rep(NA, iters)

  for (iter in 1:iters) { for (ttt in 1:thin) {
    # update beta and xi
    # update beta
    beta.update <- updateBeta(y=y, theta.star=theta.star, alpha=alpha,
                              z=z, z.star=z.star, beta=beta,
                              beta.m=beta.m, beta.s=beta.s, x.beta=x.beta,
                              xi=xi, x=x, cur.lly=cur.lly, acc=acc.beta,
                              att=att.beta, mh=mh.beta, thresh=thresh)
    beta     <- beta.update$beta
    x.beta   <- beta.update$x.beta
    z        <- beta.update$z
    z.star   <- beta.update$z.star
    cur.lly  <- beta.update$cur.lly
    att.beta <- beta.update$att
    acc.beta <- beta.update$acc

    if (iter < burn / 2) {
      mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta,
                            nattempts=beta.attempts)
      acc.beta  <- mh.update$acc
      att.beta  <- mh.update$att
      mh.beta   <- mh.update$mh
    }

    # update xi
    xi.update <- updateXi(y=y, theta.star=theta.star, alpha=alpha, z=z,
                          z.star=z.star, x.beta=x.beta, xi=xi, xi.m=xi.m,
                          xi.s=xi.s, cur.lly=cur.lly,
                          acc=acc.xi, att=att.xi, mh=mh.xi, thresh=thresh)
    xi      <- xi.update$xi
    z       <- xi.update$z
    z.star  <- xi.update$z.star
    cur.lly <- xi.update$cur.lly
    att.xi  <- xi.update$att
    acc.xi  <- xi.update$acc

    if (iter < burn / 2) {
      mh.update <- mhUpdate(acc=acc.xi, att=att.xi, mh=mh.xi,
                            nattempts=xi.attempts)
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

    if (iter < burn / 2) {
      mh.update <- mhUpdate(acc=acc.alpha, att=att.alpha, mh=mh.alpha,
                            nattempts=alpha.attempts)
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

    if (iter < burn / 2) {
      mh.update <- mhUpdate(acc=acc.rho, att=att.rho, mh=mh.rho,
                            nattempts=rho.attempts)
      acc.rho   <- mh.update$acc
      att.rho   <- mh.update$att
      mh.rho    <- mh.update$mh
    }

    }  # end thin

    if ((iter > burn) & predictions) {
      yp <- rRareBinarySpat(x=x.pred, s=s.pred, knots=knots, beta=beta, xi=xi,
                            alpha=alpha, rho=rho, thresh=thresh, dw2=dw2p)
      keepers.y.pred[(iter - burn), , ] <- yp$y
    }

    # storage
    keepers.beta[iter, ] <- beta
    keepers.xi[iter]       <- xi
    keepers.a[iter, , ]    <- a
    keepers.alpha[iter]    <- alpha
    keepers.rho[iter]      <- rho

    if (iter %% update == 0) {
      acc.rate.beta  <- round(acc.beta / att.beta, 3)
      acc.rate.xi    <- round(acc.xi / att.xi, 3)
      acc.rate.alpha <- round(acc.alpha / att.alpha, 3)
      acc.rate.rho   <- round(acc.rho / att.rho, 3)

      if (iterplot) {
        if (iter < burn) {
          begin <- max(1, (iter - 2000))
        } else {
          begin <- burn
        }

        par(mfrow=c(3, 6))
        plot(keepers.beta[begin:iter, 1], type="l", main=bquote(beta[0]),
             xlab=round(mh.beta[1], 4), ylab=acc.rate.beta[1])
        if (p > 1) {
          plot(keepers.beta[begin:iter, 2], type="l", main=bquote(beta[1]),
               xlab=round(mh.beta[2], 4), ylab=acc.rate.beta[2])
        }
        if (p > 2) {
          plot(keepers.beta[begin:iter, 3], type="l", main=bquote(beta[2]),
               xlab=round(mh.beta[3], 4), ylab=acc.rate.beta[3])
        }
        plot(keepers.xi[begin:iter], type="l", main=bquote(xi),
             xlab=round(mh.xi, 4), ylab=acc.rate.xi)
        plot(keepers.alpha[begin:iter], type="l", main=bquote(alpha),
             xlab=round(mh.alpha, 4), ylab=acc.rate.alpha)
        plot(keepers.rho[begin:iter], type="l", main=bquote(rho),
             xlab=round(mh.rho, 4), ylab=acc.rate.rho)


        for (i in 1:min(nt, 2)) {
          for (j in 1:6) {
            plot(keepers.a[begin:iter, j, i], type="l",
                 main=paste("knot ", j, ", day ", i, sep=""),
                 xlab="", ylab="")
          }
        }
      }
      cat("\t Iter", iter, "\n")
    }

  } # end iter

  return.iters <- (burn + 1):iters
  results <- list(beta=keepers.beta[return.iters, ],
                  xi=keepers.xi[return.iters],
                  a=keepers.a[return.iters, , ],
                  alpha=keepers.alpha[return.iters],
                  rho=keepers.rho[return.iters],
                  y.pred=keepers.y.pred)

  return(results)
}  # end mcmc
