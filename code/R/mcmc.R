mcmc <- function(y, s, x, s.pred = NULL, x.pred = NULL,
                 beta.init = NULL, beta.m = 0, beta.s = 20,
                 xi.init = NULL, xi.m = 0, xi.s = 0.5,
                 npts = 100, knots = NULL, thresh = 0,
                 beta.tune = 0.01, xi.tune = 0.1,
                 alpha.tune = 0.1, alpha.m = 0.5, alpha.s = sqrt(1 / 12),
                 rho.tune = 0.1, logrho.m = -1, logrho.s = 2,
                 A.tune = 1, IDs = NULL,
                 beta.attempts = 50, xi.attempts = 50,
                 alpha.attempts = 200, rho.attempts = 200,
                 spatial = TRUE,
                 rho.init = 1, rho.upper = Inf, alpha.init = 0.5, a.init = 1,
                 beta.fix = FALSE, xi.fix = FALSE,   # debug
                 rho.fix = FALSE, alpha.fix = FALSE, # debug
                 xibeta.joint = FALSE,
                 iterplot=FALSE, iters=50000, burn=10000, update=100, thin=1
    ) {
  library(fields)

  # initial setup
  if (is.null(dim(y))){
    y <- matrix(y, length(y), 1)
  }

  ns <- nrow(y)
  nt <- ncol(y)

  x <- adjustX(x=x, y=y)
  np <- ncol(x)

  nknots <- dim(knots)[1]

  # predictions
  predictions <- !is.null(s.pred) & !is.null(x.pred)
  if (predictions) {
    npred <- nrow(s.pred)
    dw2p  <- as.matrix(rdist(s.pred, knots))^2
    keepers.y.pred <- array(0, c((iters - burn), npred, nt))
  } else {
    npred <- 0
    keepers.y.pred <- NULL
  }

  # get initial z
  if (is.null(beta.init)) {
    if (xi == 0) {
      beta.init <- thresh + log(-log(1 - mean(y)))
      beta <- rep(beta.init, np)
    } else {
      beta.init <- thresh + (1 - log(1 - mean(y)))^(-xi) / xi
      beta <- rep(beta.init, np)
    }
  } else {
    if (length(beta.init) == 1) {
      beta <- rep(beta.init, np)
    } else {
      beta   <- beta.init
    }
  }

  if (is.null(xi.init)) {
    xi <- 0
  } else {
    xi <- xi.init
  }

  x.beta <- getXBeta(x, ns, nt, beta)
  xt <- t(x)
  xtx.inv <- solve(xt %*% x)

  z <- getZ(xi=xi, x.beta=x.beta)

  # get the initial set of weights for the sites and knots
  rho    <- rho.init
  dw2    <- as.matrix(rdist(s, knots))^2  # dw2 is ns x nknots
  w      <- stdW(makeW(dw2, rho))         # w is ns x nknots
  if (length(a.init) > 1) {
    a <- a.init
  } else {
    a <- matrix(a.init, nknots, nt)
  }
  if (is.null(IDs)) {
    # if not specified, make the list of sites that are impacted by each knot
    IDs <- vector(mode = "list", length = nknots)
  }

  if (spatial) {
    alpha <- alpha.init
    u.beta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
    mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
    bin.width <- u.beta[-1] - u.beta[-(npts + 1)]
    cur.llps <- dPS.Rcpp(a=a, alpha=alpha, mid.points=mid.points,
                         bin.width=bin.width)
  } else {
    alpha <- 1
  }

  alpha.a <- (1 - alpha.m) * alpha.m^2 / alpha.s^2 - alpha.m
  alpha.b <- alpha.a * (1 / alpha.m - 1)
  wz.star <- getwzStarCPP(z = z, w = w, alpha = alpha)
  kernel  <- getKernelCPP(wz_star = wz.star, a = a)
#   wz.star <- getwzStar(z = z, w = w, alpha = alpha)
#   kernel  <- getKernel(wz.star = wz.star, a = a)

  # keep current likelihood values in mcmc for time savings
  cur.lly  <- logLikeY(y = y, kernel = kernel)

  # MH tuning parameters
  if (length(beta.tune) == 1) {
    acc.beta <- att.beta <- mh.beta <- rep(beta.tune, np)
  } else {
    acc.beta <- att.beta <- mh.beta <- beta.tune
  }
  acc.xi    <- att.xi    <- mh.xi    <- xi.tune
  cuts      <- exp(c(-1, 0, 1, 2, 5, 10))
  mh.a      <- rep(A.tune, 100)
  acc.a     <- att.a     <- 0 * mh.a
  acc.alpha <- att.alpha <- mh.alpha <- alpha.tune
  acc.rho   <- att.rho   <- mh.rho   <- rho.tune

  # storage
  keepers.beta  <- matrix(NA, nrow=iters, ncol=np)
  keepers.xi    <- rep(NA, iters)
  keepers.a     <- array(NA, dim=c(iters, nknots, nt))
  keepers.alpha <- rep(NA, iters)
  keepers.rho   <- rep(NA, iters)
  keepers.lly   <- rep(NA, iters)

  for (iter in 1:iters) { for (ttt in 1:thin) {

    if (xibeta.joint) {  # update beta and xi
      # we are actual sampling for p = P(Y = 0)
      xibeta.update <- updateXiBeta(y = y, alpha = alpha, z = z, w = w,
                                    wz.star = wz.star, beta = beta,
                                    kernel = kernel, a = a, x.beta = x.beta,
                                    xi = xi, x = x, xt = xt, xtx.inv = xtx.inv,
                                    xi.m = xi.m, xi.s = xi.s, cur.lly = cur.lly,
                                    xi.fix = xi.fix, beta.fix = beta.fix,
                                    acc.p = acc.beta, att.p = att.beta,
                                    mh.p = mh.beta,
                                    acc.xi = acc.xi, att.xi = att.xi,
                                    mh.xi = mh.xi, thresh = 0)

      beta     <- xibeta.update$beta
      x.beta   <- xibeta.update$x.beta
      xi       <- xibeta.update$xi
      z        <- xibeta.update$z
      wz.star  <- xibeta.update$wz.star
      kernel   <- xibeta.update$kernel
      cur.lly  <- xibeta.update$cur.lly
      att.beta <- xibeta.update$att.p
      acc.beta <- xibeta.update$acc.p
      att.xi   <- xibeta.update$att.xi
      acc.xi   <- xibeta.update$acc.xi

      if (iter < burn / 2) {
        mh.update <- mhUpdate(acc=acc.beta, att=att.beta, mh=mh.beta,
                              nattempts=beta.attempts)
        acc.beta  <- mh.update$acc
        att.beta  <- mh.update$att
        mh.beta   <- mh.update$mh

        mh.update <- mhUpdate(acc=acc.xi, att=att.xi, mh=mh.xi,
                              nattempts=xi.attempts)
        acc.xi  <- mh.update$acc
        att.xi  <- mh.update$att
        mh.xi   <- mh.update$mh
      }

    } else {  # update beta
      if (!beta.fix) {
        beta.update <- updateBeta(y = y, kernel = kernel, alpha = alpha,
                                  a = a, z = z, w = w, wz.star = wz.star,
                                  beta = beta, beta.m = beta.m, beta.s = beta.s,
                                  x.beta = x.beta, xi = xi, x = x,
                                  cur.lly = cur.lly, acc = acc.beta,
                                  att = att.beta, mh = mh.beta, thresh = thresh)
        beta     <- beta.update$beta
        x.beta   <- beta.update$x.beta
        z        <- beta.update$z
        wz.star  <- beta.update$wz.star
        kernel   <- beta.update$kernel
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
      }

      # update xi
      if (!xi.fix) {
        xi.update <- updateXi(y = y, kernel = kernel, alpha = alpha,
                              a = a, z = z, w = w, wz.star = wz.star,
                              x.beta = x.beta, xi = xi, xi.m = xi.m,
                              xi.s = xi.s, cur.lly = cur.lly, acc = acc.xi,
                              att = att.xi, mh = mh.xi, thresh = thresh)
        xi      <- xi.update$xi
        z       <- xi.update$z
        wz.star <- xi.update$wz.star
        kernel  <- xi.update$kernel
        cur.lly <- xi.update$cur.lly
        att.xi  <- xi.update$att
        acc.xi  <- xi.update$acc

        if (iter < burn / 2) {
          mh.update <- mhUpdate(acc = acc.xi, att = att.xi, mh = mh.xi,
                                nattempts = xi.attempts)
          acc.xi  <- mh.update$acc
          att.xi  <- mh.update$att
          mh.xi   <- mh.update$mh
        }
      }
    }

    if (spatial) {
      # update a - NOTE: does not use acc, att, and mh like usual
      old.a    <- a
      a.update <- updateA(y = y, kernel = kernel, a = a, alpha = alpha,
                          wz.star = wz.star, cur.lly = cur.lly,
                          cur.llps = cur.llps, mid.points = mid.points,
                          bin.width = bin.width, mh = mh.a, cuts = cuts)
      a        <- a.update$a
      kernel   <- a.update$kernel
      cur.lly  <- a.update$cur.lly
      cur.llps <- a.update$cur.llps

      # adjust the candidate standard deviations
      if (iter < burn / 2) {
        level <- get.level(old.a, cuts)
        for (j in 1:length(mh.a)) {
          acc.a[j] <- acc.a[j] + sum(old.a[level == j] != a[level == j])
          att.a[j] <- att.a[j] + sum(level == j)
          if (att.a[j] > 100) {
            if (acc.a[j] / att.a[j] < 0.3) { mh.a[j] <- mh.a[j] * 0.9 }
            if (acc.a[j] / att.a[j] > 0.6) { mh.a[j] <- mh.a[j] * 1.1 }
            acc.a[j] <- att.a[j] <- 0
          }
        }
      }

      # update alpha
      if (!alpha.fix) {
        alpha.update <- updateAlpha(y = y, kernel = kernel, a = a,
                                    alpha = alpha, z = z, w = w,
                                    wz.star = wz.star, IDs = IDs, 
                                    alpha.a = alpha.a, alpha.b = alpha.b,
                                    cur.lly = cur.lly, cur.llps = cur.llps,
                                    mid.points=mid.points, bin.width=bin.width,
                                    acc=acc.alpha, att=att.alpha, mh=mh.alpha)

        alpha     <- alpha.update$alpha
        wz.star   <- alpha.update$wz.star
        kernel    <- alpha.update$kernel
        cur.lly   <- alpha.update$cur.lly
        cur.llps  <- alpha.update$cur.llps
        att.alpha <- alpha.update$att
        acc.alpha <- alpha.update$acc

        if (iter < burn / 2) {
          mh.update <- mhUpdate(acc = acc.alpha, att = att.alpha, mh = mh.alpha,
                                nattempts = alpha.attempts)
          acc.alpha <- mh.update$acc
          att.alpha <- mh.update$att
          mh.alpha  <- mh.update$mh
        }
      }  # fi !alpha.fix

      # update rho
      if (!rho.fix) {
        rho.update <- updateRho(y = y, kernel = kernel, a = a, alpha = alpha,
                                cur.lly = cur.lly, w = w, z = z,
                                wz.star = wz.star, dw2 = dw2, rho = rho,
                                logrho.m = logrho.m, logrho.s = logrho.s,
                                rho.upper = rho.upper, acc = acc.rho,
                                att = att.rho, mh = mh.rho)
        rho     <- rho.update$rho
        w       <- rho.update$w
        wz.star <- rho.update$wz.star
        kernel  <- rho.update$kernel
        cur.lly <- rho.update$cur.lly
        att.rho <- rho.update$att
        acc.rho <- rho.update$acc

        if (iter < burn / 2) {
          mh.update <- mhUpdate(acc = acc.rho, att = att.rho, mh = mh.rho,
                                nattempts = rho.attempts)
          acc.rho   <- mh.update$acc
          att.rho   <- mh.update$att
          mh.rho    <- mh.update$mh
        }
      }  # fi rho.fix
    }

    }  # end thin

    if ((iter > burn) & predictions) {
      yp <- rRareBinarySpat(x = x.pred, s = s.pred, knots = knots,
                            beta = beta, xi = xi, alpha = alpha, rho = rho,
                            thresh = thresh, dw2 = dw2p, a = a)
      keepers.y.pred[(iter - burn), , ] <- yp$y
    }

    # storage
    keepers.beta[iter, ] <- beta
    keepers.xi[iter]     <- xi
    if (spatial) {
      keepers.a[iter, , ]  <- a
      keepers.alpha[iter]  <- alpha
      keepers.rho[iter]    <- rho
    }
    keepers.lly[iter] <- sum(logLikeY(y = y, kernel = kernel))

    if (iter %% update == 0) {
      acc.rate.beta  <- round(acc.beta / att.beta, 3)
      acc.rate.xi    <- round(acc.xi / att.xi, 3)
      if (spatial) {
        acc.rate.alpha <- round(acc.alpha / att.alpha, 3)
        acc.rate.rho   <- round(acc.rho / att.rho, 3)
      }

      if (iterplot) {
        if (iter < burn) {
          begin <- max(1, (iter - 2000))
        } else {
          begin <- burn
        }

        par(mfrow=c(2, 4))
        plot(keepers.beta[begin:iter, 1], type="l", main=bquote(beta[0]),
             xlab=round(mh.beta[1], 4), ylab=acc.rate.beta[1])
        if (np > 1) {
          plot(keepers.beta[begin:iter, 2], type="l", main=bquote(beta[1]),
               xlab=round(mh.beta[2], 4), ylab=acc.rate.beta[2])
        }
        if (np > 2) {
          plot(keepers.beta[begin:iter, 3], type="l", main=bquote(beta[2]),
               xlab=round(mh.beta[3], 4), ylab=acc.rate.beta[3])
        }
        plot(keepers.xi[begin:iter], type="l", main=bquote(xi),
             xlab=round(mh.xi, 4), ylab=acc.rate.xi)

        if (spatial) {
          plot(keepers.alpha[begin:iter], type="l", main=bquote(alpha),
               xlab=round(mh.alpha, 4), ylab=acc.rate.alpha)
          plot(keepers.rho[begin:iter], type="l", main=bquote(rho),
               xlab=round(mh.rho, 4), ylab=acc.rate.rho)

          for (i in 1:4) {
            plot(keepers.a[begin:iter, i, 1], type="l",
                 main=paste("knot ", i, ", day 1", sep=""),
                 xlab="", ylab="")
          }
        }
      }
      cat("\t Iter", iter, "\n")
    }

  } # end iter

  return.iters <- (burn + 1):iters

  if (spatial) {
    results <- list(beta = keepers.beta[return.iters, , drop = F],
                    xi = keepers.xi[return.iters],
                    a = keepers.a[return.iters, , , drop = F],
                    alpha = keepers.alpha[return.iters],
                    rho = keepers.rho[return.iters],
                    lly = keepers.lly[return.iters],
                    y.pred = keepers.y.pred)
  } else {
    results <- list(beta = keepers.beta[return.iters, , drop = F],
                    xi = keepers.xi[return.iters],
                    a = NULL, alpha = NULL, rho = NULL,
                    lly = keepers.lly[return.iters],
                    y.pred = keepers.y.pred)
  }

  return(results)
}  # end mcmc
