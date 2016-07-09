make.B <- function(d, rho) {
  b <- exp(-(d / rho)^2)
  b[b <= 1e-100] <- 0
  b2 <- rowSums(b^2)
  b <- sweep(b, 1, sqrt(b2), "/")
  if (any(is.nan(b))) {  # nan comes from 0 / 0 which means bw is very small
    cat("\t Some sites are independent from others due to small bandwidth \n")
  }
  b[is.nan(b)] <- 0
  return(b)
}

# predictions
pred.spprob <- function(mcmcoutput, X.pred, s.pred, knots,
                        start = 1, end = NULL, update = NULL) {

  if (is.null(end)) {
    end <- nrow(mcmcoutput$beta)
  }

  # bookkeeping
  np    <- nrow(s.pred)
  iters <- end - start + 1
  dp    <- as.matrix(rdist(s.pred, knots))
  prob.success <- matrix(NA, nrow=iters, ncol=np)

  beta <- mcmcoutput$beta[start:end, , drop = F]
  bw   <- mcmcoutput$bw[start:end]
  alpha <- mcmcoutput$alpha[start:end, , drop = F]

  for (i in 1:iters) {
    # beta and bandwidth come directly from output
    z.pred <- X.pred %*% beta[i, ] + make.B(dp, bw[i]) %*% alpha[i, ]
    prob.success[i, ] <- pnorm(z.pred)

    if (!is.null(update)) {
      if (i %% update == 0) {
        cat("\t Iter", i, "\n")
      }
    }
  }

  return(prob.success)
}

dic.spprob <- function(mcmcoutput, Y, X, s, knots,
                       start = 1, end = NULL, update = NULL) {
  if (is.null(end)) {
    end <- length(mcmcoutput$bw)
  }

  d <- rdist(s, knots)

  beta  <- mcmcoutput$beta[start:end, , drop = F]
  bw    <- mcmcoutput$bw[start:end]
  alpha <- mcmcoutput$alpha[start:end, , drop = F]

  betabar  <- apply(beta, 2, mean)
  bwbar    <- mean(bw)
  alphabar <- apply(alpha, 2, mean)
  Bbar     <- make.B(d, bwbar)
  BAbar    <- Bbar %*% alphabar
  XBbar    <- X %*% betabar

  dthetabar <- -2 * sum(dbinom(Y, 1, pnorm(XBbar + BAbar), log = TRUE))

  niters <- length(start:end)
  dbar   <- 0
  for (i in 1:end) {
    beta.i  <- beta[i, ]
    bw.i    <- bw[i]
    alpha.i <- alpha[i, ]
    B       <- make.B(d, bw.i)
    BA      <- B %*% alpha.i
    XB      <- X %*% beta.i

    dbar  <- dbar - 2 * sum(dbinom(Y, 1, pnorm(XB + BA), log = TRUE)) / niters
    if (!is.null(update)) {
      if (i %% update == 0) {
        cat("\t Iter", i, "\n")
      }
    }
  }

  pd  <- dbar - dthetabar
  dic <- dbar + pd

  results <- list(dic = dic, pd = pd, dbar = dbar, dthetabar = dthetabar)
  return(results)
}

# predictions
probit <- function(Y, X, s, knots, sp=NULL, Xp=NULL,
                   logbw.mn=-1, logbw.sd=2, eps=0.01, a=0.1, b=0.1,
                   iters=5000, burn=1000, thin=1, update=2, iterplot=FALSE){

    par(mfrow=c(1, 1))
    tick <- proc.time()[3]

   #Bookkeeping

    n     <- length(Y)
    p     <- ncol(X)
    m     <- nrow(knots)

    d     <- rdist(s,knots)
    if(!is.null(sp)){
      np  <- nrow(sp)
      dp  <- rdist(sp,knots)
    }

    LOW  <- ifelse(Y == 1, 0, -Inf)
    HIGH <- ifelse(Y == 1, Inf, 0)
    Y    <- ifelse(Y == 1, 2, -2)

   #initial values

    beta  <- rep(0, p)
    alpha <- rep(0, m)
    bw    <- exp(logbw.mn)
    taua  <- 1

    B     <- make.B(d, bw)
    XB    <- X %*% beta
    BA    <- B %*% alpha

   # Keep track of stuff

    PRED <- 0
    FIT  <- 0

    # trying to keep memory usage low, so only samples saving post burnin
    keep.bw             <- rep(0, iters - burn)
    keep.beta           <- matrix(0, iters - burn, p)
    keep.taua           <- rep(0, iters - burn)
    keep.alpha          <- matrix(0, iters - burn, m)
    colnames(keep.beta) <- colnames(X)

   # Preprocessing

    tXX    <- t(X) %*% X

    t_last_plot <- proc.time()[3]

    for(iter in 1:iters){

      for(ttt in 1:thin){

       # LATENT CONTINUOUS VARIABLES:

        MMM  <- XB + BA
        lp   <- pnorm(LOW, MMM, 1)
        up   <- pnorm(HIGH, MMM, 1)
        U    <- lp + (up - lp) * runif(n)
        Y    <- qnorm(U, MMM, 1)

       # BETA

        VVV  <- chol2inv(chol(tXX + diag(eps, p)))
        MMM  <- crossprod(X, Y - BA)
        beta <- VVV %*% MMM + t(chol(VVV)) %*% rnorm(p)
        XB   <- X %*% beta

       # ALPHA
        VVV   <- crossprod(B) + diag(taua, m)
        VVV   <- chol2inv(chol(VVV))
        MMM   <- crossprod(B, Y - XB)
        alpha <- VVV %*% MMM + t(chol(VVV)) %*% rnorm(m)
        BA    <- B %*% alpha

       # TAUA

        taua <- rgamma(1, m / 2 + a, sum((alpha^2)) / 2 + b)

       # BW

        canbw <- exp(rnorm(1, log(bw), 0.1))
        canB  <- make.B(d, canbw)
        canBA <- canB %*% alpha
        R     <- sum(dnorm(Y, canBA + XB, 1, log=TRUE)) -
                 sum(dnorm(Y, BA + XB, 1, log=TRUE)) +
                 dnorm(log(canbw), logbw.mn, logbw.sd, log=TRUE) -
                 dnorm(log(bw), logbw.mn, logbw.sd, log=TRUE)

        if (!is.na(R)) { if (log(runif(1))<R) {
          bw <- canbw
          B  <- canB
          BA <- canBA
        } }

      }


      #Plot current iteration
      if (iterplot) {
        if(proc.time()[3] - t_last_plot > update | iter == iters){
          plot(s, pch=19,
               col=ifelse(Y > 0, 2, 1),
               cex=2 * pnorm(XB + BA),
               main=paste(iter, "of", iters))
          t_last_plot <- proc.time()[3]
        }
      }

      # Keep track of stuff:

      if(iter > burn){
        nnn <- iters - burn
        FIT <- FIT + pnorm(XB + BA) / nnn

        if(!is.null(Xp)){
          MMM  <- Xp %*% beta + make.B(dp, bw) %*% alpha
          PRED <- PRED + pnorm(MMM) / nnn
        }

        keep.bw[iter - burn]      <- bw
        keep.beta[iter - burn, ]  <- beta
        keep.taua[iter - burn]    <- taua
        keep.alpha[iter - burn, ] <- alpha
      }

      if (iter %% update == 0) {
        cat("      Iter", iter, "of", iters, "\n")
      }
    }

    tock   <- proc.time()[3]
    output <- list(fitted=FIT,
                   pred=PRED,
                   beta=keep.beta,
                   bw=keep.bw,
                   taua=keep.taua,
                   alpha=keep.alpha,
                   minutes=(tock - tick)/60)

return(output)}


############################################################
##########       SIMULATED EXAMPLE
############################################################


# library(fields)
#
# n  <- 1000
# np <- 1000
#
# s     <- cbind(runif(n),runif(n))
# sp    <- cbind(runif(np),runif(np))
# knots <- expand.grid(seq(0,1,.1),seq(0,1,.1))
#
# X     <- matrix(1,n,1)
# Xp    <- matrix(1,np,1)
#
# alpha <- rnorm(nrow(knots))
# beta  <- -.5
# bw    <- 0.1
#
# B     <- make.B(rdist(s,knots),bw)
# Bp    <- make.B(rdist(sp,knots),bw)
#
# M     <- B%*%alpha+X%*%beta
# Z     <- M+rnorm(n)
# FIT   <- pnorm(M)
# Y     <- ifelse(Z>0,1,0)
# MP    <- Bp%*%alpha+Xp%*%beta
# PRED  <- pnorm(MP)
#
# fit<-probit(Y,X,s,knots,Xp=Xp,sp=sp)
#
# par(mfrow=c(2,2))
#
# plot(FIT,fit$fitted)
# abline(0,1,col=2)
#
# plot(PRED,fit$pred)
# abline(0,1,col=2)
#
# plot(fit$beta,type="l")
# abline(beta,0,col=2)
#
# plot(fit$bw,type="l")
# abline(bw,0,col=2)

# slow <- function(X, taua, m) {
#   return (t(X) %*% X + taua * diag(m))
# }
#
# med <- function(X, taua, m) {
#   return (crossprod(X) + taua * diag(m))
# }
#
# fast <- function(X, taua, m) {
#   return (crossprod(X) + diag(taua, m))
# }
#
# fast1 <- function(X, taua, m) {
#   return (chol2inv(chol(crossprod(X) + diag(taua, m))))
# }
#
# fast2 <- function(X, taua, m) {
#   VVV <- crossprod(X) + diag(taua, m)
#   return (chol2inv(chol(VVV)))
# }
#
# slowa <- function(X, res) {
#   return (t(X) %*% res)
# }
#
# fasta <- function(X, res) {
#   return (crossprod(X, res))
# }
#
# microbenchmark(slow(X = B, taua = 0.25, m = ncol(B)),
#                med(X = B, taua = 0.25, m = ncol(B)),
#                fast(X = B, taua = 0.25, m = ncol(B)))
#
# microbenchmark(slowa(X = B, res = res),
#                fasta(X = B, res = res))
#
# microbenchmark(fast1(X = X, taua = 0.25, m = ncol(X)),
#                fast2(X = X, taua = 0.25, m = ncol(X)))