make.B <- function(d, rho){
  exp(-(d / rho)^2)
}

#Main MCMC function

probit <- function(Y, X, s, knots, sp=NULL, Xp=NULL,
                   logbw.mn=-1, logbw.sd=2, eps=0.01, a=0.1, b=0.1,
                   iters=5000, burn=1000, thin=1, update=2){

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

    keep.bw             <- rep(0, iters)
    keep.beta           <- matrix(0, iters, p)
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

        VVV  <- solve(tXX + eps * diag(p))
        MMM  <- t(X) %*% (Y - BA)
        beta <- VVV %*% MMM + t(chol(VVV)) %*% rnorm(p)
        XB   <- X %*% beta

       # ALPHA

        VVV   <- solve(t(B) %*% B + taua * diag(m))
        MMM   <- t(B) %*% (Y-XB)
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

        if(log(runif(1))<R){
          bw <- canbw
          B  <- canB
          BA <- canBA
        }

      }


      #Plot current iteration

      if(proc.time()[3] - t_last_plot > update | iter == iters){
        plot(s, pch=19,
             col=ifelse(Y > 0, 2, 1),
             cex=2 * pnorm(XB + BA),
             main=paste(iter, "of", iters))
        t_last_plot <- proc.time()[3]
      }

      # Keep track of stuff:

      if(iter > burn){
        nnn <- iters - burn
        FIT <- FIT + pnorm(XB + BA) / nnn

        if(!is.null(Xp)){
          MMM  <- Xp %*% beta + make.B(dp, bw) %*% alpha
          PRED <- PRED + pnorm(MMM) / nnn
        }

      }

      keep.bw[iter]    <- bw
      keep.beta[iter, ] <- beta
    }

    tock   <- proc.time()[3]
    output <- list(fitted=FIT,
                   pred=PRED,
                   beta=keep.beta,
                   bw=keep.bw,
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
