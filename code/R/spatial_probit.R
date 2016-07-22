make.B <- function(d, rho) {
  b <- exp(-(d / rho)^2)
  b[b <= 1e-6] <- 0
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
                        start = 1, thin = 1, end = NULL, update = NULL) {

  if (is.null(end)) {
    end <- nrow(mcmcoutput$beta)
  }

  # bookkeeping
  np     <- nrow(s.pred)
  iters  <- end - start + 1
  dp     <- as.matrix(rdist(s.pred, knots))
  y.pred <- matrix(NA, nrow = floor(iters / thin), ncol=np)

  beta <- mcmcoutput$beta[start:end, , drop = F]
  bw   <- mcmcoutput$bw[start:end]
  alpha <- mcmcoutput$alpha[start:end, , drop = F]

  for (i in 1:iters) {
    if (i %% thin == 0) {
      # beta and bandwidth come directly from output
      z.pred <- X.pred %*% beta[i, ] + make.B(dp, bw[i]) %*% alpha[i, ]
      prob.success <- pnorm(z.pred)
      y.pred[(i / thin), ]  <- rbinom(n = np, size = 1, prob = prob.success)
    }
    
    if (!is.null(update)) {
      if (i %% update == 0) {
        cat("\t Iter", i, "\n")
      }
    }
  }

  return(y.pred)
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
                   logbw.mn=-1, logbw.sd=2, bw.lower = NULL, bw.upper = NULL,
                   eps=0.01, a=0.1, b=0.1,
                   iters=5000, burn=1000, thin=1, update=2, 
                   keep.burn = FALSE, iterplot=FALSE){

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
    
    plot.ys <- sample(which(Y == 1), size = 3)
    plot.ys <- c(plot.ys, sample(which(Y == 0), size = 3))
    
    LOW  <- ifelse(Y == 1, 0, -Inf)
    HIGH <- ifelse(Y == 1, Inf, 0)
    Y    <- ifelse(Y == 1, 2, -2)

   #initial values
    if (is.null(bw.lower)) {
      bw.lower <- 1e-6
    }
    if(is.null(bw.upper)) {
      bw.upper <- max(d)
    }
    beta  <- rep(0, p)
    alpha <- rep(0, m)
    bw    <- (bw.upper - bw.lower) / 2
    taua  <- 1

    B     <- make.B(d, bw)
    XB    <- X %*% beta
    BA    <- B %*% alpha

   # Keep track of stuff

    PRED <- 0
    FIT  <- 0

    # trying to keep memory usage low, so only samples saving post burnin
    keep.bw     <- rep(0, iters)
    keep.beta   <- matrix(0, iters, p)
    keep.taua   <- rep(0, iters)
    keep.alpha  <- matrix(0, iters, m)
    prob.y      <- matrix(0, iters, n)
    colnames(keep.beta) <- colnames(X)

   # Preprocessing
    att.bw <- acc.bw <- MH.bw <- 0.1

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
        att.bw <- att.bw + 1
        bw.star <- transform$logit(bw, lower = bw.lower, upper = bw.upper)
        canbw.star <- rnorm(1, bw.star, MH.bw)
        canbw <- transform$inv.logit(canbw.star, lower = bw.lower,
                                     upper = bw.upper)
        # canbw <- exp(rnorm(1, log(bw), MH.bw))
        canB  <- make.B(d, canbw)
        canBA <- canB %*% alpha
        R     <- sum(dnorm(Y, canBA + XB, 1, log=TRUE)) -
                 sum(dnorm(Y, BA + XB, 1, log=TRUE)) +
                 log(canbw - bw.lower) + log(bw.upper - canbw) - # Jacobian
                 log(bw - bw.lower) - log(bw.upper - bw)
                 # dnorm(log(canbw), logbw.mn, logbw.sd, log=TRUE) -
                 # dnorm(log(bw), logbw.mn, logbw.sd, log=TRUE)

        if (!is.na(R)) { if (log(runif(1))<R) {
          acc.bw <- acc.bw + 1
          bw <- canbw
          B  <- canB
          BA <- canBA
        } }

      }
      
      if (iter < burn / 2) {
        this.update <- mhUpdate(acc = acc.bw, att = att.bw, MH = MH.bw,
                                nattempts = 200)
        att.bw <- this.update$att
        acc.bw <- this.update$acc
        MH.bw  <- this.update$MH
      }

      # Keep track of stuff:
      keep.bw[iter] <- bw
      keep.beta[iter, ] <- beta
      keep.taua[iter]   <- taua
      keep.alpha[iter, ]  <- alpha
      prob.y[iter, ] <- pnorm(XB + BA)
        
      #Plot current iteration
      if (iter %% update == 0) {
        cat("      Iter", iter, "of", iters, "\n")
        
        if (iter < burn) {
          start <- max(0, iter - 5000) 
        } else {
          start <- burn + 1
        }
        
        if (iterplot) {
        par(mfrow = c(3, 3))
        plot(keep.beta[start:iter, 1], type = "l", main = "beta")
        plot(keep.taua[start:iter], type = "l", main = "taua")
        plot(keep.bw[start:iter], type = "l", main = "bw",
             xlab = paste("acc.rate = ", round(acc.bw / att.bw, 3), sep = ""),
             ylab = paste("MH = ", round(MH.bw, 3)))
        
        
        for (i in 1:length(plot.ys)) {
          plot(prob.y[start:iter, plot.ys[i]], type = "l", 
               main = paste("P(Y", plot.ys[i], " = 1)"),
               ylim = c(0, 1))
        }
        # if(proc.time()[3] - t_last_plot > update | iter == iters){
        #   plot(s, pch=19,
        #        col=ifelse(Y > 0, 2, 1),
        #        cex=2 * pnorm(XB + BA),
        #        main=paste(iter, "of", iters))
        #   t_last_plot <- proc.time()[3]
        # }
        }
      }
      
      
      if(iter > burn){
        nnn <- iters - burn
        FIT <- FIT + pnorm(XB + BA) / nnn

        if(!is.null(Xp)){
          MMM  <- Xp %*% beta + make.B(dp, bw) %*% alpha
          PRED <- PRED + pnorm(MMM) / nnn
        }

        # keep.bw[iter - burn]      <- bw
        # keep.beta[iter - burn, ]  <- beta
        # keep.taua[iter - burn]    <- taua
        # keep.alpha[iter - burn, ] <- alpha
      }


    }
    
    if (!keep.burn) {
      return.iters <- (burn + 1):iters
    } else {
      return.iters <- 1:iters
    }
    
    tock   <- proc.time()[3]
    output <- list(fitted=FIT,
                   pred=PRED,
                   beta=keep.beta[return.iters, , drop = FALSE],
                   bw=keep.bw[return.iters],
                   taua=keep.taua[return.iters],
                   alpha=keep.alpha[return.iters, , drop = FALSE],
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