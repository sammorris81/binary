#############################################
#
# Required inputs:
#
#  Y        := n-vector of binary responses
#  s        := nx2 lat/long matrix
#
# Optional inputs 
#
#  X        := n x p matrix of covariates
#  nknots   := number of knots in each direction (there are nknots^2 knots)
#
# Model
#
#  Prob(Y=1) = expit(X%*%beta + W%*%beta)
#  W         = make.W(s,nknots) are kernel functions
#  beta      ~ CAR with mean mu, precision tau, and dependence rho 
#  alpha[j]  ~ N(0,1/eps)
#  mu        ~ N(0,eps)
#  tau       ~ gamma(a,b)
#  rho       ~ beta(ar,br)
#
# Output:
#
#  alpha     := posterior draws of alpha
#  hyper     := posterior draws of mu, 1/sqrt(tau), and rho
#  beta      := posterior mean of beta
#  fitted    := posterior mean of Prob(Y=1)
#
#############################################

source("hmc_logit.R") 

spatial_logit <- function(Y, s, knots = NULL, X = NULL,
                          nknots = 10, MH = c(0.5, 0.1),
                          a = 0.1, b = 0.1, eps = 0.01, ar = 5, br = 1,
                          iters = 2000, burn = 1000, iterplot = FALSE, 
                          update = 100){
  library(fields)
  library(emulator)
  library(spam)
  
  #Bookkeeping
  
  n    <- length(Y)
  
  # Set up adjacency between knots
  if (is.null(knots)) {
    W     <- make.W(s, nknots)
    ADJ   <- W$ADJ
    knots <- W$knots
    L     <- nrow(knots)
    W     <- W$W
  } else {
    W     <- make.W.knots(s, knots)
    ADJ   <- W$ADJ
    knots <- knots
    L     <- nrow(knots)
    W     <- W$W
  }
  nknots <- nrow(knots)
  M     <- rowSums(ADJ)
  adj1  <- row(ADJ)[ADJ == 1]
  adj2  <- col(ADJ)[ADJ == 1]
  low   <- adj1 < adj2
  adj1  <- adj1[low]
  adj2  <- adj2[low]
  
  M2       <- sum(M)
  ADJ1     <- rowSums(ADJ)
  ADJ2     <- sum(ADJ1)
  
  # Initial value and discretized prior for rho
  rho      <- 0.9
  rho.grid <- qbeta(seq(1 / 1000, 1 - 1 / 1000, length = 1000), ar, br)
  C        <- sweep(ADJ, 1, sqrt(M), "/")
  C        <- sweep(C, 2, sqrt(M), "/")
  d        <- eigen(as.spam(C))$values
  logd     <- rho.grid
  for (j in 1:length(rho.grid)) {
    logd[j] <- sum(log(1 - rho.grid[j] * d))
  }
  
  
  
  # Initial values and MH covariance for the fixed effects
  
  p    <- 0
  Xa   <- 0
  if (!is.null(X)) {
    if (is.vector(X)) {X <- matrix(X, n, 1)}
    p     <- ncol(X)
    beta  <- glm(Y ~ X, family = "binomial")
    if (p == 1) {
      beta <- beta$coef[-1]
      Xb   <- X * beta
    }
    if (p > 1) {
      Cb    <- summary(beta)$cov.un
      Pb    <- t(chol(Cb[-1, -1]))
      beta  <- beta$coef[-1]
      Xb    <- X %*% beta
    }
  } else {
    X    <- matrix(1, n, 1)
    p    <- 1
    beta <- glm(Y ~ X, family = "binomial")
    beta <- beta$coef[1]
    Xb   <- X %*% beta
  }
  
  # Initial values for the spatial random effects
  
  K <- round(n / L) + 1
  phat <- rep(0, n)
  for (j in 1:n) {
    ddd     <- (s[j, 1] - s[, 1])^2 + (s[j, 2] - s[, 2])^2
    phat[j] <- mean(Y[rank(ddd) <= K])
  }
  phat      <- ifelse(phat < 0.01, 0.01, phat)
  phat      <- ifelse(phat > 0.99, 0.99, phat)
  wahat     <- logit(phat) - Xb
  tWW       <- crossprod(W)
  ddd       <- diag(tWW)
  diag(tWW) <- ddd + max(ddd) / 1000
  
  alpha     <- as.vector(solve(tWW, t(W) %*% wahat))
  tau       <- 1 / var(alpha)
  mu        <- mean(alpha)
  rm(tWW, phat, wahat, ddd)
  
  # More bookkeeping
  Wa    <- W %*% alpha
  curp  <- expit(Xb + Wa)
  curll <- sum(dbinom(Y, 1, curp, log = TRUE))
  
  # Set up the storage
  
  keep.hyper <- matrix(0, iters, 3)
  keep.beta  <- matrix(0, iters, max(c(1, p)))
  keep.alpha <- matrix(0, iters, nknots)
  MNA        <- 0
  VARA       <- 0
  fitted     <- 0
  
  colnames(keep.hyper) <- c("alpha_mu", "alpha_sd", "alpha_rho")
  if (p > 0) {colnames(keep.beta) <- colnames(X)}
  
  
  # Metropolis tuning
  beta.mh  <- list(acc = 0, att = 0, eps = MH[1], attempts = 50, 
                   target.l = 0.5, target.u = 0.8)
  alpha.mh <- list(acc = 0, att = 0, eps = MH[2], attempts = 50,
                   target.l = 0.5, target.u = 0.8)
  
  # Start the timer
  
  tick    <- proc.time()[3]
  
  # MCMC!
  
  for (iter in 1:iters) {
    
    # Update beta (Random-walk Metropolis)
    
#     if (p > 0) {
#       beta.mh$att <- beta.mh$att + 1
#       if (p == 1) {
#         canb    <- rnorm(1, beta, beta.mh$eps)
#         canXb   <- X * canb
#       }
#       if (p > 1) {
#         canb    <- beta + as.vector(Pb %*% rnorm(p, 0, beta.mh$eps))
#         canXb   <- X %*% canb
#       }
#       canp     <- make.prob(canXb + Wa)
#       canll    <- sum(dbinom(Y, 1, canp, log=TRUE))
#       R        <- canll - curll +
#         sum(dnorm(canb, 0, 1 / sqrt(eps), log = TRUE)) -
#         sum(dnorm(beta, 0, 1 / sqrt(eps), log=TRUE))
#       if (!is.na(R)) {if (log(runif(1)) < R) {
#         beta.mh$acc <- beta.mh$acc + 1
#         beta   <- canb
#         Xb     <- canXb
#         curll  <- canll
#         curp   <- canp
#       }}   
#     }
    # beta.mh$att <- beta.mh$att + 1
    # others <- list(Y = Y, X = X, Wa = Wa, beta.mn = rep(0, p), 
    #                beta.sd = 1 / sqrt(eps))
    # HMCout <- logitHMC(neg_log_post_beta.logit, neg_log_post_grad_beta.logit, 
    #                    beta, epsilon = beta.mh$eps, L = 10, others)
    # if (HMCout$accept) {
    #   beta.mh$acc <- beta.mh$acc + 1
    #   beta <- HMCout$q
    #   Xb   <- X %*% beta
    #   curp <- make.prob(Xb + Wa)
    #   curll <- sum(dbinom(Y, 1, curp, log = TRUE))
    # }
    
    # Update alpha (HMC)
    Q       <- as.spam(diag(M) - rho * ADJ)
    alpha.mh$att  <- alpha.mh$att + 1
    others  <- list(Y = Y, Q = Q, Xb = Xb, W = W, mu = mu, tau = tau)
    HMCout  <- logitHMC(neg_log_post_alpha.logit, neg_log_post_grad_alpha.logit, 
                        alpha, epsilon = alpha.mh$eps, L = 10, others)
    if (HMCout$accept) {
      alpha.mh$acc <- alpha.mh$acc + 1
      alpha  <- HMCout$q
      Wa     <- W %*% alpha
      curp   <- make.prob(Xb + Wa)
      curll  <- sum(dbinom(Y, 1, curp, log = TRUE))
    }
    
    
    # Update the CAR hyperparameters
    
    VVV <- tau * (M2 - rho * ADJ2) + eps
    MMM <- tau * sum((M - rho * ADJ1) * alpha)
    mu  <- rnorm(1, MMM / VVV, 1 / sqrt(VVV))
    
    R   <- alpha - mu
    SS  <- sum(M * R^2) - rho * 2 * sum(R[adj1] * R[adj2])
    tau <- rgamma(1, L / 2 + a, SS / 2 + b)
    
    R   <- 0.5 * logd + tau * sum(R[adj1] * R[adj2]) * rho.grid
    rho <- sample(rho.grid, 1, prob = exp(R - max(R)))
    
    # Keep track of stuff
    
    if (iter > burn) {
      MNA    <- MNA    + alpha / (iters - burn)
      VARA   <- VARA   + alpha * alpha / (iters - burn)
      fitted <- fitted + curp / (iters - burn)
    }
    
    keep.alpha[iter, ] <- alpha
    keep.hyper[iter, ] <- c(mu, 1 / sqrt(tau), rho)
    if (p > 0) {keep.beta[iter, ] <- beta}
    
    # Metropolis tuning
    
    if (iter < burn) {
      eps.update <- epsUpdate(beta.mh)
      beta.mh$att   <- eps.update$att
      beta.mh$acc   <- eps.update$acc
      beta.mh$eps   <- eps.update$eps
      
      eps.update <- epsUpdate(alpha.mh)
      alpha.mh$att  <- eps.update$att
      alpha.mh$acc  <- eps.update$acc
      alpha.mh$eps  <- eps.update$eps
    }
    
    # Plot the current iteration
    
    if (iter %% update == 0) {
      cat("      Iter", iter, "of", iters, "\n")
      if (iterplot) {
        par(mfrow = c(ifelse(p > 0, 3, 2), 3))
        for (j in 1:5) {
          plot(keep.alpha[1:iter, j], type = "l", main = paste("a", j),
               xlab = round(alpha.mh$acc / alpha.mh$att, 2), 
               ylab = round(alpha.mh$eps, 4))
        }
        plot(keep.beta[1:iter, 1], type = "l", main = bquote(beta[0]),
             xlab = round(beta.mh$acc / beta.mh$att, 2), 
             ylab = round(beta.mh$eps, 4))
        plot(keep.hyper[1:iter, 1], type = "l")
        plot(keep.hyper[1:iter, 2], type = "l")
        plot(keep.hyper[1:iter, 3], type = "l")
        
      }
    }
  }
  tock <- proc.time()[3]
  
  out <-list(knots = knots, W = W, ADJ = ADJ,
             fitted = fitted,
             alpha.mn = MNA, alpha.var = VARA - MNA^2,
             alpha = keep.alpha, beta = keep.beta, hyper = keep.hyper,
             acc.rate = c(beta.mh$acc / beta.mh$att, alpha.mh$acc / alpha.mh$att), 
             minutes = (tock - tick) / 60)
  
  return(out)
}

logit <- function(p) {log(p / (1 - p))}
expit <- function(x) {1 / (1 + exp(-x))}

make.prob<-function(Xb){
  expit(ifelse(Xb > 10, 10, Xb))
}

neg_log_post_alpha.logit <- function(alpha, o) {
  prob <- make.prob(o$Xb + o$W %*% alpha)
  R    <- alpha - o$mu
  SS   <- sum(R * (o$Q %*% R))
  ll   <- sum(dbinom(o$Y, 1, prob, log = TRUE)) - 0.5 * o$tau * SS
  return(-ll)
}

neg_log_post_grad_alpha.logit <- function(alpha, o) {
  prob <- make.prob(o$Xb + o$W %*% alpha)
  grad <- t(o$W) %*% (prob - o$Y) + o$tau * o$Q %*% (alpha - o$mu)
  return(grad)
}

neg_log_post_beta.logit <- function(beta, o) {
  prob <- make.prob(o$X %*% beta + o$Wa)
  ll <- sum(dbinom(o$Y, 1, prob, log = TRUE))
  ll <- ll + sum(dnorm(beta, o$beta.mn, o$beta.sd, log = TRUE))
  return(-ll)
}

neg_log_post_grad_beta.logit <- function(beta, o) {
  prob <- make.prob(o$X %*% beta + o$Wa)
  grad <- t(o$X) %*% (prob - o$Y) + (beta - o$beta.mn) / o$beta.sd
  return(grad)
}

# update mh settings
epsUpdate <- function(param, lower = 0.8, higher = 1.2) {
  acc <- param$acc
  att <- param$att
  eps <- param$eps
  attempts <- param$attempts
  
  acc.rate     <- acc / att
  these.update <- att > attempts
  these.low    <- (acc.rate < param$target.l) & these.update
  these.high   <- (acc.rate > param$target.u) & these.update
  
  eps[these.low]  <- eps[these.low] * lower
  eps[these.high] <- eps[these.high] * higher
  
  acc[these.update] <- 0
  att[these.update] <- 0
  
  results <- list(att = att, acc = acc, eps = eps)
  return(results)
}

make.W <- function(s, nknots, buffer = 0) {
  r1     <- range(s[, 1])
  r2     <- range(s[, 2])
  r1.min <- floor(r1[1])
  r1.max <- ceiling(r1[2])
  r2.min <- floor(r2[1])
  r2.max <- ceiling(r2[2])
  k1     <- seq(r1.min - buffer * (r1.max - r1.min), r1.max + buffer * (r1.max - r1.min), length = nknots)
  k2     <- seq(r2.max - buffer * (r2.max - r2.min), r2.max + buffer * (r2.max - r2.min), length = nknots)
  knots  <- expand.grid(k1, k2)
  dsk    <- (rdist(s[, 1], knots[, 1]) / (r1.max - r1.min))^2 +
            (rdist(s[, 2], knots[, 2]) / (r2.max - r2.min))^2
  W      <- dnorm(dsk, 0, 0.5)
  W      <- sweep(W, 1, sqrt(rowSums(W^2)), "/")
  loc    <- expand.grid(1:nknots, 1:nknots)
  ADJ    <- ifelse(as.matrix(rdist(loc)) == 1, 1, 0)
  return(list(W = W, ADJ = ADJ, knots = knots))
}

make.W.knots <- function(s, knots) {
  nknots.x <- length(unique(knots[, 1]))
  nknots.y <- length(unique(knots[, 2]))
  r1 <- range(knots[, 1])
  r2 <- range(knots[, 2])
  dsk    <- (rdist(s[, 1], knots[, 1]) / (r1[2] - r1[1]))^2 +
            (rdist(s[, 2], knots[, 2]) / (r2[2] - r2[1]))^2
  W      <- dnorm(dsk, 0, 0.5)
  W      <- sweep(W, 1, sqrt(rowSums(W^2)), "/")
  loc    <- expand.grid(1:nknots.x, 1:nknots.y)
  ADJ    <- ifelse(as.matrix(rdist(loc)) == 1, 1, 0)
  return(list(W = W, ADJ = ADJ))
}

pred.splogit <- function(mcmcoutput, X.pred = NULL, s.pred, knots,
                         start = 1, end = NULL, update = NULL) {
  
  if (is.null(end)) {
    end <- nrow(mcmcoutput$beta)
  }
  
  
  # bookkeeping
  np    <- nrow(s.pred)
  iters <- end - start + 1
  dp    <- as.matrix(rdist(s.pred, knots))
  prob.success <- matrix(NA, nrow=iters, ncol=np)
  
  beta  <- mcmcoutput$beta[start:end, , drop = F]
  alpha <- mcmcoutput$alpha[start:end, , drop = F]
  
  W <- make.W.knots(s.pred, knots)$W
  if (is.null(X.pred)) {
    X.pred <- matrix(1, np, 1)
  }
  
  for (i in 1:iters) {
    # beta and bandwidth come directly from output
    z.pred <- X.pred %*% beta[i, ] + W %*% alpha[i, ]
    prob.success[i, ] <- pnorm(z.pred)
    
    if (!is.null(update)) {
      if (i %% update == 0) {
        cat("\t Iter", i, "\n")
      }
    }
  }
  
  return(prob.success)
}

