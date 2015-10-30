mcmc.gev.HMC <- function(y, s, x, knots = NULL, 
                         beta.init = NULL, beta.mn = 0, beta.sd = 20, 
                         beta.eps = 0.1, beta.attempts = 50, beta.fix = FALSE,
                         xi.init = NULL, xi.mn = 0, xi.sd = 0.5, 
                         xi.eps = 0.1, xi.attempts = 50, xi.fix = TRUE,
                         a.init = 10, a.eps = 0.2, a.attempts = 100, 
                         a.fix = FALSE,
                         b.init = 0.5, b.eps = 0.2, b.attempts = 100,
                         b.fix = FALSE,
                         alpha.init = 0.5, alpha.eps = 0.1, alpha.fix = FALSE,
                         rho.init = 1, logrho.mn = -1, logrho.sd = 2, 
                         rho.eps = 0.1, rho.fix = FALSE,
                         threads = 1, iterplot=FALSE, iters=50000, burn=10000, 
                         update=100, thin=1, thresh = 0) {

  library(fields)
  
  # pre-processing for data list
  if (is.null(dim(y))){  # always want a matrix for the responses
    y <- matrix(y, length(y), 1)
  }
  
  if (is.null(knots)) {
    print("using s for knots")
    knots <- s
  }
  
  x <- adjustX(x=x, y=y)  # turn into a stacked matrix (ns * nt x np)
  
  ns     <- nrow(y)
  nt     <- ncol(y)
  np     <- ncol(x)
  nknots <- dim(knots)[1]
  nkt    <- nknots * nt
  
  # create data list
  data  <- list(y = y, x = x, s = s, knots = knots)
  
  # list for miscellaneous things we want to store
  # distance squared
  dw2    <- as.matrix(rdist(s, knots))^2  # dw2 is ns x nknots
  dw2[dw2 < 1e-6] <- 0
  
  others <- list(A.cutoff = max(sqrt(dw2)), thresh = 0, dw2 = dw2)
  
  # initialize MCMC params
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
  
  if (length(a.init) > 1) {
    a <- a.init
  } else {
    a <- matrix(a.init, nknots, nt)
  }
  a.init <- a  # used to make sure all random effects are moving after updateA
  if (is.null(IDs)) {
    # if not specified, make the list of sites that are impacted by each knot
    if (is.null(A.cutoff)) {
      A.cutoff <- 2 * max(dw2)
      print("Including all sites for every knot")
    }
    IDs <- getIDs(dw2, A.cutoff)
  }
  
  if (length(b.init) > 1) {
    b <- b.init
  } else {
    b <- matrix(b.init, nknots, nt)
  }
  
  alpha <- alpha.init
  
  # get the initial set of weights for the sites and knots
  rho    <- rho.init
  
  # create lists for MCMC parameters
  beta  <- list(cur = beta.init, att = 0, acc = 0, eps = beta.eps, 
                mn = beta.mn, sd = beta.sd, attempts = beta.attempts)
  xi    <- list(cur = xi.init, att = 0, acc = 0, eps = xi.eps, 
                mn = xi.mn, sd = xi.sd, attempts = xi.attempts)
  a     <- list(cur = a.init, att = 0, acc = 0, eps = a.eps, infinite = 0,
                attempts = a.attempts)
  b     <- list(cur = b.init, att = 0, acc = 0, eps = b.eps, infinite = 0,
                attempts = b.attempts)
  alpha <- list(cur = alpha.init, att = 0, acc = 0, eps = alpha.eps, 
                infinite = 0, attempts = alpha.attempts)
  rho   <- list(cur = rho.init, att = 0, acc = 0, eps = rho.eps, 
                mn = logrho.mn, sd = logrho.sd, attempts = rho.attempts)
  
  # get calculated list
  calc <- list()
  calc$x.beta <- getXBeta(y = data$y, x = data$x, beta = beta$cur)
  calc$z      <- getZ(xi = xi$cur, x.beta = calc$x.beta, thresh = others$thresh)
  calc$lz     <- log(calc$z)
  calc$w      <- getW(rho = rho$cur, dw2 = others$dw2, A.cutoff = others$A.cutoff)
  calc$lw     <- log(calc$w)
  calc$w.star <- getWStar(alpha = alpha$cur, w = calc$w)
  calc$aw     <- getAW(a = a$cur, w.star = calc$w.star)
  calc$theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
  
  # storage for MCMC
  storage.a     <- array(NA, dim=c(niters, nknots, nt))
  storage.b     <- array(NA, dim=c(niters, nknots, nt))
  storage.alpha <- rep(NA, niters)
  storage.beta  <- rep(NA, niters)
  storage.rho   <- rep(NA, niters)
  storage.prob  <- array(NA, dim = c(niters, ns, nt))
    
  
  tic.1 <- proc.time()
  for (i in 1:niters) {
    beta$att <- beta$att + 1
    MHout <- updateBeta(data = data, beta = beta, xi = xi, alpha = alpha, 
                        calc = calc, others = others)
    if (MHout$accept) {
      beta$acc     <- beta$acc + 1
      beta$cur     <- MHout$q
      calc$x.beta  <- getXBeta(y = data$y, x = data$x, beta = beta$cur)
      calc$z       <- getZ(xi = xi$cur, x.beta = calc$x.beta, thresh = others$thresh)
      calc$lz      <- log(calc$z)
      calc$theta   <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
    }
    if (beta$att > 100) {
      beta.rate <- beta$acc / beta$att
      if (beta.rate < 0.20) {
        beta$eps <- beta$eps * 0.8
      } else if (beta.rate > 0.60) {
        beta$eps <- beta$eps * 1.2
      }
      beta$acc <- beta$att <- 0
    }
    
    # random effect
    a$att <- a$att + 1
    q <- log(a$cur)
    HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q, 
                   epsilon = a$eps, L = 20, 
                   data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                   rho = rho, calc = calc, others = others, this.param = "a")
    if (HMCout$accept) {
      a$acc <- a$acc + 1
      a$cur <- exp(HMCout$q)
      calc$aw  <- getAW(a = a$cur, w.star = calc$w.star)
      calc$theta <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
    }
    a$infinite <- a$infinite + HMCout$infinite
    if (a$infinite > 50) {
      print("reducing a$eps")
      a$eps <- a$eps * 0.8
      a$infinite <- 0
    }
    
    q <- transform$logit(b$cur)
    b$att <- b$att + 1
    HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q, epsilon = b$eps, 
                   L = 10, 
                   data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                   rho = rho, calc = calc, others = others, this.para = "b")
    if (HMCout$accept) {
      b$acc <- b$acc + 1
      b$cur <- transform$inv.logit(HMCout$q)
    }
    b$infinite <- b$infinite + HMCout$infinite
    if (b$infinite > 50) {
      print("reducing b$eps")
      b$eps <- b$eps * 0.8
      b$infinite <- 0
    }
    
    # spatial dependence
    alpha$att <- alpha$att + 1
    q <- transform$logit(alpha$cur)
    HMCout  <- HMC(neg_log_post_alpha, neg_log_post_grad_alpha, q, 
                   epsilon = alpha$eps, L = 10, 
                   data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                   rho = rho, calc = calc, others = others, this.param = "alpha")
    if (HMCout$accept) {
      alpha$acc   <- alpha$acc + 1
      alpha$cur   <- transform$inv.logit(HMCout$q)
      calc$w.star <- getWStar(alpha = alpha$cur, w = calc$w)
      calc$aw     <- getAW(a = a$cur, w.star = calc$w.star)
      calc$theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
    }
    
    alpha$infinite <- alpha$infinite + HMCout$infinite
    if (alpha$infinite > 50) {
      print("reducing alpha$eps")
      alpha$eps <- alpha$eps * 0.8
      alpha$infinite <- 0
    }
    
    rho$att <- rho$att + 1
    MHout <- updateRho(data = data, a = a, alpha = alpha, rho = rho, calc = calc, 
                       others = others)
    if (MHout$accept) {
      rho$acc     <- rho$acc + 1
      rho$cur     <- MHout$q
      calc$w      <- getW(rho = rho$cur, dw2 = others$dw2, 
                          A.cutoff = others$A.cutoff)
      calc$lw     <- log(calc$w)
      calc$w.star <- getWStar(alpha = alpha$cur, w = calc$w)
      calc$aw     <- getAW(a = a$cur, w.star = calc$w.star)
      calc$theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
    }
    
    if (rho$att > 100) {
      rho.rate <- rho$acc / rho$att
      if (rho.rate < 0.20) {
        rho$eps <- rho$eps * 0.8
      } else if (rho.rate > 0.60) {
        rho$eps <- rho$eps * 1.2
      }
      rho$acc <- rho$att <- 0
    }
    
    if (iter < burn / 2) {
      mh.update <- mhUpdate(acc = acc.rho, att = att.rho, mh = mh.rho,
                            nattempts = rho.attempts)
      acc.rho   <- mh.update$acc
      att.rho   <- mh.update$att
      mh.rho    <- mh.update$mh
    }
    
    storage.a[i, , ] <- a$cur
    storage.b[i, , ] <- b$cur
    storage.alpha[i] <- alpha$cur
    storage.beta[i]  <- beta$cur
    storage.rho[i]   <- rho$cur
    storage.prob[i, , ] <- 1 - exp(-calc$theta)
    
    if (i %% update == 0) {
      print(paste("iter:", i, "of", niters, sep=" "))
      if (iterplot) {
        start <- max(i - 5000, 1)
        end   <- i
        par(mfrow=c(4, 5))
        plot.idx <- seq(1, 18, by = 2)
        for (idx in plot.idx){
          plot(log(storage.a[1:i, idx, 1]), type = "l", 
               main = round(log(gen$a[idx, 1]), 2), 
               xlab = round(a$acc / a$att, 3))
        }
        plot.idx <- seq(1, 16, by = 2)
        for (idx in plot.idx){
          plot(storage.b[1:i, idx, 1], type = "l", 
               xlab = round(b$acc / b$att, 3))
        }
        
        plot(storage.beta[start:end], type = "l", main = bquote(beta[0]),
             xlab = round(beta$acc / beta$att, 2), ylab = round(beta$eps, 3))
        plot(storage.alpha[start:end], type = "l", main = bquote(alpha),
             xlab = round(alpha$acc / alpha$att, 2), ylab = "")
        plot(storage.rho[start:end], type = "l", main = bquote(rho),
             xlab = round(rho$acc / rho$att, 2), ylab = round(rho$eps, 3))
      }
    }
  }
  toc <- proc.time()
  
  return.iters <- (burn + 1):iters
  results <- list(beta = storage.beta[return.iters, , drop = F],
                  xi = storage.xi[return.iters],
                  a = storage.a[return.iters, , , drop = F],
                  b = storage.b[return.iters, , , drop = F],
                  alpha = storage.alpha[return.iters],
                  rho = storage.rho[return.iters],
                  A.cutoff = A.cutoff)
}