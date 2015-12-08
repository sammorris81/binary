source("./hmc_gev.R")
source("./aux_gev.R")
source("./update_gev.R")
source("./nll_post_gev.R")

spatial_GEV <- function(y, s, x, knots = NULL, 
                        beta.init = NULL, beta.mn = 0, beta.sd = 20, 
                        beta.eps = 0.1, beta.attempts = 50,
                        xi.init = NULL, xi.mn = 0, xi.sd = 0.5, 
                        xi.eps = 0.1, xi.attempts = 50, xi.fix = TRUE,
                        a.init = 10, a.eps = 0.2, a.attempts = 100, 
                        a.cutoff = NULL,
                        b.init = 0.5, b.eps = 0.2, b.attempts = 100,
                        alpha.init = 0.5, alpha.eps = NULL, alpha.attempts = 50,
                        a.alpha.joint = TRUE,
                        rho.init = 1, logrho.mn = -2, logrho.sd = 1, 
                        rho.eps = 0.1, rho.attempts = 50,
                        threads = 1, iterplot = FALSE, iters = 50000, 
                        burn = 10000, update = 100, thin = 1, thresh = 0) {

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
  
  if (is.null(a.cutoff)) {
    a.cutoff <- max(sqrt(dw2)) + 0.001
    print("Including all sites for every knot")
  }
  
  # IDs is assigned after creating A
  others <- list(a.cutoff = a.cutoff, thresh = 0, dw2 = dw2)
  
  # if not specified, make the list of sites that are impacted by each knot
  others$IDs <- getIDs(dw2, others$a.cutoff)
  
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
                mn = beta.mn, sd = beta.sd, attempts = beta.attempts,
                target.l = 0.25, target.u = 0.5)
  xi    <- list(cur = xi.init, att = 0, acc = 0, eps = xi.eps, 
                mn = xi.mn, sd = xi.sd, attempts = xi.attempts,
                target.l = 0.25, target.u = 0.5, fix = xi.fix)
  a     <- list(cur = a.init, att = 0, acc = 0, eps = a.eps, infinite = 0,
                attempts = a.attempts, target.l = 0.5, target.u = 0.8)
  b     <- list(cur = b.init, att = 0, acc = 0, eps = b.eps, infinite = 0,
                attempts = b.attempts, target.l = 0.5, target.u = 0.8)
  alpha <- list(cur = alpha.init, att = 0, acc = 0, eps = alpha.eps, 
                infinite = 0, attempts = alpha.attempts, 
                target.l = 0.5, target.u = 0.8)
  rho   <- list(cur = rho.init, att = 0, acc = 0, eps = rho.eps, 
                mn = logrho.mn, sd = logrho.sd, attempts = rho.attempts,
                target.l = 0.25, target.u = 0.5)
  
  # get calculated list
  calc <- list()
  calc$x.beta <- getXBeta(y = data$y, x = data$x, beta = beta$cur)
  calc$z      <- getZ(xi = xi$cur, x.beta = calc$x.beta, thresh = others$thresh)
  calc$lz     <- log(calc$z)
  calc$w      <- getW(rho = rho$cur, dw2 = others$dw2, a.cutoff = others$a.cutoff)
  calc$lw     <- log(calc$w)
  calc$w.star <- getWStarIDs(alpha = alpha$cur, w = calc$w, IDs = others$IDs)
  calc$aw     <- getAW(a = a$cur, w.star = calc$w.star)
  calc$theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
  
  # the stepsize for a and alpha should be standardized to account for the 
  # different magnitudes in the gradient
  if (is.null(alpha$eps)) {
    temp.q <- c(as.vector(log(a$cur)), transform$inv.logit(alpha$cur))
    temp.grad <- neg_log_post_grad_a_alpha(
      q = temp.q, data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
      rho = rho, calc = calc, others = others
    )
    alpha$eps <- a$eps * mean(temp.grad[1:nkt]) / tail(temp.grad, 1)
  } 
  
  # storage for MCMC
  storage.beta  <- matrix(NA, iters, np)
  storage.xi    <- rep(NA, iters)
  storage.a     <- array(NA, dim=c(iters, nknots, nt))
  storage.b     <- array(NA, dim=c(iters, nknots, nt))
  storage.alpha <- rep(NA, iters)
  storage.rho   <- rep(NA, iters)
  storage.prob  <- array(NA, dim = c(iters, ns, nt))
    
  tic <- proc.time()[3]
  for (iter in 1:iters) {
    
    beta$att <- beta$att + 1
    MHout <- updateBeta(data = data, beta = beta, xi = xi, alpha = alpha, 
                        calc = calc, others = others)
    if (MHout$accept) {
      beta$acc     <- beta$acc + 1
      beta$cur     <- MHout$q
      calc$x.beta  <- getXBeta(y = data$y, x = data$x, beta = beta$cur)
      calc$z       <- getZ(xi = xi$cur, x.beta = calc$x.beta, 
                           thresh = others$thresh)
      calc$lz      <- log(calc$z)
      calc$theta   <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
    }
    
    if (!(xi$fix)) {
      xi$att <- xi$att + 1
      MHout <- updateXi(data = data, xi = xi, alpha = alpha, calc = calc, 
                        others = others)
      if (MHout$accept) {
        xi$acc     <- xi$acc + 1
        xi$cur     <- MHout$q
        calc$z     <- getZ(xi = xi$cur, x.beta = calc$x.beta, 
                          thresh = others$thresh)
        calc$lz    <- log(calc$z)
        calc$theta <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
      }
    }
    
    if (a.alpha.joint) {
      a$att <- a$att + 1
      alpha$att <- alpha$att + 1
      
      q <- c(as.vector(log(a$cur)), transform$logit(alpha$cur))
      epsilon <- c(rep(a$eps, nkt), alpha$eps)
      HMCout <- gevHMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q, 
                    epsilon = epsilon, L = 20, 
                    data = data, beta = beta, xi = xi, a = a, b = b, 
                    alpha = alpha, rho = rho, calc = calc, others = others, 
                    this.param = "a_alpha")
      if (HMCout$accept) {
        a$acc     <- a$acc + 1
        a$cur     <- matrix(exp(HMCout$q[1:nkt]), nknots, nt)
        alpha$acc <- alpha$acc + 1
        alpha$cur <- transform$inv.logit(tail(HMCout$q, 1))
        calc$w.star <- getWStarIDs(alpha = alpha$cur, w = calc$w, IDs = others$IDs)
        calc$aw     <- getAW(a = a$cur, w.star = calc$w.star)
        calc$theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
      }
      
      a$infinite <- a$infinite + HMCout$infinite[1]
      alpha$infinite <- alpha$infinite + HMCout$infinite[2]
      
    } else {  # update a and alpha separately
      # random effect
      a$att <- a$att + 1
      q <- log(a$cur)
      HMCout  <- gevHMC(neg_log_post_a, neg_log_post_grad_a, q, 
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
      
      # spatial dependence
      alpha$att <- alpha$att + 1
      q <- transform$logit(alpha$cur)
      HMCout  <- gevHMC(neg_log_post_alpha, neg_log_post_grad_alpha, q, 
                     epsilon = alpha$eps, L = 10, 
                     data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                     rho = rho, calc = calc, others = others, this.param = "alpha")
      if (HMCout$accept) {
        alpha$acc   <- alpha$acc + 1
        alpha$cur   <- transform$inv.logit(HMCout$q)
        calc$w.star <- getWStarIDs(alpha = alpha$cur, w = calc$w, IDs = others$IDs)
        calc$aw     <- getAW(a = a$cur, w.star = calc$w.star)
        calc$theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
      }
      alpha$infinite <- alpha$infinite + HMCout$infinite
    }
    
    # auxiliary variable
    q <- transform$logit(b$cur)
    b$att <- b$att + 1
    HMCout  <- gevHMC(neg_log_post_b, neg_log_post_grad_b, q, epsilon = b$eps, 
                   L = 10, 
                   data = data, beta = beta, xi = xi, a = a, b = b, alpha = alpha, 
                   rho = rho, calc = calc, others = others, this.para = "b")
    if (HMCout$accept) {
      b$acc <- b$acc + 1
      b$cur <- transform$inv.logit(HMCout$q)
    }
    b$infinite <- b$infinite + HMCout$infinite
  
    # kernel bandwidth
    rho$att <- rho$att + 1
    MHout <- updateRho(data = data, a = a, alpha = alpha, rho = rho, calc = calc, 
                       others = others)
    if (MHout$accept) {
      rho$acc     <- rho$acc + 1
      rho$cur     <- MHout$q
      calc$w      <- getW(rho = rho$cur, dw2 = others$dw2, 
                          a.cutoff = others$a.cutoff)
      calc$lw     <- log(calc$w)
      calc$w.star <- getWStarIDs(alpha = alpha$cur, w = calc$w, others$IDs)
      calc$aw     <- getAW(a = a$cur, w.star = calc$w.star)
      calc$theta  <- getTheta(alpha = alpha$cur, z = calc$z, aw = calc$aw)
    }
    
    # adjustments to eps
    if (iter < burn / 2) {
      eps.update <- epsUpdate(beta)
      beta$att   <- eps.update$att
      beta$acc   <- eps.update$acc
      beta$eps   <- eps.update$eps
      
      if (!xi$fix) {
        eps.update <- epsUpdate(xi)
        xi$att     <- eps.update$att
        xi$acc     <- eps.update$acc
        xi$eps     <- eps.update$eps
      }
      
      if (a$infinite > 10) {
        print("reducing a$eps")
        a$eps <- a$eps * 0.8
        a$infinite <- 0
      }
      eps.update <- epsUpdate(a)
      a$att      <- eps.update$att
      a$acc      <- eps.update$acc
      a$eps      <- eps.update$eps
      
      if (b$infinite > 10) {
        print("reducing b$eps")
        b$eps <- b$eps * 0.8
        b$infinite <- 0
      }
      eps.update <- epsUpdate(b)
      b$att      <- eps.update$att
      b$acc      <- eps.update$acc
      b$eps      <- eps.update$eps
      
      if (alpha$infinite > 10) {
        print("reducing alpha$eps")
        alpha$eps <- alpha$eps * 0.8
        alpha$infinite <- 0
      }
      eps.update <- epsUpdate(alpha)
      alpha$att  <- eps.update$att
      alpha$acc  <- eps.update$acc
      alpha$eps  <- eps.update$eps
      
      eps.update <- epsUpdate(rho)
      rho$att    <- eps.update$att
      rho$acc    <- eps.update$acc
      rho$eps    <- eps.update$eps
    }
    
    storage.beta[iter, ] <- beta$cur
    storage.xi[iter]     <- xi$cur
    storage.a[iter, , ]  <- a$cur
    storage.b[iter, , ]  <- b$cur
    storage.alpha[iter]  <- alpha$cur
    storage.rho[iter]    <- rho$cur
    storage.prob[iter, , ] <- 1 - exp(-calc$theta)
    
    if (iter %% update == 0 | iter %in% c(621, 622, 623, 624, 625, 626, 627, 628, 629)) {
      cat("      Iter", iter, "of", iters, "\n")
      if (iterplot) {
        start <- max(iter - 5000, 1)
        end   <- iter
        par(mfrow=c(4, 5))
        plot.idx <- seq(1, 16, by = 2)
        for (idx in plot.idx){
          plot(log(storage.a[start:end, idx, 1]), type = "l", 
               main = paste("a", idx), 
               xlab = round(a$acc / a$att, 2), ylab = round(a$eps, 4))
        }
        plot(log(storage.a[start:end, 81, 1]), type = "l",
             main = paste("a", idx),
             xlab = round(a$acc / a$att, 2), ylab = round(a$eps, 4))
        plot.idx <- seq(1, 16, by = 2)
        for (idx in plot.idx){
          plot(storage.b[start:end, idx, 1], type = "l", 
               main = paste("b", idx),
               xlab = round(b$acc / b$att, 2), ylab = round(b$eps, 4))
        }
        
        plot(storage.beta[start:end], type = "l", main = bquote(beta[0]),
             xlab = round(beta$acc / beta$att, 2), ylab = round(beta$eps, 4))
        plot(storage.alpha[start:end], type = "l", main = bquote(alpha),
             xlab = round(alpha$acc / alpha$att, 2), ylab = round(alpha$eps, 4))
        plot(storage.rho[start:end], type = "l", main = bquote(rho),
             xlab = round(rho$acc / rho$att, 2), ylab = round(rho$eps, 4))
      }
    }
  }
  toc <- proc.time()[3]
  
  return.iters <- (burn + 1):iters
  results <- list(beta = storage.beta[return.iters, , drop = F],
                  xi = storage.xi[return.iters],
                  a = storage.a[return.iters, , , drop = F],
                  b = storage.b[return.iters, , , drop = F],
                  alpha = storage.alpha[return.iters],
                  rho = storage.rho[return.iters],
                  a.cutoff = a.cutoff,
                  minutes = (toc - tic) / 60)
}