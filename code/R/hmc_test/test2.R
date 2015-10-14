rm(list=ls())
source("../hmc_aux.R", chdir = TRUE)

openblas.set.num.threads(2)

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 2000
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots  <- nknotsx * nknotsy
nkt     <- nknots * nt
rho.t   <- 0.25
alpha.t <- 0.2
s <- cbind(runif(ns), runif(ns))
knots <- expand.grid(seq(0, 1, length = nknotsx), seq(0, 1, length = nknotsy))
dw2 <- rdist(s, knots)
w.t <- stdW(makeW(dw2 = dw2, rho = rho.t))
z.t <- matrix(rgev(n = ns * nt, 1, 1, 1), ns, nt)
a.t <- matrix(rPS(nknots * nt, alpha = alpha.t), nknots, nt)
wz.t <- getwzCPP(z = z.t, w = w.t)
theta.t <- getThetaCPP(wz= wz.t, a_star = a.t^alpha.t, alpha = alpha.t)
y.t <- matrix(rbinom(ns * nt, size = 1, prob = -expm1(-theta.t)), ns, nt) 

niters <- 30000
storage.a <- array(NA, dim=c(niters, nknots, nt))
storage.b <- array(NA, dim=c(niters, nknots, nt))
storage.alpha <- rep(NA, niters)
q.a     <- matrix(log(10), nknots, nt)
q.b     <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = 0.5, wz = wz.t, 
               a = matrix(100, nknots, nt), b = matrix(0.5, nknots, nt))

# find vector for a and alpha
q.ps <- as.vector(c(q.a, 0))

set.seed(200)
tic <- proc.time()
for (i in 1:niters) {
  HMCout  <- HMC(neg_log_post_a_alpha, neg_log_post_grad_a_alpha, q.ps, 
                 epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    q.ps         <- HMCout$q
    others$a     <- matrix(exp(q.ps[1:nkt]), nknots, nt)
    others$alpha <- transform$inv.logit(tail(q.ps, 1))
  }

  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q.b, epsilon=0.005, L=10, others)
  if (HMCout$accept) {
    others$b <- transform$inv.logit(HMCout$q)
    q.b      <- HMCout$q
  }

  storage.a[i, , ] <- log(others$a)
  storage.b[i, , ] <- others$b
  storage.alpha[i] <- others$alpha
  if (i %% 500 == 0) {
    par(mfrow = c(3, 3))
    plot.idx <- c(4, 14, 24)
    print(paste("iter:", i, "of", niters, sep=" "))
    for (idx in plot.idx) {
      main.title <- paste("log(a[", idx, "])", sep = "")
      plot(storage.a[1:i, idx, ], type = "l", main = main.title)
    }
    for (idx in plot.idx) {
      main.title <- paste("b[", idx, "]", sep = "")
      plot(storage.b[1:i, idx, ], type = "l", main = main.title)
    }
    plot(storage.alpha[1:i], type = "l", main = "alpha")
  }
}
toc <- proc.time()
save.image(file = "alpha-2.RData")