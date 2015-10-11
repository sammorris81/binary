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
nknots <- nknotsx * nknotsy
rho.t <- 0.25
alpha.t <- 0.5
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
q.a <- matrix(log(100), nknots, nt)
q.b <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = q.b, a = q.a)

set.seed(200)
tic <- proc.time()
for (i in 1:niters) {
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q.a, epsilon=0.001, L=10, others)
  if (HMCout$accept) {
    q.a      <- HMCout$q
    others$a <- HMCout$q
  }
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q.b, epsilon=0.005, L=10, others)
  if (HMCout$accept) {
    others$b <- HMCout$q
    q.b <- HMCout$q
  }
  storage.a[i, , ] <- others$a
  storage.b[i, , ] <- others$b
  if (i %% 500 == 0) {
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc <- proc.time()
save.image(file = "alpha-5.RData")