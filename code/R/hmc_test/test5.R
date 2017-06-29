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
storage.alpha <- rep(NA, niters)
q.a     <- matrix(log(100), nknots, nt)
q.b     <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = 0.5, wz = wz.t, 
               a = matrix(100, nknots, nt), b = matrix(0.5, nknots, nt))

# find a good starting point for alpha
q.alpha.temp <- seq(-10, 10, 0.01)
ll.temp <- rep(Inf, length(q.alpha.temp))
for (i in 1:length(q.alpha.temp)) {
  ll.temp[i] <- neg_log_post_alpha(q = q.alpha.temp[i], others = others)
}

# pick the starting value that minimizes the neg log likelihood
q.alpha <- q.alpha.temp[which(ll.temp == min(ll.temp))]
others$alpha <- transform$inv.logit(q.alpha)

set.seed(200)
tic <- proc.time()
for (i in 1:niters) {
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q.a, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    others$a <- exp(HMCout$q)
    q.a      <- HMCout$q
  }
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q.b, epsilon=0.005, L=10, others)
  if (HMCout$accept) {
    others$b <- transform$inv.logit(HMCout$q)
    q.b      <- HMCout$q
  }
  HMCout  <- HMC(neg_log_post_alpha, neg_log_post_grad_alpha, q.alpha, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    others$alpha <- transform$inv.logit(HMCout$q)
    q.alpha      <- HMCout$q
  }
  storage.a[i, , ] <- others$a
  storage.b[i, , ] <- others$b
  storage.alpha[i] <- others$alpha
  if (i %% 500 == 0) {
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc <- proc.time()
save.image(file = "alpha-5.RData")