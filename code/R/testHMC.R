rm(list=ls())
source("./hmc_aux.R")

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 400
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
rho.t <- 0.25
alpha.t <- 0.8
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

grad(func = neg_log_post_a, x = matrix(log(a.t)), others = others)
neg_log_post_grad_a(matrix(log(a.t)), others = others)

library(microbenchmark)
microbenchmark(grad(func = neg_log_post_a, x = matrix(log(a.t)), others = others), 
               neg_log_post_grad_a(matrix(log(a.t)), others = others))

for (i in 1:niters) {
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q.a, epsilon=0.001, L=2, others)
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
    par(mfrow=c(3, 4))
    plot(storage.a[1:i, 1, 1], type = "l")
    plot(storage.a[1:i, 3, 1], type = "l")
    plot(storage.a[1:i, 5, 1], type = "l")
    plot(storage.a[1:i, 7, 1], type = "l")
    plot(storage.a[1:i, 9, 1], type = "l")
    plot(storage.a[1:i, 11, 1], type = "l")
    plot(storage.a[1:i, 13, 1], type = "l")
    plot(storage.a[1:i, 15, 1], type = "l")
    plot(storage.a[1:i, 17, 1], type = "l")
    plot(storage.a[1:i, 21, 1], type = "l")
    plot(storage.a[1:i, 23, 1], type = "l")
    plot(storage.a[1:i, 25, 1], type = "l")
  }
}


plot(storage.b[1:i, 1, 1], type = "l")
plot(storage.b[1:i, 3, 1], type = "l")
plot(storage.b[1:i, 5, 1], type = "l")
plot(storage.b[1:i, 7, 1], type = "l")
plot(storage.b[1:i, 9, 1], type = "l")
plot(storage.b[1:i, 11, 1], type = "l")
plot(storage.b[1:i, 13, 1], type = "l")
plot(storage.b[1:i, 15, 1], type = "l")
plot(storage.b[1:i, 17, 1], type = "l")
plot(storage.b[1:i, 21, 1], type = "l")
plot(storage.b[1:i, 23, 1], type = "l")
plot(storage.b[1:i, 25, 1], type = "l")
# for running through HMC.R line by line
others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = q.b, a = q.a)
U <- neg_log_post_a
grad_U <- neg_log_post_grad_a
current_q <- q.a
epsilon=0.05
L=10

b.t <- matrix(runif(nknots * nt), nknots, nt)
q <- rbind(log(a.t), transform$logit(b.t))

others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = b.t, a = a.t)
neg_log_post(q, others)

q.a <- log(a.t)
neg_log_post_a(q = q.a, others)
neg_log_post_grad_a(q.a, others)

q.b <- transform$logit(b.t)
neg_log_post_b(q.b, others)
neg_log_post_grad_b(q.b, others)


f1 <- function(q, others) {
  # extract from the list
  y <- others$y
  alpha <- others$alpha
  
  nt <- ncol(y)
  nknots <- length(q)
  a  <- exp(q)
  b  <- transform$inv.logit(others$b)
  
  alpha1m <- 1 - alpha
  
  # start with the log prior
  # Remember: q = log(a)
  lc <- logc(b = b, alpha = alpha)
  ll <- sum(-alpha / alpha1m * q - exp(lc) * a^(-alpha / alpha1m))
  
  return(-ll)
}

gradf1 <- function(q, others) {
  # extract from the list
  y <- others$y
  alpha <- others$alpha
  
  nt <- ncol(y)
  nknots <- length(q)
  a  <- exp(q)
  b  <- transform$inv.logit(others$b)
  
  alpha1m <- 1 - alpha
  lc <- logc(b = b, alpha = alpha)
  grad <- -alpha / (alpha1m) * (1 + exp(lc) * a^(-alpha / alpha1m))
  
  return(-grad)
}

gradf1(matrix(log(a.t)), others) / grad(func = f1, x = matrix(log(a.t)), others = others)
gradf1(matrix(log(a.t)), others)
grad(func = f1, x = matrix(log(a.t)), others = others)
