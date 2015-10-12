rm(list=ls())
source("./hmc_aux.R")

# Test out the functions
library(fields)
library(evd)
set.seed(200)
ns <- 1500
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
storage.a.8 <- array(NA, dim=c(niters, nknots, nt))
storage.b.8 <- array(NA, dim=c(niters, nknots, nt))
q.a <- matrix(log(10), nknots, nt)
q.b <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = q.b, a = q.a)

set.seed(200)
tic.1 <- proc.time()
for (i in 1:niters) {
  HMCout  <- HMC(neg_log_post_a, neg_log_post_grad_a, q.a, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    q.a      <- HMCout$q
    others$a <- HMCout$q
  }
  HMCout  <- HMC(neg_log_post_b, neg_log_post_grad_b, q.b, epsilon=0.01, L=10, others)
  if (HMCout$accept) {
    others$b <- HMCout$q
    q.b <- HMCout$q
  }
  storage.a.8[i, , ] <- others$a
  storage.b.8[i, , ] <- others$b
  if (i %% 100 == 0) {
    par(mfrow=c(3, 4))
    plot.idx <- seq(1, 23, by = 2)
    for (idx in plot.idx){
      plot(storage.b.8[1:i, idx, 1], type = "l",  main = round(log(a.t[idx, 1]), 2))
    }
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.1 <- proc.time()
save.image(file = "alpha8.RData")


set.seed(200)
ns <- 10
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
storage.a.5 <- array(NA, dim=c(niters, nknots, nt))
storage.b.5 <- array(NA, dim=c(niters, nknots, nt))
q.a <- matrix(log(10), nknots, nt)
q.b <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = q.b, a = q.a)

set.seed(200)
tic.2 <- proc.time()
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
  storage.a.5[i, , ] <- others$a
  storage.b.5[i, , ] <- others$b
  if (i %% 500 == 0) {
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.2 <- proc.time()
save.image(file = "alpha5.RData")

set.seed(200)
ns <- 2000
nt <- 1
nknotsx <- 5
nknotsy <- 5
nknots <- nknotsx * nknotsy
rho.t <- 0.25
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
storage.a.2 <- array(NA, dim=c(niters, nknots, nt))
storage.b.2 <- array(NA, dim=c(niters, nknots, nt))
q.a <- matrix(log(100), nknots, nt)
q.b <- matrix(0, nknots, nt)
others <- list(y = y.t, alpha = alpha.t, wz = wz.t, b = q.b, a = q.a)

set.seed(200)
tic.3 <- proc.time()
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
  storage.a.2[i, , ] <- others$a
  storage.b.2[i, , ] <- others$b
  if (i %% 500 == 0) {
    print(paste("iter:", i, "of", niters, sep=" "))
  }
}
toc.3 <- proc.time()
save.image(file = "alpha2.RData")


par(mfrow=c(3, 4))
plot.idx <- seq(1, 23, by = 2)
for (i in plot.idx){
  plot(storage.a[, i, 1], type = "l",  main = round(log(a.t[1, 1]), 2))
}

for (i in plot.idx){
  plot(storage.b[, i, 1], type = "l")
}

toc - tic



# neg_log_post_a(q.a, others)
# 
# grad(func = neg_log_post_a, x = q.a, others = others)
# neg_log_post_grad_a(q.a, others = others)
# 
# grad(func = neg_log_post_b, x = q.b, others = others)
# neg_log_post_grad_b(q.b, others = others)

# library(microbenchmark)
# microbenchmark(grad(func = neg_log_post_a, x = matrix(log(a.t)), others = others), 
#                neg_log_post_grad_a(matrix(log(a.t)), others = others))

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
  alpha <- others$alpha
  
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

library(numDeriv)
gradf1(matrix(log(a.t)), others) / grad(func = f1, x = matrix(log(a.t)), others = others)
gradf1(matrix(log(a.t)), others)
grad(func = f1, x = matrix(log(a.t)), method = "complex", others = others)
log(a.t)

xplot <- seq(-2, 3, length = 1000)
yplot <- rep(0, length(xplot))
for (i in 1:length(xplot)) {
  q <- matrix(xplot[i], 1, 1)
  yplot[i] <- neg_log_post_a(q, others)
}
plot(xplot, yplot, type = "l")
neg_log_post_grad_a(q = log(a.t), others)
grad(func = neg_log_post_a, x = matrix(0, 1, 1), others = others)
neg_log_post_grad_a(q = matrix(0, 1, 1), others)

part1 <- -alpha.t / (1 - alpha.t)
part2 <- -(sin(alpha.t * pi * 0.5) / sin(pi * 0.5))^(1 / (1 - alpha.t)) * sin((1 - alpha.t) * pi * 0.5) / sin(alpha.t * pi * 0.5)
part2 <- part2 * alpha.t / (1- alpha.t)
part2 <- part2 * exp(2)^(-alpha.t / (1 - alpha.t))

# are the gradients reasonable
a.t[3] <- 1
neg_log_post_grad_a(q = log(a.t), others)

xplot <- seq(-4, 4, length = 1000)
yplot <- rep(0, length(xplot))
for (i in 1:length(yplot)) {
  q <- matrix(xplot[i], 1, 1)
  yplot[i] <- neg_log_post_b(q, others)
}
plot(xplot, yplot, type = "l")

which(yplot == min(yplot))
xplot[926]
grad(func = neg_log_post_b, x = matrix(3.4, 1, 1), others = others)
grad(func = neg_log_post_b, x = matrix(3.5, 1, 1), others = others)
neg_log_post_grad_b(matrix(3.4, 1, 1), others = others)


long <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  alpha1m <- 1 - alpha
  
  u <- (sin(alpha * pi * b) / sin(pi * b))
  grad <- u^(-1 / alpha1m) / sin(alpha1m * pi * b) * ((
    1 / (alpha1m * sin(alpha * pi * b)^2) * (
      u^(1 / alpha1m) * (
        alpha * pi *  cos(alpha * pi * b) / sin(pi * b) -
          sin(alpha * pi * b) * pi * cos(pi * b) / sin(pi * b)^2
        ) * sin(pi * b) * sin(alpha1m * pi * b)
      ) + 
      u^(1 / alpha1m) * alpha1m * pi * cos(alpha1m * pi * b) / sin(alpha * pi * b) - 
      u^(1 / alpha1m) * sin(alpha1m * pi * b) * alpha * pi * cos(alpha * pi * b) / sin(alpha * pi * b)^2
      ) * sin(alpha * pi * b))
  return (-grad)
}

short <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  alpha1m <- 1 - alpha
  
  u <- sin(alpha * pi * b) / sin(pi * b)
  
  grad <- alpha * pi / (alpha1m * tan(alpha * pi * b)) - u * pi * cos(pi * b) / (alpha1m * sin(alpha * pi * b)) +
    alpha1m * pi / tan(alpha1m * pi * b) - alpha * pi / tan(alpha * pi * b)
return (-grad)
}

shorter <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  alpha1m <- 1 - alpha
  
  apb <- alpha * pi * b
  pb  <- pi * b
  a1mpb <- alpha1m * pi * b
  
  u <- sin(apb) / sin(pb)
  
  grad <- alpha * pi / (alpha1m * tan(apb)) - u * pi * cos(pb) / (alpha1m * sin(apb)) +
    alpha1m * pi / tan(a1mpb) - alpha * pi / tan(apb)
  return (-grad)
}

long2 <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  
  # repeated quantities
  alpha1m <- 1 - alpha
  apb <- alpha * pi * b
  pb  <- pi * b
  a1mpb <- alpha1m * pi * b
  
  u <- (sin(apb) / sin(pb))
  
  grad <- u^(1 / alpha1m) * (alpha * pi * cos(apb) / sin(pb) - 
                       sin(apb) * pi * cos(pb) / sin(pb)^2) * 
    sin(pb) * sin(a1mpb) * a^(-alpha / alpha1m) / (alpha1m * sin(apb)^2)
  
  return(grad)
}

long3 <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  
  # repeated quantities
  alpha1m <- 1 - alpha
  apb <- alpha * pi * b
  pb  <- pi * b
  a1mpb <- alpha1m * pi * b
  
  u <- (sin(apb) / sin(pb))
  
  grad <- u^(1 / alpha1m) * alpha1m * pi * cos(a1mpb) * a^(-alpha / alpha1m) / sin(apb)
  return(grad)
}

long4 <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  
  # repeated quantities
  alpha1m <- 1 - alpha
  apb <- alpha * pi * b
  pb  <- pi * b
  a1mpb <- alpha1m * pi * b
  
  u <- (sin(apb) / sin(pb))
  
  grad <- -u^(1 / alpha1m) * sin(a1mpb) * a^(-alpha / alpha1m) * alpha * pi * cos(apb) / 
    sin(apb)^2
  
  return(grad)
}

short2 <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  
  # repeated quantities
  alpha1m <- 1 - alpha
  apb <- alpha * pi * b
  pb  <- pi * b
  a1mpb <- alpha1m * pi * b
  
  u <- (sin(apb) / sin(pb))
  
  grad <- u^(1 / alpha1m) * (alpha * pi * cos(apb) - u * pi * cos(pb)) * sin(a1mpb) * a^(-alpha / alpha1m) / (alpha1m * sin(apb)^2)
  return(grad)
}

short3 <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  
  # repeated quantities
  alpha1m <- 1 - alpha
  apb <- alpha * pi * b
  pb  <- pi * b
  a1mpb <- alpha1m * pi * b
  
  u <- (sin(apb) / sin(pb))
  
  grad <- u^(1 / alpha1m) * alpha1m * pi * cos(a1mpb) * a^(-alpha / alpha1m) / sin(apb)
  return(grad)
}

short4 <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  
  # repeated quantities
  alpha1m <- 1 - alpha
  apb <- alpha * pi * b
  pb  <- pi * b
  a1mpb <- alpha1m * pi * b
  
  u <- (sin(apb) / sin(pb))
  grad <- -u^(1 / alpha1m) * sin(a1mpb) * a^(-alpha / alpha1m) * alpha * pi / (sin(apb) * tan(apb))
  
  return(grad)
}

short.2 <- function(q, others) {
  a <- exp(others$a)
  alpha <- others$alpha
  b <- transform$inv.logit(q)
  
  # repeated quantities
  alpha1m <- 1 - alpha
  apb <- alpha * pi * b
  pb  <- pi * b
  a1mpb <- alpha1m * pi * b
  
  u <- (sin(apb) / sin(pb))
  
  grad <- u^(1 / alpha1m) * a^(-alpha / alpha1m) * (
    (alpha * pi * cos(apb) - u * pi * cos(pb)) * sin(a1mpb) / (alpha1m * sin(apb)) +
      alpha1m * pi * cos(a1mpb) - alpha * pi * sin(a1mpb) / tan(apb)
  ) / sin(apb)
  
  return(grad)
}

grad(func = neg_log_post_b, x = matrix(0, 1, 1), others = others)
neg_log_post_grad_b(q = matrix(0, 1, 1), others = others)
 
grad(func = neg_log_post_b, x = matrix(3, 1, 1), others = others)
neg_log_post_grad_b(q = matrix(3, 1, 1), others = others)

grad(func = neg_log_post_b, x = matrix(3.4, 1, 1), others = others)
neg_log_post_grad_b(q = matrix(3.4, 1, 1), others = others)

grad(func = neg_log_post_b, x = matrix(3.5, 1, 1), others = others)
neg_log_post_grad_b(q = matrix(3.5, 1, 1), others = others)

grad(func = neg_log_post_b, x = matrix(3.41, 1, 1), others = others)
neg_log_post_grad_b(q = matrix(3.41, 1, 1), others = others)

grad(func = neg_log_post_b, x = matrix(3.4, 1, 1), others = others)
neg_log_post_grad_b(q = matrix(3.4, 1, 1), others = others)
