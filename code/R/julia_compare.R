ns <- 1000
nt <- 365
y <- matrix(rbinom(ns * nt, size=1, prob=0.1), ns, nt)
theta.star <- matrix(rchisq(ns * nt, 3), ns, nt)
z <- matrix(rchisq(ns * nt, 2), ns, nt)
alpha <- 0.5
system.time(logLikeY(y, theta.star, alpha, z))