source("auxfunctions.R")
source("updateModel.R")
set.seed(10)
ns   <- 200
nt   <- 1
s    <- cbind(runif(ns), runif(ns))
x <- array(1, dim=c(ns, nt, 3))
x[, , 2] <- s[, 1]
x[, , 3] <- s[, 2]

xi.t <- 0.1
beta.t <- c(5, 2, 0)

y <- rBinaryRareInd(x, beta=beta.t, xi=xi.t)


