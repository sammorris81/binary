library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)
# library(microbenchmark)
library(mvtnorm)
library(Rcpp)
library(numDeriv)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp(file = "../../code/R/pairwise.cpp")

source("../../code/R/auxfunctions.R", chdir = TRUE)
source("../../code/R/updateModel.R")
source("../../code/R/mcmc.R")
source("../../code/R/probit.R", chdir=T)

bs5.1.20 <- matrix(NA, 4, 10)

done <- c(1, 5, 9, 13, 17)
for (i in 1:10) {
  set <- (i - 1) * 2 + 1
  if (set %in% done) {
    file <- paste("sim-results/sim-knots-", set, "-5-1.RData", sep = "")
    print(paste("start: set", set))
    load(file)
    bs5.1.20[1, i] <- BrierScore(post.prob.gev, y.i.p)
    bs5.1.20[2, i] <- BrierScore(post.prob.log, y.i.p)
    bs5.1.20[3, i] <- BrierScore(post.prob.pro, y.i.p)
    file <- paste("sim-results/sim-knots-", set, "-5-1-r.RData", sep = "")
    load(file)
    bs5.1.20[4, i] <- BrierScore(post.prob.gev, y.i.p)
  }
}


done <- c(1, 6, 11, 16)
bs5.2.20 <- matrix(NA, 4, 10)
for (i in 1:10) {
  set <- (i - 1) * 2 + 1
  if (set %in% done) {
    file <- paste("sim-results/sim-knots-", set, "-5-2.RData", sep = "")
    print(paste("start: set", set))
    load(file)
    bs5.2.20[1, i] <- BrierScore(post.prob.gev, y.i.p)
    bs5.2.20[2, i] <- BrierScore(post.prob.log, y.i.p)
    bs5.2.20[3, i] <- BrierScore(post.prob.pro, y.i.p)
    file <- paste("sim-results/sim-knots-", set, "-5-2-r.RData", sep = "")
    load(file)
    bs5.2.20[4, i] <- BrierScore(post.prob.gev, y.i.p)
  }
}

bs5.2.10 <- matrix(NA, 4, 10)
idx <- 1
for (i in 1:10) {
  set <- i * 2
  if (set %in% done) {
    file <- paste("sim-results/sim-knots-", set, "-5-2.RData", sep = "")
    print(paste("start: set", set))
    load(file)
    bs5.2.10[1, i] <- BrierScore(post.prob.gev, y.i.p)
    bs5.2.10[2, i] <- BrierScore(post.prob.log, y.i.p)
    bs5.2.10[3, i] <- BrierScore(post.prob.pro, y.i.p)
    file <- paste("sim-results/sim-knots-", set, "-5-2-r.RData", sep = "")
    load(file)
    bs5.2.10[4, i] <- BrierScore(post.prob.gev, y.i.p)
  }
}

rownames(bs5.1.20) <- rownames(bs5.2.20) <- rownames(bs5.2.10) <- c("gev", "log", "pro", "gev2")
bs5.1.20 <- cbind(bs5.1.20, rowMeans(bs5.1.20, na.rm = TRUE))
bs5.2.20 <- cbind(bs5.2.20, rowMeans(bs5.2.20, na.rm = TRUE))
bs5.2.10 <- cbind(bs5.2.10, rowMeans(bs5.2.10, na.rm = TRUE))

round(bs5.1.20 * 100, 2)
round(bs5.2.20 * 100, 2)
round(bs5.2.10 * 100, 2)


# > round(bs5.1.20 * 100, 2)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# gev  4.11   NA 4.69   NA 3.86   NA 4.42   NA 4.93    NA  4.40
# log  3.13   NA 4.31   NA 3.86   NA 4.30   NA 4.92    NA  4.10
# pro  3.08   NA 4.36   NA 3.05   NA 4.50   NA 4.72    NA  3.94
# gev2   NA   NA 4.47   NA 3.41   NA 4.32   NA   NA    NA  4.07
# 
# 
# > round(bs5.2.20 * 100, 2)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# gev  4.30   NA   NA   NA   NA 4.22   NA   NA   NA    NA  4.26
# log  3.38   NA   NA   NA   NA 3.22   NA   NA   NA    NA  3.30
# pro  2.87   NA   NA   NA   NA 1.84   NA   NA   NA    NA  2.36
# gev2 3.99   NA   NA   NA   NA 3.08   NA   NA   NA    NA  3.54
# 
# > round(bs5.2.10 * 100, 2)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# gev    NA   NA 2.45   NA   NA   NA   NA 3.37   NA    NA  2.91
# log    NA   NA 3.15   NA   NA   NA   NA 1.29   NA    NA  2.22
# pro    NA   NA 1.08   NA   NA   NA   NA 1.20   NA    NA  1.14
# gev2   NA   NA 2.25   NA   NA   NA   NA 3.31   NA    NA  2.78