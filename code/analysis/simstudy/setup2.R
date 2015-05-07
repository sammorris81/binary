#########################################################################
# A small-scale simulation study to determine when spatial GEV link
# performs better than
#
# data settings:
#   All: s in [0, 6] x [0, 6]
#   1: GEV link
#      a: alpha = 0.3, 144 knots, 1% rareness, bw = 0.5, xi = 0.25
#      b: alpha = 0.7, 144 knots, 1% rareness, bw = 0.5, xi = 0.25
#      c: alpha = 0.3, 144 knots, 5% rareness, bw = 0.5, xi = 0.25
#      d: alpha = 0.7, 144 knots, 5% rareness, bw = 0.5, xi = 0.25
#   2: Logit link
#      a: rho = 3, 144 knots, 1% rareness
#      b: rho = 1, 144 knots, 1% rareness
#      c: rho = 3, 144 knots, 5% rareness
#      d: rho = 1, 144 knots, 5% rareness
#
# methods:
#   1: Spatial logit   (spbayes)
#   2: Spatial probit
#   3: Spatial GEV (our method)
#   4: Spatial GEV (our method with alpha and rho fixed)
#########################################################################

rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)

source("../../R/auxfunctions.R", chdir=T)
source("../../R/updateModel.R", chdir=T)
source("../../R/mcmc.R", chdir=T)
source("../../R/probit.R", chdir=T)

ns        <- 4000
nt        <- 1
np        <- 1
nsets     <- 10
nlinks    <- 2   # using Max-stable and logit
nsettings <- 4   # a, b, c, d from above

knots.1 <- knots.2 <- seq(from = 0.2, to = 5.8, length = 12)
knots <- as.matrix(expand.grid(knots.1, knots.2))
nknots <- nrow(knots)

X <- matrix(1, nrow=ns, ncol=np)

y         <- array(NA, dim=c(ns, nt, nsets, nsettings * nlinks))
a.t       <- array(NA, dim=c(nknots, nt, nsets, nsettings))
w.tilde.t <- array(NA, dim=c(ns, nt, nsets, nsettings))

set.seed <- 7483 # site
s <- cbind(runif(ns, 0, 6), runif(ns, 0, 6))

# max-stable % rareness is set with prob.ms
prob.ms  <- c(0.01, 0.01, 0.05, 0.05)
int.ms   <- matrix(NA, nsets, nsettings)
alpha.ms <- c(0.3, 0.7, 0.3, 0.7)
xi.ms    <- 0.25
rho.ms   <- 6 / 11

# to generate logistic data with around 1% rareness, need intercept = -log(99)
# to generate logistic data with around 5% rareness, need intercept = -log(19)
int.logit   <- c(-log(99), -log(99), -log(19), -log(19))
rho.logit   <- c(3, 1, 3, 1)
sigsq.logit <- 1
nu.logit    <- 0.5

set.seed <- 3181 # data
for (i in 1:nsets) {
  for (j in 1:nlinks) {
    for (k in 1:nsettings) {
      setting <- (j - 1) * nsettings + k
      if (j == 1) {  # max-stable
        data <- rRareBinarySpat(x = X, s = s, knots = knots, beta = 0, xi.ms,
                                alpha = alpha.ms[k], rho.ms,
                                prob.success = prob.ms[k])
        a.t[, , i, k] <- data$a
        int.ms[i, k] <- data$thresh
      } else {  # logit
        data <- rLogitSpat(x = X, s = s, knots = knots, beta = int.logit[k],
                           rho = rho.logit[k], sigma.sq = sigsq.logit,
                           nu = nu.logit)
        w.tilde.t[, , i, k] <- data$w.tilde
      }
      y[, , i, setting] <- data$y
    }
  }
  cat("Set", i, "finished \n")
}

save.image(file = "simdata2.RData")

# look at a few of the datasets
# y[, , set, setting]
load("simdata2.RData")
ntrain <- 3000
ntest  <- 1000
obs <- c(rep(T, ntrain), rep(F, ntest))

par(mfrow=c(2, 2))
these.train <- which(y[obs, 1, 1, 1] == 1)
these.test  <- which(y[!obs, 1, 1, 1] == 1)
s.o <- s[obs, ]
s.p <- s[!obs, ]
plot(s.o[these.train, ], pch=22, bg = "dodgerblue1", col = "dodgerblue4",
     xlab = "", ylab = "", xlim = c(0, 6), ylim = c(0, 6),
     main = bquote(paste(pi, "= 0.01, ", alpha, " = 0.3", sep="")))
points(s.p[these.test, ], pch=22, bg="firebrick1", col="firebrick4")
points(knots)

these.train <- which(y[obs, 1, 5, 1] == 1)
these.test  <- which(y[!obs, 1, 5, 1] == 1)
s.o <- s[obs, ]
s.p <- s[!obs, ]
plot(s.o[these.train, ], pch=22, bg = "dodgerblue1", col = "dodgerblue4",
     xlab = "", ylab = "", xlim = c(0, 6), ylim = c(0, 6),
     main = bquote(paste(pi, "= 0.01, ", alpha, " = 0.3", sep="")))
points(s.p[these.test, ], pch=22, bg="firebrick1", col="firebrick4")
points(knots)

these.train <- which(y[obs, 1, 7, 1] == 1)
these.test  <- which(y[!obs, 1, 7, 1] == 1)
s.o <- s[obs, ]
s.p <- s[!obs, ]
plot(s.o[these.train, ], pch=22, bg = "dodgerblue1", col = "dodgerblue4",
     xlab = "", ylab = "", xlim = c(0, 6), ylim = c(0, 6),
     main = bquote(paste(pi, "= 0.01, ", alpha, " = 0.3", sep="")))
points(s.p[these.test, ], pch=22, bg="firebrick1", col="firebrick4")
points(knots)

these.train <- which(y[obs, 1, 10, 1] == 1)
these.test  <- which(y[!obs, 1, 10, 1] == 1)
s.o <- s[obs, ]
s.p <- s[!obs, ]
plot(s.o[these.train, ], pch=22, bg = "dodgerblue1", col = "dodgerblue4",
     xlab = "", ylab = "", xlim = c(0, 6), ylim = c(0, 6),
     main = bquote(paste(pi, "= 0.01, ", alpha, " = 0.3", sep="")))
points(s.p[these.test, ], pch=22, bg="firebrick1", col="firebrick4")
points(knots)


# 5% rareness
these.train <- which(y[obs, 1, 1, 3] == 1)
these.test  <- which(y[!obs, 1, 1, 3] == 1)
s.o <- s[obs, ]
s.p <- s[!obs, ]
plot(s.o[these.train, ], pch=22, bg = "dodgerblue1", col = "dodgerblue4",
     xlab = "", ylab = "", xlim = c(0, 6), ylim = c(0, 6),
     main = bquote(paste(pi, "= 0.05, ", alpha, " = 0.3", sep="")))
points(s.p[these.test, ], pch=22, bg="firebrick1", col="firebrick4")
points(knots)

these.train <- which(y[obs, 1, 5, 3] == 1)
these.test  <- which(y[!obs, 1, 5, 3] == 1)
s.o <- s[obs, ]
s.p <- s[!obs, ]
plot(s.o[these.train, ], pch=22, bg = "dodgerblue1", col = "dodgerblue4",
     xlab = "", ylab = "", xlim = c(0, 6), ylim = c(0, 6),
     main = bquote(paste(pi, "= 0.05, ", alpha, " = 0.3", sep="")))
points(s.p[these.test, ], pch=22, bg="firebrick1", col="firebrick4")
points(knots)

these.train <- which(y[obs, 1, 7, 3] == 1)
these.test  <- which(y[!obs, 1, 7, 3] == 1)
s.o <- s[obs, ]
s.p <- s[!obs, ]
plot(s.o[these.train, ], pch=22, bg = "dodgerblue1", col = "dodgerblue4",
     xlab = "", ylab = "", xlim = c(0, 6), ylim = c(0, 6),
     main = bquote(paste(pi, "= 0.05, ", alpha, " = 0.3", sep="")))
points(s.p[these.test, ], pch=22, bg="firebrick1", col="firebrick4")
points(knots)

these.train <- which(y[obs, 1, 10, 3] == 1)
these.test  <- which(y[!obs, 1, 10, 3] == 1)
s.o <- s[obs, ]
s.p <- s[!obs, ]
plot(s.o[these.train, ], pch=22, bg = "dodgerblue1", col = "dodgerblue4",
     xlab = "", ylab = "", xlim = c(0, 6), ylim = c(0, 6),
     main = bquote(paste(pi, "= 0.05, ", alpha, " = 0.3", sep="")))
points(s.p[these.test, ], pch=22, bg="firebrick1", col="firebrick4")
points(knots)

# get extremal coefficient as function of distance
rm(list=ls())
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)

source("../../R/auxfunctions.R", chdir=T)
source("../../R/updateModel.R", chdir=T)
source("../../R/mcmc.R", chdir=T)
source("../../R/probit.R", chdir=T)

ns        <- 4000
nt        <- 1
np        <- 1
nsets     <- 200
nsettings <- 4   # a, b, c, d from above

knots.1 <- knots.2 <- seq(from = 0.2, to = 5.8, length = 12)
knots <- as.matrix(expand.grid(knots.1, knots.2))
nknots <- nrow(knots)

X <- matrix(1, nrow=ns, ncol=np)

y         <- array(NA, dim=c(ns, nt, nsets, nsettings))
a.t       <- array(NA, dim=c(nknots, nt, nsets, nsettings))
w.tilde.t <- array(NA, dim=c(ns, nt, nsets, nsettings))

set.seed <- 7483 # site
s <- cbind(runif(ns, 0, 6), runif(ns, 0, 6))

d <- rdist(s)
diag(d) <- 0
bins <- seq(0, max(d), length=100)
xplot <- (bins[-1] + bins[-100]) / 2
acc <- att <- array(0, dim=c(length(bins) - 1, nsets, 4))

# max-stable % rareness is set with prob.ms
prob.ms  <- c(0.01, 0.01, 0.05, 0.05)
int.ms   <- matrix(NA, nsets, nsettings)
alpha.ms <- c(0.3, 0.7, 0.3, 0.7)
xi.ms    <- 0.25
rho.ms   <- 6 / 11

library(doMC)
registerDoMC(4)

set.seed <- 3181 # data
for (j in 1:nsettings) {
  for (i in 1:nsets) {
    data <- rRareBinarySpat(x = X, s = s, knots = knots, beta = 0, xi.ms,
                            alpha = alpha.ms[j], rho.ms,
                            prob.success = prob.ms[j])
    a.t[, , i, j] <- data$a
    int.ms[i, j] <- data$thresh
    y[, , i, j] <- data$y
  }
  cat("Seting", j, "finished \n")
}

results <- foreach (setting = 1:4) %dopar% {
  acc <- att <- matrix(0, nrow=length(bins) - 1, ncol=nsets)
# acc <- att <- array(0, dim=c(length(bins) - 1, nsets, 4))
# for (setting in 1:2) {
  for (set in 1:nsets) {
    these <- which(y[, , set, setting] == 1)
    for (j in these) {
      for (i in 1:ns) {
        if (i != j) {
          dij <- d[i, j]  # get the distance
          this.bin <- max(which(dij > bins))  # figure out which bin it belongs to
          att[this.bin, set] <- att[this.bin, set] + 1
          # att[this.bin, set, setting] <- att[this.bin, set, setting] + 1
          acc[this.bin, set] <- acc[this.bin, set] + y[i, , set, setting]
          # acc[this.bin, set, setting] <- acc[this.bin, set, setting] +
                                         # y[i, , set, setting]
        }
      }
      # print(paste("Setting ", setting, ", set ", set, ", site ", j, sep=""))
    }
    # print(paste("Setting ", setting, ", set ", set, sep=""))
  }
  results <- list(acc=acc, att=att)
  return(results)
}

save(results, file="simulatedchi.RData")
load("simulatedchi.RData")

acc <- att <- array(0, dim=c(length(bins) - 1, nsets, 4))
for (setting in 1:4) {
  acc[, , setting] <- results[[setting]]$acc
  att[, , setting] <- results[[setting]]$att
}


theta <- 2 - acc / att
par(mfrow=c(2, 2))
for (setting in 1:1) {
  plot(xplot, theta[, 1, setting], type = "l",
       main = bquote(paste(vartheta, "(h)")),
       xlab = "h", ylab = paste("setting", setting), ylim = c(1, 2))
  for (line in 2:200) {
    lines(xplot, theta[, line, 1])
  }
}

mean.theta <- apply(theta, c(1, 3), mean, na.rm=TRUE)
for (setting in 1:1) {
  plot(xplot, mean.theta[, setting], type = "l",
       main = bquote(paste("mean ", vartheta, "(h)")),
       xlab = "h", ylab = paste("setting", setting), ylim = c(1, 2))
}
