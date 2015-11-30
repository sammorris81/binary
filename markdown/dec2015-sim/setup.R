#### Exploring a couple of hotspots
rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)
# library(microbenchmark)  # comment out for beowulf
library(mvtnorm)
library(Rcpp)
library(numDeriv)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

source("../../code/R/spatial_gev.R", chdir = TRUE)
source("../../code/R/spatial_logit.R", chdir = TRUE)
source("../../code/R/spatial_probit.R", chdir = TRUE)

set.seed(7483)  # site
ns <- c(650, 1300, 650, 1300, 650, 1300)  # 500 train and 1000 train
nsettings <- length(ns)  # storing y in a list

### gev settings

gev.alpha <- 0.3
gev.rho   <- 0.05
gev.xi    <- 0
gev.prob  <- 0.05

### logit settings

log.gamma <- 0.99  # proportion of variability due to spatial
log.var   <- 1
log.rho   <- 0.025
log.prob  <- 0.03

### hotspot settings. generates around 5% 1s.

nhotspots <- 9
p <- 0.60  # P(Y=1|hot spot)
q <- 0.01  # P(Y=1|background)
r <- 0.05  # Hot spot radius

simdata <- vector(mode = "list", length = nsettings)

nsets <- 100
set.seed(3282)  # data

for (setting in 1:nsettings) {
  simdata[[setting]]$y <- matrix(data = NA, nrow = ns[setting], ncol = nsets)
  simdata[[setting]]$s <- array(runif(2 * nsets * ns[setting]), 
                                dim = c(ns[setting], 2, nsets))
  simdata[[setting]]$x <- matrix(1, ns[setting], 1)
  simdata[[setting]]$hotspots <- vector(mode = "list", length = nsets)
  
  for (set in 1:nsets) {
    ### GEV generation
#     if (setting == 1 | setting == 2) {
#       data <- rRareBinarySpat(x = simdata[[setting]]$x, 
#                               s = simdata[[setting]]$s[, , set], 
#                               knots = simdata[[setting]]$s[, , set], 
#                               beta = 0, xi = gev.xi, alpha = gev.alpha, 
#                               rho = gev.rho, prob.success = gev.prob)
#       
#       simdata[[setting]]$y[, set] <- data$y
#     } 
    
    ### logit generation
    if (setting == 3 | setting == 4) {
      d <- as.matrix(rdist(simdata[[setting]]$s[, , set]))
      diag(d) <- 0
      Sigma <- simple.cov.sp(D=d, sp.type="exponential", 
                             sp.par=c(log.var, log.rho), error.var=0, 
                             finescale.var=0)
      data <- transform$logit(log.prob) + t(chol(Sigma)) %*% rnorm(ns[setting])
      data <- rbinom(n = ns[setting], size = 1, prob = transform$inv.logit(data))
      simdata[[setting]]$y[, set] <- data
    }
    
    ### hotspot generation
    if (setting == 5 | setting == 6) {
      k  <- rpois(1, nhotspots) + 1
      c1 <- runif(k)
      c2 <- runif(k)
      
      hot <- rep(FALSE, ns[setting])
      for (j in 1:k) {
        d   <- sqrt((simdata[[setting]]$s[, 1, set] - c1[j])^2 + 
                    (simdata[[setting]]$s[, 2, set] - c2[j])^2) 
        hot <- ifelse(d < r, TRUE, hot)
      }
      simdata[[setting]]$y[, set] <- rbinom(ns[setting], 1, ifelse(hot, p, q))
      simdata[[setting]]$hotspots[[set]] <- cbind(c1, c2)
    }
    
    if (set %% 20 == 0) {
      print(paste("  Dataset ", set, " finished", sep = ""))
    }
  }
  print(paste("Setting ", setting, " finished", sep = ""))
}

# quartz(width = 10, height = 7)
par(mfrow = c(3, 2))

setting <- 1
set     <- 1
# plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
#      main = "Setting 1: 5%, ns = 1000")
plot(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] != 1), , set], pch = 21, 
     col = "dodgerblue4", bg = "dodgerblue1", 
     xlab = "", ylab = "", 
     main = paste("Setting ", setting, ": GEV - 10%, ns = ", ns[setting], sep = ""))
points(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] == 1), , set], pch = 21, 
       col = "firebrick4", bg = "firebrick1")

setting <- 2
set     <- 1
# plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
#      main = "Setting 1: 5%, ns = 1000")
plot(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] != 1), , set], pch = 21, 
     col = "dodgerblue4", bg = "dodgerblue1", 
     xlab = "", ylab = "", 
     main = paste("Setting ", setting, ": GEV - 10%, ns = ", ns[setting], sep = ""))
points(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] == 1), , set], pch = 21, 
       col = "firebrick4", bg = "firebrick1")

setting <- 3
set     <- 3
# plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
#      main = "Setting 1: 5%, ns = 1000")
plot(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] != 1), , set], pch = 21, 
     col = "dodgerblue4", bg = "dodgerblue1", 
     xlab = "", ylab = "", 
     main = paste("Setting ", setting, ": Logit - 10%, ns = ", ns[setting], sep = ""))
points(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] == 1), , set], pch = 21, 
       col = "firebrick4", bg = "firebrick1")

setting <- 4
set     <- 3
# plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
#      main = "Setting 1: 5%, ns = 1000")
plot(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] != 1), , set], pch = 21, 
     col = "dodgerblue4", bg = "dodgerblue1", 
     xlab = "", ylab = "", 
     main = paste("Setting ", setting, ": Logit - 10%, ns = ", ns[setting], sep = ""))
points(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] == 1), , set], pch = 21, 
       col = "firebrick4", bg = "firebrick1")

setting <- 5
set     <- 1
# plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
#      main = "Setting 1: 5%, ns = 1000")
plot(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] != 1), , set], pch = 21, 
     col = "dodgerblue4", bg = "dodgerblue1", 
     xlab = "", ylab = "", 
     main = paste("Setting ", setting, ": Hotspot - 10%, ns = ", ns[setting], sep = ""))
points(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] == 1), , set], pch = 21, 
       col = "firebrick4", bg = "firebrick1")

setting <- 6
set     <- 1
# plot(knots, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "", 
#      main = "Setting 1: 5%, ns = 1000")
plot(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] != 1), , set], pch = 21, 
     col = "dodgerblue4", bg = "dodgerblue1", 
     xlab = "", ylab = "", 
     main = paste("Setting ", setting, ": Hotspot - 10%, ns = ", ns[setting], sep = ""))
points(simdata[[setting]]$s[which(simdata[[setting]]$y[, set] == 1), , set], pch = 21, 
       col = "firebrick4", bg = "firebrick1")

mean(simdata[[1]]$y)
mean(simdata[[2]]$y)
mean(simdata[[3]]$y)
mean(simdata[[4]]$y)

save(simdata, gev.rho, gev.prob, log.rho, log.prob, file = "simdata10.RData")

# # setting trial 1
# nhotspots <- c(5, 5, 3, 3)
# knots   <- expand.grid(seq(0, 1, length=21), seq(0, 1, length=21))
# p <- 0.95  # P(Y=1|hot spot)
# q <- 0.01  # P(Y=1|background)
# r <- 0.05  # Hot spot radius

# # setting trial 2
# nhotspots <- c(7, 7, 3, 3)
# knots   <- expand.grid(seq(0, 1, length=21), seq(0, 1, length=21))
# p <- 0.65  # P(Y=1|hot spot)
# q <- 0.01  # P(Y=1|background)
# r <- 0.05  # Hot spot radius

# # setting trial 3
# nhotspots <- c(5, 5, 2, 2)
# knots   <- expand.grid(seq(0, 1, length=13), seq(0, 1, length=13))
# p <- 0.400  # P(Y=1|hot spot)
# q <- 0.005  # P(Y=1|background)
# r <- 0.083  # Hot spot radius
#        [,1]   [,2]   [,3]
# [1,] 0.0427 0.0418 0.0519
# [2,] 0.0377 0.0364 0.0462
# [3,] 0.0208 0.0205 0.0257
# [4,] 0.0234 0.0234 0.0285


library(fields)
library(SpatialTools)
log.var  <- c(1 ,3, 5, 7, 9, 11)
log.rho  <- 0.025
log.prob <- c(0.05, 0.03, 0.01, 0.005, 0.005, 0.005)
ns       <- 1300

s <- cbind(runif(ns), runif(ns))
d <- as.matrix(rdist(s))
diag(d) <- 0

par(mfrow = c(3, 2))
for (i in 1:length(log.var)) {
  Sigma <- simple.cov.sp(D=d, sp.type="exponential", 
                         sp.par=c(log.var[i], log.rho), error.var=0,
                         finescale.var=0)
  data <- transform$logit(log.prob[i]) + t(chol(Sigma)) %*% rnorm(ns)
  y <- rbinom(n = ns, size = 1, prob = transform$inv.logit(data))
  
  plot(s[which(y != 1), ], pch = 21,
       col = "dodgerblue4", bg = "dodgerblue1", cex = 1.5,
       xlab = "", ylab = "",
       main = bquote(paste("Logit - 5%, ns = ", .(ns), " , ", 
                           sigma^2, "=", .(log.var[i]),  
                           ", Rareness = ", .(round(100 * mean(y), 2)), "%", 
                           sep = "")))
  
  points(s[which(y == 1), ], pch = 21, cex = 1.5,
         col = "firebrick4", bg = "firebrick1")
}