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

source("../../../code/R/spatial_gev.R", chdir = TRUE)
source("../../../code/R/spatial_logit.R", chdir = TRUE)
source("../../../code/R/spatial_probit.R", chdir = TRUE)

set.seed(7483)  # site
ns <- c(1650, 2300, 1650, 2300, 1650, 2300)  # 650 train and 1300 train
nsettings <- length(ns)  # storing y in a list

################################################################################
### gev settings
################################################################################
gev.alpha  <- 0.3
gev.rho    <- 0.05  # 1.5 x knot spacing
gev.xi     <- 0
gev.prob   <- 0.05
gev.thresh <- -log(-log(1 - gev.prob))  # thresh = -Intercept
knots <- as.matrix(expand.grid(x = seq(0, 1, length = 21), 
                               y = seq(0, 1, length = 21)))


################################################################################
### logit settings
###   We use log.prob to give a set-specific intercept for the logit function. 
###   Originally, we set the intercept based on logit(p), but in doing so, the 
###   generation of the binomial RV resulted in high variability in the rareness 
###   from set to set. The average rareness was around 5%, but we ended up with 
###   some sets that had 1% rareness and some with 16%.
################################################################################
log.var    <- 10
log.rho    <- 0.05
log.prob   <- 0.05  # used to set the intercept for the xbeta
log.thresh <- transform$logit(log.prob)
log.error  <- 0  # let the bernoulli r.v. take care of this noise

################################################################################
### hotspot settings. generates around 5% 1s.
###   I tried using a fixed radius for all the sets and it yields pretty high 
###   variability in %1s. So, I changed to a set-specific radius for the size of
###   the hotspot. If there are a higher number of hotspots, then the radius 
###   decreases to account for the fact that there are more hotspot locations. 
###   This has the benefit of keeping the %1s around 5% for all datasets, but 
###   reducing the variability quite a bit. With a fixed radius, we get an
###   average of around 5% rareness, but have some sets with around 2% rareness
###   and other sets with around 18% rareness. 
################################################################################
nhotspots <- 2
p <- 0.95   # P(Y=1|hot spot)
q <- 0.001  # P(Y=1|background)
r <- 0.08  # Hot spot radius
# hot.prob <- 0.045

# nhotspots <- 3
# p <- 0.85   # P(Y=1|hot spot)
# q <- 0.005  # P(Y=1|background)
# r <- 0.05   # Hot spot radius

##############################################
### generate the data
##############################################
simdata <- vector(mode = "list", length = nsettings)

nsets <- 100
set.seed(3282)  # data

for (setting in 1:nsettings) {
  simdata[[setting]]$y <- matrix(data = NA, nrow = ns[setting], ncol = nsets)
  simdata[[setting]]$thresh <- rep(NA, nsets)
  simdata[[setting]]$s <- array(runif(2 * nsets * ns[setting]), 
                                dim = c(ns[setting], 2, nsets))
  simdata[[setting]]$x <- matrix(1, ns[setting], 1)
  simdata[[setting]]$hotspots <- vector(mode = "list", length = nsets)
  
  for (set in 1:nsets) {
   
    ### GEV generation
    if (setting == 1 | setting == 2) {
      data <- rRareBinarySpat(x = simdata[[setting]]$x, 
                              s = simdata[[setting]]$s[, , set], 
                              knots = knots, beta = 0, xi = gev.xi, 
                              alpha = gev.alpha, rho = gev.rho, 
                              prob.success = gev.prob, thresh = gev.thresh)
      
      simdata[[setting]]$y[, set]    <- data$y
      simdata[[setting]]$thresh[set] <- data$thresh
    } 
    
    ### logit generation
    if (setting == 3 | setting == 4) {
      d <- as.matrix(rdist(simdata[[setting]]$s[, , set]))
      diag(d) <- 0
      Sigma <- simple.cov.sp(D = d, sp.type = "exponential", 
                             sp.par = c(log.var, log.rho), 
                             error.var = log.error, finescale.var = 0)
      
      data <- transform$logit(log.prob) + t(chol(Sigma)) %*% rnorm(ns[setting])
      # thresh <- quantile(data, probs = (1 - log.prob))
      data <- log.thresh + data
      data <- rbinom(n = ns[setting], size = 1, prob = transform$inv.logit(data))
      
      simdata[[setting]]$y[, set]    <- data
      simdata[[setting]]$thresh[set] <- log.thresh
    }
    
    ### hotspot generation
    if (setting == 5 | setting == 6) {
      k  <- rpois(1, nhotspots) + 1
      hotspots <- cbind(runif(k), runif(k))
      d <- rdist(simdata[[setting]]$s[, , set], hotspots)
      
      # get the radius for the hotspots.
      #   1. Look at the distance to the closes knot for all sites
      #   2. Set the hotspot radius to the quantile of the minimum distances 
      #      that corresponds to the desired rareness / P(Y = 1|in hotspot)
      # r <- quantile(apply(d, 1, min), probs = hot.prob / p)
      
      hot <- rowSums(d <= r) > 0

      simdata[[setting]]$y[, set] <- rbinom(ns[setting], 1, ifelse(hot, p, q))
      simdata[[setting]]$hotspots[[set]] <- hotspots
      simdata[[setting]]$thresh[set] <- r
    }
    
    if (set %% 20 == 0) {
      print(paste("  Dataset ", set, " finished", sep = ""))
    }
  }
  print(paste("Setting ", setting, " finished", sep = ""))
}

################################################################################
# Split observations for testing and training. Using a stratified sample here to 
# further reduce the variability in rareness across training and testing sites.
################################################################################
for (setting in 1:nsettings) {
  ntest  <- 1000
  ntrain <- ns[setting] - ntest
  for (set in 1:nsets) {
    # simplify the notation
    y <- simdata[[setting]]$y[, set]
    s <- simdata[[setting]]$s[, , set]
    
    # get the sample breakdown of 1s and 0s for testing
    phat <- mean(y)
    ntest.1 <- floor(ntest * phat)
    ntest.0 <- ntest - ntest.1
    idx.1   <- sample(which(simdata[[setting]]$y[, set] == 1), size = ntest.1)
    idx.0   <- sample(which(simdata[[setting]]$y[, set] == 0), size = ntest.0)
    
    # reorder y and s so the train are the first observations
    test  <- c(idx.1, idx.0)
    train <- (1:ns[setting])[-test]
    
    simdata[[setting]]$y[, set]   <- y[c(train, test)]
    simdata[[setting]]$s[, , set] <- s[c(train, test), ]
  }
}

for (i in 1:6) {
  print(mean(simdata[[i]]$y))
}

#### Plot datasets for different settings with highest and lowest rareness
dev.new(width = 12, height = 9)
par(mfrow = c(3, 4))
settings <- c("GEV", "GEV", "Logit", "Logit", "Hotspot", "Hotspot")

for (setting in 1:length(settings)) {
  if (setting %in% c(1, 3, 5)) {
    end <- 650
  } else {
    end <- 1300
  }
  
  sets <- tail(order(colMeans(simdata[[setting]]$y[1:end, ])), 2)
  for(set in sets) {
    plot(simdata[[setting]]$s[which(simdata[[setting]]$y[1:end, set] != 1), , set], 
         pch = 21, cex = 1.5, col = "dodgerblue4", bg = "dodgerblue1", 
         xlab = "", ylab = "", 
         main = paste(settings[setting], ": ", 
                      round(100 * mean(simdata[[setting]]$y[1:end, set]), 2), 
                      "%, ns = ", ns[setting] - 1000, sep = ""))
    points(simdata[[setting]]$s[which(simdata[[setting]]$y[1:end, set] == 1), , set], 
           pch = 21, cex = 1.5, col = "firebrick4", bg = "firebrick1")
  }
}
dev.print(device = pdf, file = "five-high.pdf")
dev.off()

dev.new(width = 12, height = 9)
par(mfrow = c(3, 4))
settings <- c("GEV", "GEV", "Logit", "Logit", "Hotspot", "Hotspot")

for (setting in 1:length(settings)) {
  if (setting %in% c(1, 3, 5)) {
    end <- 650
  } else {
    end <- 1300
  }
  
  sets <- order(colMeans(simdata[[setting]]$y[1:end, ]))[1:2]
  
  for (set in sets) {
    plot(simdata[[setting]]$s[which(simdata[[setting]]$y[1:end, set] != 1), , set], 
         pch = 21, cex = 1.5, col = "dodgerblue4", bg = "dodgerblue1", 
         xlab = "", ylab = "", 
         main = paste(settings[setting], ": ", 
                      round(100 * mean(simdata[[setting]]$y[1:end, set]), 2), 
                      "%, ns = ", ns[setting] - 1000, sep = ""))
    points(simdata[[setting]]$s[which(simdata[[setting]]$y[1:end, set] == 1), , set], 
           pch = 21, cex = 1.5, col = "firebrick4", bg = "firebrick1")
  }
}
dev.print(device = pdf, file = "five-low.pdf")
dev.off()

save(simdata, gev.rho, gev.prob, log.rho, log.prob, file = "simdata.RData")

# for processing over many machines at once
sets.remain <- matrix(TRUE, nsets, nsettings)
write.table(x = sets.remain, file = "./sim-control/sets-remain.txt")
system(paste("scp ./sim-control/sets-remain.txt samorris@hpc.stat.ncsu.edu:~/",
       "repos-git/rare-binary/markdown/dec2015/third-run-sim/sim-control/", 
       sep = ""))

# ns <- 1300
# p <- 0.01
# s <- cbind(runif(ns), runif(ns))
# d <- as.matrix(rdist(s))
# diag(d) <- 0
# Sigma <- simple.cov.sp(D=d, sp.type="exponential", 
#                        sp.par=c(7, 0.075), error.var=1, 
#                        finescale.var=0)
# # data <- log(p / (1 - p)) + t(chol(Sigma)) %*% rnorm(ns)
# data <- t(chol(Sigma)) %*% rnorm(ns)
# data <- data - quantile(data, probs = 0.97)  # uses a little padding
# data <- rbinom(n = ns, size = 1, prob = 1 / (1 + exp(-data)))
# 
# plot(s[which(data != 1), ], pch = 21, cex = 1.5, 
#      col = "dodgerblue4", bg = "dodgerblue1", xlab = "", ylab = "", 
#      main = paste("P(Y = 1) = ", round(mean(data), 4), sep = ""))
# points(s[which(data == 1), ], pch = 21, cex = 1.5, 
#        col = "firebrick4", bg = "firebrick1")
# 
# z <- log(p / (1 - p)) + t(chol(Sigma)) %*% rnorm(ns)
# hist(z)
# mean(z > 0)

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


# library(fields)
# library(SpatialTools)
# log.var  <- c(1 ,3, 5, 7, 9, 11)
# log.rho  <- 0.025
# log.prob <- c(0.05, 0.03, 0.01, 0.005, 0.005, 0.005)
# ns       <- 1300
# 
# s <- cbind(runif(ns), runif(ns))
# d <- as.matrix(rdist(s))
# diag(d) <- 0
# 
# par(mfrow = c(3, 2))
# for (i in 1:length(log.var)) {
#   Sigma <- simple.cov.sp(D=d, sp.type="exponential", 
#                          sp.par=c(log.var[i], log.rho), error.var=0,
#                          finescale.var=0)
#   data <- transform$logit(log.prob[i]) + t(chol(Sigma)) %*% rnorm(ns)
#   y <- rbinom(n = ns, size = 1, prob = transform$inv.logit(data))
#   
#   plot(s[which(y != 1), ], pch = 21,
#        col = "dodgerblue4", bg = "dodgerblue1", cex = 1.5,
#        xlab = "", ylab = "",
#        main = bquote(paste("Logit - 5%, ns = ", .(ns), " , ", 
#                            sigma^2, "=", .(log.var[i]),  
#                            ", Rareness = ", .(round(100 * mean(y), 2)), "%", 
#                            sep = "")))
#   
#   points(s[which(y == 1), ], pch = 21, cex = 1.5,
#          col = "firebrick4", bg = "firebrick1")
# }