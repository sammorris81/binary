rm(list = ls())
source("./package_load.R", chdir = TRUE)

set.seed(7483)  # site
# ns <- c(1650, 2300, 1650, 2300, 1650, 2300)  # 650 train and 1300 train
# ns <- c(1100, 1250, 1100, 1250, 1100, 1250)  # 100 train and 250 train
ns <- c(100, 250)
samp.types  <- c("clu", "srs")
gen.methods <- c("gev", "logistic", "hotspot")
nmethods <- length(gen.methods)  # storing y in a list

################################################################################
### gev settings
################################################################################
# gev.alpha  <- 0.3
gev.alpha  <- 0.4
# gev.rho    <- 0.05  # 1.5 x knot spacing
gev.rho    <- 0.025
gev.xi     <- 0
gev.prob   <- 0.05
gev.thresh <- -log(-log(1 - gev.prob))  # thresh = -Intercept
knots <- as.matrix(expand.grid(x = seq(1 / 60, 59 / 60, length = 30), 
                               y = seq(1 / 60, 59 / 60, length = 30)))

################################################################################
### logit settings
###   We use log.prob to give a set-specific intercept for the logit function.
###   Originally, we set the intercept based on logit(p), but in doing so, the
###   generation of the binomial RV resulted in high variability in the rareness
###   from set to set. The average rareness was around 5%, but we ended up with
###   some sets that had 1% rareness and some with 16%.
################################################################################
log.var    <- 10
# log.rho    <- 0.05
log.rho    <- 0.025
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
p <- 0.90   # P(Y=1|hot spot)
q <- 0.005  # P(Y=1|background)
r <- 0.07  # Hot spot radius
# hot.prob <- 0.045

# nhotspots <- 3
# p <- 0.85   # P(Y=1|hot spot)
# q <- 0.005  # P(Y=1|background)
# r <- 0.05   # Hot spot radius

##############################################
### generate the data
##############################################
simdata <- vector(mode = "list", length = nmethods)

nsets <- 100
set.seed(3282)  # data
s.grid <- as.matrix(expand.grid(seq(0, 1, length = 100), 
                                seq(0, 1, length = 100)))
d <- rdist(s.grid)
diag(d) <- 0
Sigma <- simple.cov.sp(D = d, sp.type = "exponential",
                       sp.par = c(log.var, log.rho),
                       error.var = log.error, finescale.var = 0)
t.Sigma.chol <- t(chol(Sigma))

for (method in 1:nmethods) {
  simdata[[method]]$y.grid <- matrix(data = NA, nrow = nrow(s.grid), ncol = nsets)
  simdata[[method]]$thresh <- rep(NA, nsets)
  simdata[[method]]$x <- matrix(1, nrow(s.grid), 1)
  simdata[[method]]$hotspots <- vector(mode = "list", length = nsets)
  
  for (set in 1:nsets) {
    ### GEV generation
    if (method == 1) {
      nobs <- 0
      while (nobs < 100 | nobs > 700) {
        data <- rRareBinarySpat(x = simdata[[method]]$x, s = s.grid,
                                knots = knots, beta = 0, xi = gev.xi,
                                alpha = gev.alpha, rho = gev.rho,
                                prob.success = gev.prob, thresh = gev.thresh)
        
        simdata[[method]]$y.grid[, set] <- data$y
        simdata[[method]]$thresh[set]   <- data$thresh
        
        nobs <- sum(data$y)
      }
    }
    
    ### logit generation
    if (method == 2) {
      nobs <- 0
      while (nobs < 100 | nobs > 700) {
        data <- transform$logit(log.prob) + t.Sigma.chol %*% rnorm(nrow(s.grid))
        data <- log.thresh + data
        data <- rbinom(n = nrow(s.grid), size = 1, prob = transform$inv.logit(data))
        
        simdata[[method]]$y.grid[, set] <- data
        simdata[[method]]$thresh[set]   <- log.thresh
        
        nobs <- sum(data)
      }
    }
    
    ### hotspot generation
    if (method == 3) {
      nobs <- 0
      while (nobs < 100 | nobs > 700) {
        k  <- rpois(1, nhotspots) + 1
        hotspots <- cbind(runif(k), runif(k))
        d <- rdist(s.grid, hotspots)
        
        # get the radius for the hotspots.
        #   1. Look at the distance to the closes knot for all sites
        #   2. Set the hotspot radius to the quantile of the minimum distances
        #      that corresponds to the desired rareness / P(Y = 1|in hotspot)
        # r <- quantile(apply(d, 1, min), probs = hot.prob / p)
        
        hot <- rowSums(d <= r) > 0
        this.y <- rbinom(nrow(s.grid), 1, ifelse(hot, p, q))
        
        simdata[[method]]$hotspots[[set]] <- hotspots
        simdata[[method]]$thresh[set] <- r
        
        nobs <- sum(this.y)
        simdata[[method]]$y.grid[, set] <- this.y
      }
    }
    
    if (set %% 20 == 0) {
      print(paste("  Dataset ", set, " finished", sep = ""))
    }
  }
  print(paste("Method ", gen.methods[method], " finished", sep = ""))
}
save(s.grid, simdata, file = "simdata-grid.RData")

################################################################################
# Split observations for testing and training. Using a stratified sample here to
# further reduce the variability in rareness across training and testing sites.
################################################################################
for (method in 1:3) {
  for (n.idx in 1:length(ns)) {
    n <- ns[n.idx]
    this.clu <- vector(mode = "list", length = nsets)
    this.srs <- vector(mode = "list", length = nsets)
    clu.name <- paste("clu.lst.", method, ".", n, sep = "")
    srs.name <- paste("srs.lst.", method, ".", n, sep = "")
    for (set in 1:nsets) {
      set.seed(n.idx * 1000 + set)
      nobs <- 0
      this.y <- simdata[[method]]$y.grid[, set]
      
      while (nobs < 3) {
        # keep repeating the sampling until there are at least 3 observations
        these.train <- sort(sample(length(this.y), n))
        y.o <- this.y[these.train]
        these.cluster <- y.o == 1
        these.cluster.ids <- these.train[these.cluster]
        for (i in 1:length(these.cluster.ids)) {
          # get the location of the cell where y == 1
          this.cell <- get.arr.idx(these.cluster.ids[i], nrows = 100)
          this.row  <- this.cell[1]
          this.col  <- this.cell[2]
          
          # account for the boundary
          neighbors.row <- c(this.row + 1, this.row - 1)
          neighbors.col <- c(this.col + 1, this.col - 1)
          neighbors.row <- neighbors.row[neighbors.row > 0 & neighbors.row < 100]
          neighbors.col <- neighbors.col[neighbors.col > 0 & neighbors.col < 100]
          for (j in 1:length(neighbors.row)) {
            these.train <- c(these.train,
                             get.idx(row = neighbors.row[j], col = this.col, 
                                     nrows = 100))
          }
          for (j in 1:length(neighbors.col)) {
            these.train <- c(these.train, 
                             get.idx(row = this.row, col = neighbors.col[j], 
                                     nrows = 100))
          }
        }
        these.train <- sort(unique(these.train))
        y.o <- this.y[these.train]
        nobs <- sum(y.o)
      }
      this.clu[[set]] <- these.train
      
      set.seed(n.idx * 1000 + set)
      nobs <- 0
      nsamp <- length(these.train)
      while (nobs < 3) {
        these.train <- sort(sample(length(this.y), nsamp))
        y.o <- this.y[these.train]
        nobs <- sum(y.o)
      }
      this.srs[[set]] <- these.train
    }
    
    assign(clu.name, this.clu)
    assign(srs.name, this.srs)
  }
}

save(simdata, s.grid, 
     srs.lst.1.100, srs.lst.2.100, srs.lst.3.100, 
     srs.lst.1.250, srs.lst.2.250, srs.lst.3.250,
     clu.lst.1.100, clu.lst.2.100, clu.lst.3.100,
     clu.lst.1.250, clu.lst.2.250, clu.lst.3.250,
     file = "simdata-grid.RData")

for (i in 1:6) {
  if (i %in% c(1, 3, 5)) {
    print(mean(simdata[[i]]$y[1:100, ]))
  } else {
    print(mean(simdata[[i]]$y[1:250, ]))
  }
}

#### Plot datasets for different settings with highest and lowest rareness
dev.new(width = 12, height = 9)
par(mfrow = c(3, 4))
settings <- c("GEV", "GEV", "Logit", "Logit", "Hotspot", "Hotspot")

for (setting in 1:length(settings)) {
  end <- ns[setting]
  
  sets <- tail(order(colMeans(simdata[[setting]]$y[1:end, ])), 2)
  for(set in sets) {
    plot(simdata[[setting]]$s[which(simdata[[setting]]$y[1:end, set] != 1), , set],
         pch = 21, cex = 1, col = "dodgerblue4", bg = "dodgerblue1",
         xlab = "", ylab = "",
         main = paste(settings[setting], ": ",
                      round(100 * mean(simdata[[setting]]$y[1:end, set]), 2),
                      "%, ns = ", ns[setting] - 1000, sep = ""))
    points(simdata[[setting]]$s[which(simdata[[setting]]$y[1:end, set] == 1), , set],
           pch = 21, cex = 1, col = "firebrick4", bg = "firebrick1")
  }
}
dev.print(device = pdf, file = "five-high.pdf")
dev.off()

dev.new(width = 12, height = 9)
par(mfrow = c(3, 4))
settings <- c("GEV", "GEV", "Logit", "Logit", "Hotspot", "Hotspot")

for (setting in 1:length(settings)) {
  end <- ns[setting]
  
  sets <- order(colMeans(simdata[[setting]]$y[1:end, ]))[1:2]
  
  for (set in sets) {
    plot(simdata[[setting]]$s[which(simdata[[setting]]$y[1:end, set] != 1), , set],
         pch = 21, cex = 1, col = "dodgerblue4", bg = "dodgerblue1",
         xlab = "", ylab = "",
         main = paste(settings[setting], ": ",
                      round(100 * mean(simdata[[setting]]$y[1:end, set]), 2),
                      "%, ns = ", ns[setting] - 1000, sep = ""))
    points(simdata[[setting]]$s[which(simdata[[setting]]$y[1:end, set] == 1), , set],
           pch = 21, cex = 1, col = "firebrick4", bg = "firebrick1")
  }
}
dev.print(device = pdf, file = "five-low.pdf")
dev.off()

save(simdata, gev.rho, gev.prob, log.rho, log.prob, file = "simdata.RData")

# for processing over many machines at once
nsets <- 50
nsettings <- 12
sets.remain <- matrix(TRUE, nsets, nsettings, byrow = TRUE)
write.table(x = sets.remain, file = "./sim-control-2/sets-remain.txt")
system(paste("scp ./sim-control-2/sets-remain.txt samorris@hpc.stat.ncsu.edu:~/",
             "repos-git/rare-binary/code/analysis/simstudy/sim-control-2/",
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