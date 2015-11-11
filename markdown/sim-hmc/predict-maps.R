# load packages and source files
rm(list=ls())
options(warn=2)
library(fields)
library(evd)
# library(spBayes)
library(fields)
library(SpatialTools)
# library(microbenchmark)  # comment out for beowulf
library(mvtnorm)
library(Rcpp)
library(numDeriv)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

source("./combine-tables.R")
source("../../code/R/spatial_gev.R", chdir = TRUE)
source("../../code/R/spatial_logit.R", chdir = TRUE)
source("../../code/R/spatial_probit.R", chdir = TRUE)

# get the datasets
load("./simdata.RData")

# extract the relevant setting from simdata
y <- simdata[[setting]]$y
s <- simdata[[setting]]$s
x <- simdata[[setting]]$x

# some precalculated values for quicker pairwise evaluation
dw2     <- rdist(s, knots)
d       <- rdist(s)
diag(d) <- 0

# sites at which to make predictions
s.p    <- expand.grid(seq(0, 1, length = 51), seq(0, 1, length = 51))
X.p    <- matrix(1, nrow(s.p), 1)

iters <- 20000; burn <- 10000; update <- 1000; thin <- 1

settings <- c(1:4)
gev.map <- array(0, dim=c(iters - burn, nrow(s.p), length(settings)))
pro.map <- array(0, dim=c(iters - burn, nrow(s.p), length(settings)))
log.map <- array(0, dim=c(iters - burn, nrow(s.p), length(settings)))
for (setting in settings) {
  set <- which(results[[setting]][, 1] - results[[setting]][, 2] == 
                 min(results[[setting]][, 1] - results[[setting]][, 2], na.rm = TRUE))
  filename <- paste("sim-results/", setting, "-", set, ".RData", sep = "")
  load(filename)
  
  cat("Start gev \n")
  gev.map[, , setting] <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                                     s.pred = s.p, knots = knots, start = 1, 
                                     end = iters - burn, update = update)
  
  cat("Start probit \n")
  pro.map[, , setting] <- pred.spprob(mcmcoutput = fit.probit, X.pred = X.p, 
                                      s.pred = s.p, knots = knots, start = 1, 
                                      end = iters - burn, update = update)
  
  cat("Start logit \n")
  log.map[, , setting] <- pred.splogit(mcmcoutput = fit.logit, 
                                       s.pred = s.p, knots = knots, start = 1, 
                                       end = iters - burn, update = update)
}

gev.summary <- array(0, dim = c(nrow(s.p), 5, length(settings)))
pro.summary <- array(0, dim = c(nrow(s.p), 5, length(settings)))
log.summary <- array(0, dim = c(nrow(s.p), 5, length(settings)))
for (setting in settings) {
  cat("Start setting", setting, "\n")
  gev.summary[, 1, setting] <- apply(gev.map[, , setting], 2, quantile, probs = 0.025)
  gev.summary[, 2, setting] <- apply(gev.map[, , setting], 2, median)
  gev.summary[, 3, setting] <- apply(gev.map[, , setting], 2, quantile, probs = 0.975)
  gev.summary[, 4, setting] <- apply(gev.map[, , setting], 2, mean)
  gev.summary[, 5, setting] <- apply(gev.map[, , setting], 2, sd)
  
  pro.summary[, 1, setting] <- apply(pro.map[, , setting], 2, quantile, probs = 0.025)
  pro.summary[, 2, setting] <- apply(pro.map[, , setting], 2, median)
  pro.summary[, 3, setting] <- apply(pro.map[, , setting], 2, quantile, probs = 0.975)
  pro.summary[, 4, setting] <- apply(pro.map[, , setting], 2, mean)
  pro.summary[, 5, setting] <- apply(pro.map[, , setting], 2, sd)
  
  log.summary[, 1, setting] <- apply(log.map[, , setting], 2, quantile, probs = 0.025)
  log.summary[, 2, setting] <- apply(log.map[, , setting], 2, median)
  log.summary[, 3, setting] <- apply(log.map[, , setting], 2, quantile, probs = 0.975)
  log.summary[, 4, setting] <- apply(log.map[, , setting], 2, mean)
  log.summary[, 5, setting] <- apply(log.map[, , setting], 2, sd)
}

save(s.p, gev.summary, pro.summary, log.summary, file = "maps.RData")

# generate plots of posterior summaries
for (setting in settings) {
  dev.new(width = 7, height = 7)
  par(mfrow = c(2, 2))
  quilt.plot(x = s.p[, 1], y = s.p[, 2], z = gev.summary[, 3, setting], nrow = 51, ncol = 51, main = "GEV")
  quilt.plot(x = s.p[, 1], y = s.p[, 2], z = pro.summary[, 3, setting], nrow = 51, ncol = 51, main = "Probit")
  quilt.plot(x = s.p[, 1], y = s.p[, 2], z = log.summary[, 3, setting], nrow = 51, ncol = 51, main = "Logit")
  filename <- paste("plots/post-med-", setting, ".pdf", sep = "")
  dev.print(device = pdf, file = filename)
  dev.off
}