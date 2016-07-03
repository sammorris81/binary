# get the datasets
load(paste("./", species, ".RData", sep = ""))

# get the correct y, x, and s
if (species == "cattle_egret") {
  y.o <- cattle_egret[cv.idx[[cv]]]
  y.p <- cattle_egret[-cv.idx[[cv]]]
  seed.base <- 1000
} else if (species == "common_nighthawk") {
  y.o <- common_nighthawk[cv.idx[[cv]]]
  y.p <- common_nighthawk[-cv.idx[[cv]]]
  seed.base <- 2000
} else if (species == "vesper_sparrow") {
  y.o <- vesper_sparrow[cv.idx[[cv]]]
  y.p <- vesper_sparrow[-cv.idx[[cv]]]
  seed.base <- 3000
} else if (species == "western_bluebird") {
  y.o <- western_bluebird[cv.idx[[cv]]]
  y.p <- western_bluebird[cv.idx[[cv]]]
  seed.base <- 4000
} else {
  stop("incorrect species selected")
}

if (knot.percent == 10) {
  knots <- knots.10[[cv]]
  seed.knot <- 100
} else if (knot.percent == 15) {
  knots <- knots.15[[cv]]
  seed.knot <- 200
} else if (knot.percent == 20) {
  knots <- knots.20[[cv]]
  seed.knot <- 300
} else {
  stop("only can use knots.percent = 10, 15, or 20")
}

# extract info about simulation settings
ns     <- length(y.o)
npred  <- length(y.p)
nt     <- 1
nknots <- nrow(knots)

# scale sites so in [0, 1] x [0, 1] (or close)
s.min <- apply(s, 2, min)
s.max <- apply(s, 2, max)
s.range <- c(diff(range(s[, 1])), diff(range(s[, 2])))
s.scale <- s
s.scale[, 1] <- (s[, 1] - s.min[1]) / max(s.range)
s.scale[, 2] <- (s[, 2] - s.min[2]) / max(s.range)
knots[, 1] <- (knots[, 1] - s.min[1]) / max(s.range)
knots[, 2] <- (knots[, 2] - s.min[2]) / max(s.range)

y.o <- matrix(y.o, ns, nt)
s.o <- s.scale[cv.idx[[cv]], ]
X.o <- matrix(1, nrow(s.o), 1)
y.p <- matrix(y.p, npred, nt)
s.p <- s.scale[-cv.idx[[cv]], ]
X.p <- matrix(1, nrow(s.p), 1)

####################################################################
#### Start MCMC setup: Most of this is used for the spBayes package
####################################################################
iters <- 25000; burn <- 15000; update <- 1000; thin <- 1; iterplot <- FALSE
# iters <- 10000; burn <- 5000; update <- 100; thin <- 1; iterplot <- TRUE
n.report     <- 10
batch.length <- 100
n.batch      <- floor(iters / batch.length)
verbose      <- TRUE
tuning       <- list("phi" = 0.1, "sigma.sq" = 0.2, "beta" = 1, "w" = 5)
starting     <- list("phi" = 3/0.5, "sigma.sq" = 50, "beta" = 0, "w" = 0)
priors       <- list("beta.norm" = list(0, 100),
                     "phi.unif" = c(0.1, 1e4), "sigma.sq.ig" = c(1, 1))
cov.model <- "exponential"
timings   <- rep(NA, 3)
# with so many knots, adaptive is time prohibitive
amcmc     <- list("n.batch" = n.batch, "batch.length" = batch.length,
                  "accept.rate" = 0.35)

# storage for some of the results
scores <- matrix(NA, 3, 2)  # place to store brier scores and auc
rownames(scores) <- c("gev", "probit", "logit")
colnames(scores) <- c("bs", "auc")

timings <- rep(NA, 3)

# start the simulation
set.seed(seed.base + seed.knot + cv * 10)

rho.init.pcl <- 0.05  
dw2.o     <- rdist(s.o, knots)^2
d.o       <- as.matrix(rdist(s.o))
diag(d.o) <- 0
max.dist  <- 0.20

#### spatial GEV
cat("  Start gev \n")

# using pairwise estimates as starting points for rho, alpha, and beta. also 
# using the the pairwise estimate of alpha as the mean of the prior 
# distribution along with a standard deviation of 0.05 to allow for some 
# variability, but hopefully also some better convergence w.r.t. alpha.
# we set max.dist to 0.15 in order to only consider pairs of sites that 
# are relatively close to one another.
cat("    Start pairwise fit \n")
fit.pcl <- tryCatch(
  fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                    alpha.init = 0.5, rho.init = rho.init.pcl,
                    xi.fix = TRUE, alpha.fix = FALSE,
                    rho.fix = FALSE, beta.fix = TRUE,
                    y = y.o, dw2 = dw2.o, d = d.o,
                    cov = X.o, method = "BFGS",
                    max.dist = max.dist,  
                    alpha.min = 0.1, alpha.max = 0.9,
                    threads = 2),
  error = function(e) {
    fit.rarebinaryCPP(beta.init = 0, xi.init = 0,
                      alpha.init = 0.5, rho.init = rho.init.pcl,
                      xi.fix = TRUE, alpha.fix = FALSE,
                      rho.fix = FALSE, beta.fix = TRUE,
                      y = y.o, dw2 = dw2.o, d = d.o,
                      cov = X.o, method = "Nelder-Mead",
                      max.dist = max.dist,
                      alpha.min = 0.1, alpha.max = 0.9,
                      threads = 2)
  }
)

cat("    Finish pairwise fit \n")

cat("    Start mcmc fit \n")
mcmc.seed <- 6262 + seed.base + seed.knot + cv * 10
set.seed(mcmc.seed)

alpha.mn <- fit.pcl$par[1]

# for numerical stability with the current set of starting values for the a 
# terms. if alpha is too small, the algorithm has a very hard time getting 
# started.
if (alpha.mn < 0.3) {
  alpha.init <- 0.3  
} else {
  alpha.init <- alpha.mn
}

rho.init <- fit.pcl$par[2]
beta.init <- fit.pcl$beta

fit.gev <- spatial_GEV(y = y.o, s = s.o, x = X.o, knots = knots, 
                       beta.init = log(-log(1 - mean(y.o))),
                       beta.mn = 0, beta.sd = 10,
                       beta.eps = 0.1, beta.attempts = 50, 
                       xi.init = 0, xi.mn = 0, xi.sd = 0.5, xi.eps = 0.01, 
                       xi.attempts = 50, xi.fix = FALSE, 
                       a.init = 1, a.eps = 0.2, a.attempts = 50, 
                       a.cutoff = 0.2, b.init = 0.5, b.eps = 0.2, 
                       b.attempts = 50, 
                       alpha.init = alpha.init, alpha.attempts = 50, 
                       alpha.mn = alpha.mn, alpha.sd = 0.05,
                       a.alpha.joint = TRUE, alpha.eps = 0.001,
                       rho.init = rho.init, logrho.mn = -2, logrho.sd = 1, 
                       rho.eps = 0.1, rho.attempts = 50, threads = 1, 
                       iters = iters, burn = burn, 
                       update = update, iterplot = iterplot,
                       # update = 10, iterplot = TRUE,
                       thin = thin, thresh = 0)

cat("    Start mcmc predict \n")
post.prob.gev <- pred.spgev(mcmcoutput = fit.gev, x.pred = X.p,
                            s.pred = s.p, knots = knots,
                            start = 1, end = iters - burn, update = update)
timings[1] <- fit.gev$minutes

bs.gev <- BrierScore(post.prob.gev, y.p)
post.prob.gev.med <- apply(post.prob.gev, 2, median)
roc.gev <- roc(y.p ~ post.prob.gev.med)
auc.gev <- roc.gev$auc

print(bs.gev * 100)


upload.pre <- paste("samorris@hpc.stat.ncsu.edu:~/repos-git/rare-binary/code/",
                    "analysis/", sep = "")
if (do.upload) {
  upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  system(upload.cmd)
}
save(blah, results.file)