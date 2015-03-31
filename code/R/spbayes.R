library(spBayes)
library(fields)
library(SpatialTools)

## Not run: 
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(1)

n <- 100
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- as.matrix(cbind(1, rnorm(n)))

B <- as.matrix(c(1,5))
p <- length(B)

sigma.sq <- 2
tau.sq <- 0.1
phi <- 3/0.5

###########################
##Predictive process model
###########################
# note phi is 1 / rho from our typical representation
# right now using phi = 3 / 0.5 or rho = 1 / 6
# phi.unif of U(2/3, 1e6) gives rho.unif of U(0.000001, 1.5)
# requires that nu.unif is positive

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
R <- simple.cov.sp(D=D, sp.type="matern", sp.par=c(1, 1/phi), error.var=0, smoothness=0.4, finescale.var=0)
w <- rmvn(1, rep(0,n), sigma.sq*R)  # spatial component
y <- rnorm(n, X%*%B + w, sqrt(tau.sq))

n.samples <- 2000

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
                 "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))
cov.model <- "exponential"

n.report <- 500
verbose <- TRUE

# fits the model at all locations, no dimension reduction
m.1 <- spLM(y~X-1, coords=coords, starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)


starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1, "nu"=0.5)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1, "nu"=0.1)
priors <- list("beta.norm"=list(rep(0,p), diag(1000,p)),
               "phi.unif"=c(0.5, 1e6), "sigma.sq.ig"=c(1, 1),
               "tau.sq.ig"=c(1, 1), "nu.unif"=c(1e-6, 10))
cov.model <- "matern"

# implements dimension reduction
m <- spLM(y~X-1, coords=coords, knots=c(6,6,0.1), starting=starting,
          tuning=tuning, priors=priors, cov.model=cov.model,
          n.samples=n.samples, verbose=verbose, n.report=n.report)


# move to binary to check timing
hist(y)
y.bin <- ifelse(y > 6, 1, 0)
mean(y.bin) # 0.16

m <- spGLM(y.bin~X-1, family="binomial", coords=coords, knots=c(6, 6, 0.1), 
           starting=starting, tuning=tuning, priors=priors, cov.model=cov.model,
           n.samples=n.samples, verbose=verbose, n.report=n.report)


burn.in <- 0.5*n.samples

round(summary(window(m.1$p.beta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
round(summary(window(m.2$p.beta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)

round(summary(window(m.1$p.theta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
round(summary(window(m.2$p.theta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)

## End(Not run)
