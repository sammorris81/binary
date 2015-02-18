#### get z from covariates
getZ <- function(xi, x.beta) {
  z <- -(1 - xi * x.beta)^(-1 / xi)
  return(z)
}

a2A <- function(FAC, a, alpha) {
  # theta is nxnF
  # s is nFxnt
  # alpha in (0,1)
  W <- FAC^(1 / alpha)
  if (length(a) == 1) {xxx <- W * a}
  if (length(a) > 1) {xxx <- W %*% a}
  return(xxx)
}

# get the kernel weighting
make.fac <- function(dw2, gamma) {
  rho2 <- exp(gamma)^2
  fac <- exp(-0.5 * dw2 / rho2)
  return(fac)
}

# standardize the kernel weights
fac2FAC <- function(x, single=F) {
  if (single) {x <- x / sum(x)}
  if (!single) {x <- sweep(x, 1, rowSums(x), "/")}
  return(x)
}



logd <- function(theta, v) {
  sum(log(theta) - theta * v)
}

#### Used to set the standard deviation for the candidate distribution
#### for the A terms in the random effect. When log(A) is large means
#### sd is smaller, and log(A) small means sd is larger.
get.level <- function(logs, cuts) {
  sum(logs > cuts) + 1
}

logdet <- function(X) {
  determinant(X)$modulus
}

trunc <- function(x, eps=0.1) {
  x <- ifelse(x < eps, eps, x)
  x <- ifelse(x > 1-eps, 1-eps, x)
  return(x)
}