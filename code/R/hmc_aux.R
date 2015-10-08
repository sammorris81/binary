source("./HMC.R")

neg_log_post <- function(a, b, others) {
  # others is a list
  #   y: data
  #   alpha: spatial dependence
  #   wz: kernel weights
  
  lc <- logc(b = b, alpha = others$alpha)
  
  ll <- -log(others$alpha) - log(a) / (1 - others$alpha) - log(1 - others$alpha) + 
  
  theta <- getThetaCPP(wz = others$wz, a_star = a^others$alpha, 
                       alpha = others$alpha)
  
  
  ll <- ll + sum(dbinom(o$y, 1, prob = exp(-theta), log = TRUE))
  
  # numerical stability issue. originally was using
  return (-ll)
}

logc <- function(b, alpha) {
  results <- alpha * log(sin(alpha * pi * b)) / (1 - alpha) - 
             log(sin(pi * b)) / (1 - alpha) + 
             log(sin(1 - alpha) * pi * b)
  
  return (results)
}

ns <- 2000
nknots <- 900
temp.w <- matrix(rnorm(ns * nknots), ns, nknots)
temp.z <- rnorm(ns)
temp.a <- rnorm(nknots)

microbenchmark(sweep(temp.w, 1, temp.z, "/"), temp.w / temp.z)

alpha <- 0.5
temp.a.star <- temp.a^alpha
microbenchmark(sweep(temp.w, 2, temp.a^alpha, "*"), 
               t(t(temp.w) * temp.a^alpha),
               sweep(temp.w, 2, temp.a.star, "*"),
               t(t(temp.w) * temp.a.s))
