neg_log_post <- function(a, b, others) {
  # o is a list
  #   y: data
  #   alpha: spatial dependence
  #   wz: kernel weights
  
  tem
  
  ll <- sum(dbinom(o$y, 1, prob = prob, log = TRUE))
  
  # numerical stability issue. originally was using
  
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
