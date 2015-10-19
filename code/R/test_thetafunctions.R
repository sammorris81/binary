nknots <- 5
nt <- 2
ns <- 4
alpha <- 0.4
temp.a <- abs(matrix(rnorm(nknots * nt), nknots, nt))
temp.w <- abs(matrix(rnorm(ns * nknots), ns, nknots))
  
aw <- matrix(0, nrow = ns, ncol = nt)
for (t in 1:nt) {
  for (i in 1:ns) {
    for (l in 1:nknots) {
      aw[i, t] <- aw[i, t] + temp.a[l, t] * temp.w[i, l]^(1 / alpha)
    }
  }
}

getawCPP(w = temp.w, a_star = temp.a^alpha, alpha = alpha) / aw
  