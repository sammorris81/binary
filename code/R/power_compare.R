rm(list=ls())
options(warn=2)
library(microbenchmark)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
Sys.setenv("PKG_CXXFLAGS"="-march=native -fopenmp")  # what if we change the march option
Sys.setenv("PKG_LIBS"="-march=native -fopenmp")
sourceCpp(file = "./power.cpp", rebuild = TRUE)

alpha <- 0.1

set.seed(2087)
temp.1 <- matrix(abs(rnorm(10000)), 100, 100)
temp.1 <- temp.1^(alpha)
mean(temp.1)

set.seed(2087)
temp.2 <- matrix(abs(rnorm(10000)), 100, 100)
temp.2 <- powCPP(temp.2, alpha)
mean(temp.2)

microbenchmark(temp.1^alpha, powCPP(temp.2, alpha))
