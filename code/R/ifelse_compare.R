rm(list=ls())
options(warn=2)
library(microbenchmark)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
Sys.setenv("PKG_CXXFLAGS"="-march=native -fopenmp")  # what if we change the march option
Sys.setenv("PKG_LIBS"="-march=native -fopenmp")
sourceCpp(file = "./ifelse.cpp", rebuild = TRUE)

set.seed(2087)
temp <- matrix(rnorm(10000), 100, 100)
temp <- ifelse(abs(temp) < 2, 0, temp)
mean(temp)

set.seed(2087)
temp <- matrix(rnorm(10000), 100, 100)
temp <- ifelsematCPP(temp, 2)
mean(temp)

microbenchmark(ifelse(abs(temp) < 2, 0, temp), ifelsematCPP(temp, 2))

set.seed(2087)
temp <- rnorm(10000)
temp <- ifelse(abs(temp) < 2, 0, temp)
mean(temp)

set.seed(2087)
temp <- rnorm(10000)
temp <- ifelsevecCPP(temp, 2)
mean(temp)

microbenchmark(ifelse(abs(temp) < 2, 0, temp), ifelsematCPP(temp, 2))

