# load packages and source files
rm(list=ls())
options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)
# library(microbenchmark)  # comment out for beowulf
library(mvtnorm)
library(Rcpp)
library(numDeriv)
library(pROC)
library(compiler)
enableJIT(3)

source("../../R/spatial_gev.R", chdir = TRUE)
source("../../R/spatial_logit.R", chdir = TRUE)
source("../../R/spatial_probit.R", chdir = TRUE)
load("plant_inventory.RData")

if (Sys.info()["nodename"] == "cwl-mth-sam-001") {
  Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
  Sys.setenv("PKG_LIBS" = "-fopenmp")
  openblas.set.num.threads(1)
  do.upload <- TRUE
} else if (Sys.info()["sysname"] == "Darwin") {
  do.upload <- TRUE
} else {
  do.upload <- FALSE
}