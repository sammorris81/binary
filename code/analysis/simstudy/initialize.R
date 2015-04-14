options(warn=2)
library(fields)
library(evd)
library(spBayes)
library(fields)
library(SpatialTools)
library(MCMCpack)

rm(dPS.Rcpp)  # included to force recompile on local machine
source("../../R/auxfunctions.R", chdir=T)
source("../../R/updateModel.R", chdir=T)
source("../../R/mcmc.R", chdir=T)
source("../../R/probit.R", chdir=T)