# load packages and source files
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

source("../../../code/R/spatial_gev.R", chdir = TRUE)
source("../../../code/R/spatial_logit.R", chdir = TRUE)
source("../../../code/R/spatial_probit.R", chdir = TRUE)

if (Sys.info()["nodename"] == "cwl-mth-sam-001") {
  Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
  Sys.setenv("PKG_LIBS"="-fopenmp")
  openblas.set.num.threads(1)
  do.upload <- TRUE
} else if (Sys.info()["sysname"] == "Darwin") {
  do.upload <- TRUE
} else {
  do.upload <- FALSE
}

upload.pre <- paste("samorris@hpc.stat.ncsu.edu:~/repos-git/rare-binary/code/",
                    "analysis/simstudy/sim-tables/", sep = "")

# the directory in the sim study that controls what sets are left
if (do.upload) {
  control.dir <- paste("samorris@hpc.stat.ncsu.edu:~/repos-git/rare-binary/",
                       "code/analysis/simstudy/sim-control/", sep = "")
} else {
  control.dir <- "./sim-control/"
}

# we need this for locking while other scripts are selecting their set
lock.command <- paste("touch ~/repos-git/rare-binary/code/analysis/simstudy/",
                      "sim-control/lock.txt", sep = "")
if (do.upload) {
  lock.command.ssh <- paste("ssh samorris@hpc.stat.ncsu.edu '", lock.command,
                            "'", sep = "")
}

unlock.command <- paste("rm ~/repos-git/rare-binary/code/analysis/simstudy/",
                        "sim-control/lock.txt", sep = "")
if (do.upload) {
  unlock.command.ssh <- paste("ssh samorris@hpc.stat.ncsu.edu '",
                              unlock.command, "'", sep = "")
}