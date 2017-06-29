rm(list=ls())
source("./package_load-2.R", chdir = TRUE)

# get the datasets
load("./simdata-grid.RData")

# data setting and sets to include - written by bash script
samp.type <- "srs"
n <- 250
gen.method <- 3

source("./runsim.R", chdir = TRUE)