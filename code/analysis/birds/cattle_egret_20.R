rm(list=ls())
source(file = "./package_load.R", chdir = T)
species <- "cattle_egret"
cv <- 1
knot.percent <- 20
source(file = "./fitmodel.R")

rm(list = ls())
source(file = "./package_load.R", chdir = T)
species <- "cattle_egret"
cv <- 2
knot.percent <- 20
source(file = "./fitmodel.R")