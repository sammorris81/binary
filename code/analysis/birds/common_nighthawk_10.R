rm(list=ls())
source(file = "./package_load.R", chdir = T)
species <- "common_nighthawk"
cv <- 1
knot.percent <- 10
results.file <- paste("./cv-results/", species, "-", knot.percent, "-", cv,
                      ".RData", sep = "")
table.file   <- paste("./cv-tables/", species, "-", knot.percent, "-", cv,
                      ".RData", sep = "")
source(file = "fitmodel.R")

rm(list = ls())
source(file = "./package_load.R", chdir = T)
cv <- 2
knot.percent <- 10
species <- "cattle_egret"
results.file <- paste("./cv-results/", species, "-", knot.percent, "-", cv,
                      ".RData", sep = "")
table.file   <- paste("./cv-tables/", species, "-", knot.percent, "-", cv,
                      ".RData", sep = "")
source(file = "fitmodel.R")