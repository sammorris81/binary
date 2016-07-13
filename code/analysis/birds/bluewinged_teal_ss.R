rm(list=ls())
for (set in 1:50) {
  source(file = "./package_load.R", chdir = T)
  species <- "bluewinged_teal"
  n <- 100
  source(file = "./fitmodel_ss.R")
}
