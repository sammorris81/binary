# species.to.run <- matrix(TRUE, nrow = length(species.list), 
#                          ncol = length(knot.percents) * length(cvs))
# write.table(species.to.run, file = "species_to_run.txt")

rm(list=ls())
source(file = "./package_load.R", chdir = T)

knot.percents <- c(10, 15, 20)
cvs <- c(1, 2)
species.list <- c("bluewinged_teal", "cattle_egret", "common_grounddove", 
                  "common_nighthawk", "greater_white_goose", "lesser_goldfinch",
                  "longbilled_curlew", "mountain_bluebird", "sharpshinned_hawk",
                  "vesper_sparrow", "western_bluebird", "whiteeyed_vireo")

for (species in species.list) {
  for (knot.percent in knot.percents) {
    for (cv in cvs) {
      species.to.run <- read.table(file = "species_to_run.txt")
      this.row <- which(species.list == species)
      this.col <- (cv - 1) * length(knot.percents) + 
        which(knot.percents == knot.percent)
      if (species.to.run[this.row, this.col]) {
        print(paste("Starting: ", species, "-", knot.percent, "-", cv, 
                    sep = ""))
        species.to.run[this.row, this.col] <- FALSE
        write.table(species.to.run, file = "species_to_run.txt")
        source("./change_to_median.R")
      }
    }
  }
}