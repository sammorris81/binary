# get the data from the dataset
# Species: Sanderlings and cattle egrets
# Year: 2012
# Source: ebirds data

rm(list=ls())
birds <- read.table("cat_eg_sand_data.txt")
nobs <- nrow(birds)

cattle_egret <- sanderling <- rep(NA, nobs)
s <- cbind(birds$long, birds$lat)

for (i in 1:nrow(birds)) {
  cattle_egret[i] <- birds$cattle_egret[i] == "0"
  sanderling[i] <- birds$sanderling[i] == "0"
}