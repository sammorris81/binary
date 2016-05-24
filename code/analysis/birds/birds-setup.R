library(maps)
library(maptools)
library(ggplot2)
birds <- read.csv("./birds/2002/checklists.csv")
birds$Chordeiles_minor
birds$Bubulcus_ibis
birds$Sialia_mexicana
birds$Pooecetes_gramineus
s <- cbind(birds$LONGITUDE, birds$LATITUDE)

nobs <- nrow(birds)

cattle_egret <- common_nighthawk <- western_bluebird <- vesper_sparrow <- rep(NA, nobs)

for (i in 1:nrow(birds)) {
  cattle_egret[i] <- birds$Chordeiles_minor[i] != "0"
  common_nighthawk[i] <- birds$Bubulcus_ibis[i] != "0"
  western_bluebird[i] <- birds$Sialia_mexicana[i] != "0"
  vesper_sparrow[i] <- birds$Pooecetes_gramineus[i] != "0"
}

quartz(width = 10, height = 7)
us_map <- map("state")
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = "Cattle egret sightings in 2002", cex.main = 2)
points(s[!cattle_egret, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[cattle_egret, ], pch = 21, col = "firebrick4", bg = "firebrick1")
quartz.save(file = "../../../LaTeX/plots/cattle_egret.png", type = "png")
dev.off()

map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
points(s[!common_nighthawk, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[common_nighthawk, ], pch = 21, col = "firebrick4", bg = "firebrick1")

map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
points(s[!western_bluebird, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[western_bluebird, ], pch = 21, col = "firebrick4", bg = "firebrick1")

quartz(width = 10, height = 7)
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = "Vesper sparrow sightings in 2002", cex.main = 2)
points(s[!vesper_sparrow, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[vesper_sparrow, ], pch = 21, col = "firebrick4", bg = "firebrick1")
quartz.save(file = "../../../LaTeX/plots/vesper_sparrow.png", type = "png")
dev.off()

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
  cattle_egret[i] <- birds$cattle_egret[i] != "0"
  sanderling[i] <- birds$sanderling[i] != "0"
}

save(s, cattle_egret, sanderling, file = "birds.RData")

#### EDA
rm(list = ls())
library(maps)
library(maptools)
library(ggplot2)
load(file = "birds.RData")
set.seed(10)
these <- sample(x = length(cattle_egret), size = 150000)
s_sub <- s[these, ]

# remove weird point off the coast
long <- s_sub[, 1]
lat  <- s_sub[, 2]
exclude <- which(s_sub[, 1] > -72 & s_sub[, 1] < -70 &
                 s_sub[, 2] > 37 & s_sub[, 2] < 39, arr.ind = TRUE)
s_sub <- s_sub[-exclude, ]
these <- these[-exclude]

exclude <- which(s_sub[, 1] > -72 & s_sub[, 1] < -71 &
                 s_sub[, 2] > 34.5 & s_sub[, 2] < 35.5, arr.ind = TRUE)
s_sub <- s_sub[-exclude, ]
these <- these[-exclude]

exclude <- which(s_sub[, 1] > -74 & s_sub[, 1] < -73 &
                 s_sub[, 2] > 38 & s_sub[, 2] < 39, arr.ind = TRUE)
s_sub <- s_sub[-exclude, ]
these <- these[-exclude]

exclude <- which(s_sub[, 1] > -73 & s_sub[, 1] < -72 &
                   s_sub[, 2] > 39 & s_sub[, 2] < 40, arr.ind = TRUE)
s_sub <- s_sub[-exclude, ]
these <- these[-exclude]

cattle_egret_sub <- cattle_egret[these]
# map("state", xlim = range(s[, 1]), ylim = range(c(s[, 2], 25.1, 49.4)))
# title(main = "Cattle Egret sighting")
quartz(width = 10, height = 7)
plot(s_sub, type = 'n', main = "Subsample of 2012 Sightings",
     axes = FALSE, xlab = "", ylab = "", cex.main = 2)
points(s_sub[!cattle_egret_sub, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s_sub[cattle_egret_sub, ], pch = 21, col = "firebrick4", bg = "firebrick1")
quartz.save(file = "../../../LaTeX/plots/cattle_egret.png", type = "png")

sanderling_sub <- sanderling[these]
# map("state", xlim = range(s[, 1]), ylim = range(c(s[, 2], 25.1, 49.4)))
plot(s_sub, type = 'n', main = "Subsample of 2012 Sightings",
     axes = FALSE, xlab = "", ylab = "", cex.main = 2)
points(s_sub[!sanderling_sub, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s_sub[sanderling_sub, ], pch = 21, col = "firebrick4", bg = "firebrick1")
quartz.save(file = "../../../LaTeX/plots/sanderling.png", type = "png")

