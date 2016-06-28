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
dev.print(device = png, width = 10, height = 7,
          file = "../../../LaTeX/plots/cattle_egret.png")

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

map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = "Vesper sparrow sightings in 2002", cex.main = 2)
points(s[!vesper_sparrow, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[vesper_sparrow, ], pch = 21, col = "firebrick4", bg = "firebrick1")
dev.print(device = png, width = 10, height = 7,
          file = "../../../LaTeX/plots/vesper_sparrow.png")

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
plot(s_sub, type = 'n', main = "Subsample of 2012 Sightings",
     axes = FALSE, xlab = "", ylab = "", cex.main = 2)
points(s_sub[!cattle_egret_sub, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s_sub[cattle_egret_sub, ], pch = 21, col = "firebrick4", bg = "firebrick1")
dev.print(device = png, width = 10, height = 7,
          file = "../../../LaTeX/plots/cattle_egret.png")

sanderling_sub <- sanderling[these]
# map("state", xlim = range(s[, 1]), ylim = range(c(s[, 2], 25.1, 49.4)))
plot(s_sub, type = 'n', main = "Subsample of 2012 Sightings",
     axes = FALSE, xlab = "", ylab = "", cex.main = 2)
points(s_sub[!sanderling_sub, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s_sub[sanderling_sub, ], pch = 21, col = "firebrick4", bg = "firebrick1")
dev.print(device = png, width = 10, height = 7,
          file = "../../../LaTeX/plots/sanderling.png")

#### grid up the bird data ####
library(maps)
library(maptools)
library(ggplot2)
library(fields)
birds <- read.csv("./birds/2002/checklists.csv")
s <- cbind(birds$LONGITUDE, birds$LATITUDE)

# s is the nx2 matrix of lat/longs of the e-birds obs
# y is the n-vector of e-birds counts
# sg is the ng x 2 lat/longs of the grid centers
# Y is the counts in each grid cell

xg    <- seq(-124.75, -66.75, 0.25)
yg    <- seq( 24.50,  49.00, 0.25)

sg    <- as.matrix(expand.grid(xg, yg))
inusa <- !is.na(map.where("usa", x = sg[, 1], y = sg[, 2]))
sg    <- sg[inusa, ]
ng    <- nrow(s)

d <- rdist(s, sg)
mind  <- apply(d, 1, min)        # gives shortest distance to grid cell
cell  <- apply(d, 1, which.min)  # gives which grid cell is the closest
save(mind, cell, sg, xg, yg, file = "gridedUS.RData")

rm(list = ls())
library(maps)
library(maptools)
library(ggplot2)
library(fields)
birds <- read.csv("./birds/2002/checklists.csv")
s <- cbind(birds$LONGITUDE, birds$LATITUDE)
load("gridedUS.RData")

cattle_egret     <- birds$Chordeiles_minor != "0"
common_nighthawk <- birds$Bubulcus_ibis != "0"
western_bluebird <- birds$Sialia_mexicana != "0"
vesper_sparrow   <- birds$Pooecetes_gramineus != "0"
# our analysis is conditional on there being an observation in a cell
# in the grid vector, we're using
# 2: no observations in cell
# 1: observations exist, and species was seen
# 0: observations exist, but species not seen
cattle_egret_grid <- rep(2, nrow(sg))
cattle_egret_grid[unique(cell[!cattle_egret])] <- 0
cattle_egret_grid[unique(cell[cattle_egret])]  <- 1

common_nighthawk_grid <- rep(2, nrow(sg))
common_nighthawk_grid[unique(cell[!common_nighthawk])] <- 0
common_nighthawk_grid[unique(cell[common_nighthawk])]  <- 1

western_bluebird_grid <- rep(2, nrow(sg))
western_bluebird_grid[unique(cell[!western_bluebird])] <- 0
western_bluebird_grid[unique(cell[western_bluebird])]  <- 1

vesper_sparrow_grid <- rep(2, nrow(sg))
vesper_sparrow_grid[unique(cell[!vesper_sparrow])] <- 0
vesper_sparrow_grid[unique(cell[vesper_sparrow])]  <- 1

# us_map <- map("state")
dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = "Actual Cattle egret sightings in 2002", cex.main = 2)
points(s[!cattle_egret, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[cattle_egret, ], pch = 21, col = "firebrick4", bg = "firebrick1")
dev.print(device = pdf, width = 10, height = 7, 
          file = "plots/actual_cattle_egret.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = "Gridded Cattle egret sightings in 2002", cex.main = 2)
points(sg[cattle_egret_grid == 0, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(sg[cattle_egret_grid == 1, ], pch = 21, col = "firebrick4", bg = "firebrick1")
dev.print(device = pdf, width = 10, height = 7, 
          file = "plots/grid_cattle_egret.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = "Actual Common nighthawk sightings in 2002", cex.main = 2)
points(s[!common_nighthawk, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[common_nighthawk, ], pch = 21, col = "firebrick4", bg = "firebrick1")
dev.print(device = pdf, width = 10, height = 7, 
          file = "plots/actual_common_nighthawk.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = "Gridded Common nighthawk sightings in 2002", cex.main = 2)
points(sg[common_nighthawk_grid == 0, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(sg[common_nighthawk_grid == 1, ], pch = 21, col = "firebrick4", bg = "firebrick1")
dev.print(device = pdf, width = 10, height = 7, 
          file = "plots/grid_common_nighthawk.pdf")
dev.off()