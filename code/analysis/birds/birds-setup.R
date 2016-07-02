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

#### Look at the gridded data in comparison to actual data
rm(list = ls())
library(maps)
library(maptools)
library(ggplot2)
library(fields)
birds <- read.csv("./birds/2002/checklists.csv")
s <- cbind(birds$LONGITUDE, birds$LATITUDE)
load("griddedUS.RData")

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

us_map <- map("state")
dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = "Actual Cattle egret sightings in 2002", cex.main = 2)
points(s[!cattle_egret, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[cattle_egret, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = "plots/cattle_egret_actual.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = "Gridded Cattle egret sightings in 2002", cex.main = 2)
points(sg[cattle_egret_grid == 0, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(sg[cattle_egret_grid == 1, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = "plots/cattle_egret_grid.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = "Actual Common nighthawk sightings in 2002", cex.main = 2)
points(s[!common_nighthawk, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[common_nighthawk, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = "plots/common_nighthawk_actual.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = "Gridded Common nighthawk sightings in 2002", cex.main = 2)
points(sg[common_nighthawk_grid == 0, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(sg[common_nighthawk_grid == 1, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = "plots/common_nighthawk_grid.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = "Actual Western bluebird sightings in 2002", cex.main = 2)
points(s[!western_bluebird, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[western_bluebird, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = "plots/western_bluebird_actual.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = "Gridded Western bluebird sightings in 2002", cex.main = 2)
points(sg[western_bluebird_grid == 0, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(sg[western_bluebird_grid == 1, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = "plots/western_bluebird_grid.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = "Actual Vesper sparrow sightings in 2002", cex.main = 2)
points(s[!vesper_sparrow, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[vesper_sparrow, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = "plots/vesper_sparrow_actual.pdf")
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = "Gridded Vesper sparrow sightings in 2002", cex.main = 2)
points(sg[vesper_sparrow_grid == 0, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(sg[vesper_sparrow_grid == 1, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = "plots/vesper_sparrow_grid.pdf")
dev.off()

# For the knot selection in the simulation study, we opted for a 21 x 21 grid. 
# Because the ebirds data are not quite as homogeneous as our simulated data, 
# I'd like to opt for using a space-filling design with knots at 10%, 15% and 
# 20% of the sites. Because of some concerns with how long this might take to 
# run (the simulation study took around 5 hours to fit our model at 1300 sites), 
# I'm leaning toward 2-fold cross-validation where the holdout sample is a 
# stratified sample from the original dataset matching P(Y = 1). So for example, 
# if in the actual dataset, we have Y = 1 at 8% of the grid cells, then both the 
# testing and training will aim for around 8% 1s.

#### set up cross-validation ####
# we need separate cross-validation sets for each species because we want to 
# use a stratified sample to make sure that the proportion of 1s in the 
# testing set is similar to the training set

# first find the count of grid cells that have observations - should be the same
# for all species
s <- sg[cattle_egret_grid != 2, ]
ns <- nrow(s)

#### cattle_egret: 9.92% ####
cattle_egret <- cattle_egret_grid[cattle_egret_grid != 2]

cv.idx      <- vector(mode = "list", length = 2)
prop.ones   <- sum(cattle_egret == 1) / ns
these.ones  <- which(cattle_egret == 1)
these.zeros <- which(cattle_egret == 0)

#### Stratified cross-validation ####
# Making sure P(Y = 1) is similar for test and train
set.seed(28)  # cv
samp.ones   <- sample(these.ones)
samp.zeros  <- sample(these.zeros)
ntrain.ones <- c(ceiling(length(these.ones) / 2), floor(length(these.ones) / 2))
ntrain.zeros <- c(ceiling(length(these.zeros) / 2), 
                  floor(length(these.zeros) / 2))

cv.idx[[1]] <- c(sort(samp.ones[1:ntrain.ones[1]]), 
                 sort(samp.zeros[1:ntrain.zeros[1]]))
cv.idx[[2]] <- (1:ns)[-cv.idx[[1]]]

#### Generate the knot locations
knots.10 <- knots.15 <- knots.20 <- vector(mode = "list", length = 2)
set.seed(5668)

# Knots at 10% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.10)
  knots.10[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 15% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.15)
  knots.15[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 20% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.20)
  knots.20[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(cattle_egret, s, knots.10, knots.15, knots.20, 
     cv.idx, file = "cattle_egret.RData")

#### common_nighthawk: 7.90% ####
common_nighthawk <- common_nighthawk_grid[common_nighthawk_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(common_nighthawk == 1) / ns
these.ones  <- which(common_nighthawk == 1)
these.zeros <- which(common_nighthawk == 0)

#### Stratified cross-validation ####
# Making sure P(Y = 1) is similar for test and train
set.seed(28)  # cv
samp.ones   <- sample(these.ones)
samp.zeros  <- sample(these.zeros)
ntrain.ones <- c(ceiling(length(these.ones) / 2), floor(length(these.ones) / 2))
ntrain.zeros <- c(ceiling(length(these.zeros) / 2), 
                  floor(length(these.zeros) / 2))

cv.idx[[1]] <- c(sort(samp.ones[1:ntrain.ones[1]]), 
                 sort(samp.zeros[1:ntrain.zeros[1]]))
cv.idx[[2]] <- (1:ns)[-cv.idx[[1]]]

#### Generate the knot locations
knots.10 <- knots.15 <- knots.20 <- vector(mode = "list", length = 2)
set.seed(5668)

# Knots at 10% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.10)
  knots.10[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 15% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.15)
  knots.15[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 20% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.20)
  knots.20[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}
save(common_nighthawk, s, knots.10, knots.15, knots.20, 
     cv.idx, file = "common_nighthawk.RData")

#### western_bluebird: 6.34% ####
western_bluebird <- western_bluebird_grid[western_bluebird_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(western_bluebird == 1) / ns
these.ones  <- which(western_bluebird == 1)
these.zeros <- which(western_bluebird == 0)

#### Stratified cross-validation ####
# Making sure P(Y = 1) is similar for test and train
set.seed(28)  # cv
samp.ones   <- sample(these.ones)
samp.zeros  <- sample(these.zeros)
ntrain.ones <- c(ceiling(length(these.ones) / 2), floor(length(these.ones) / 2))
ntrain.zeros <- c(ceiling(length(these.zeros) / 2), 
                  floor(length(these.zeros) / 2))

cv.idx[[1]] <- c(sort(samp.ones[1:ntrain.ones[1]]), 
                 sort(samp.zeros[1:ntrain.zeros[1]]))
cv.idx[[2]] <- (1:ns)[-cv.idx[[1]]]

#### Generate the knot locations
knots.10 <- knots.15 <- knots.20 <- vector(mode = "list", length = 2)
set.seed(5668)

# Knots at 10% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.10)
  knots.10[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 15% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.15)
  knots.15[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 20% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.20)
  knots.20[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}
save(western_bluebird, s, knots.10, knots.15, knots.20, 
     cv.idx, file = "western_bluebird.RData")

# vesper_sparrow: 10.32%
vesper_sparrow <- vesper_sparrow_grid[vesper_sparrow_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(vesper_sparrow == 1) / ns
these.ones  <- which(vesper_sparrow == 1)
these.zeros <- which(vesper_sparrow == 0)

#### Stratified cross-validation ####
# Making sure P(Y = 1) is similar for test and train
set.seed(28)  # cv
samp.ones   <- sample(these.ones)
samp.zeros  <- sample(these.zeros)
ntrain.ones <- c(ceiling(length(these.ones) / 2), floor(length(these.ones) / 2))
ntrain.zeros <- c(ceiling(length(these.zeros) / 2), 
                  floor(length(these.zeros) / 2))

cv.idx[[1]] <- c(sort(samp.ones[1:ntrain.ones[1]]), 
                 sort(samp.zeros[1:ntrain.zeros[1]]))
cv.idx[[2]] <- (1:ns)[-cv.idx[[1]]]

#### Generate the knot locations
knots.10 <- knots.15 <- knots.20 <- vector(mode = "list", length = 2)
set.seed(5668)

# Knots at 10% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.10)
  knots.10[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 15% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.15)
  knots.15[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 20% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.20)
  knots.20[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}
save(vesper_sparrow, s, knots.10, knots.15, knots.20, 
     cv.idx, file = "vesper_sparrow.RData")