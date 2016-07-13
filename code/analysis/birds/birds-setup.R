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

seen <- birds[, 18:ncol(birds)] != 0
rate <- apply(seen, 2, sum)
rate <- rate[rate != 0]


sort(rate[which(rate < 200 & rate > 50)])
# Asio_otus - 91 - longeared_owl
# Charadrius_nivosus - 191
snowy_plover <- birds$Charadrius_nivosus != "0"
longeared_owl <- birds$Asio_otus != "0"

sort(rate[which(rate < 400 & rate > 250)])
# Charadrius_melodus - 201 (Piping plover)
# Icterus_cucullatus - 289 (Hooded oriole)
piping_plover <- birds$Charadrius_melodus != "0"
hooded_oriole <- birds$Icterus_cucullatus != "0"


sort(rate[which(rate < 1000 & rate > 500)])  # between about 1.02% and 2.04%
# Numenius_americanus - 501
# Columbina_passerina - 512
# Sialia_currucoides  - 507
# Anser_albifrons     - 544
# Chordeiles_minor    - 669

# Somewhere around 1% and 2%
cattle_egret     <- birds$Chordeiles_minor != "0"
common_nighthawk <- birds$Bubulcus_ibis != "0"
western_bluebird <- birds$Sialia_mexicana != "0"
vesper_sparrow   <- birds$Pooecetes_gramineus != "0"
longbilled_curlew <- birds$Numenius_americanus != "0"
common_grounddove <- birds$Columbina_passerina != "0"
mountain_bluebird <- birds$Sialia_currucoides != "0"
greater_white_goose <- birds$Anser_albifrons != "0"

sort(rate[which(rate < 2500 & rate > 2000)])
# Anas_discors  - 2008
# Vireo_griseus - 2131
# Accipter_striatus - 2237
# Spinus_psaltria - 2429

# Somewhere around 5%
bluewinged_teal <- birds$Anas_discors != "0"
whiteeyed_vireo <- birds$Vireo_griseus != "0"
sharpshinned_hawk <- birds$Accipiter_striatus != "0"
lesser_goldfinch <- birds$Spinus_psaltria != "0"

# our analysis is conditional on there being an observation in a cell
# in the grid vector, we're using
# 2: no observations in cell
# 1: observations exist, and species was seen
# 0: observations exist, but species not seen
cattle_egret_grid <- rep(2, nrow(sg))
cattle_egret_grid_n <- rep(0, nrow(sg))
cattle_egret_grid_y <- rep(0, nrow(sg))
cattle_egret_grid[unique(cell[!cattle_egret])] <- 0
cattle_egret_grid[unique(cell[cattle_egret])]  <- 1
for (i in 1:length(cattle_egret)) {
  cattle_egret_grid_n[cell[i]] <- cattle_egret_grid_n[cell[i]] + 1
  if (cattle_egret[i]) {
    cattle_egret_grid_y[cell[i]] <- cattle_egret_grid_y[cell[i]] + 1
  }
}
cattle_egret_grid_n <- cattle_egret_grid_n[cattle_egret_grid != 2]
cattle_egret_grid_y <- cattle_egret_grid_y[cattle_egret_grid != 2]
# phats <- (cattle_egret_grid_y / cattle_egret_grid_n)[cattle_egret_grid != 2]
# cattle_egret_grid_n[cattle_egret_grid != 2][3456]
# cattle_egret_grid_y[cattle_egret_grid != 2][3336]
# hist(phats)

common_nighthawk_grid <- rep(2, nrow(sg))
common_nighthawk_grid_n <- rep(2, nrow(sg))
common_nighthawk_grid_y <- rep(2, nrow(sg))
common_nighthawk_grid[unique(cell[!common_nighthawk])] <- 0
common_nighthawk_grid[unique(cell[common_nighthawk])]  <- 1
for (i in 1:length(common_nighthawk)) {
  common_nighthawk_grid_n[cell[i]] <- common_nighthawk_n[cell[i]] + 1
  if (common_nighthawk[i]) {
    common_nighthawk_grid_y[cell[i]] <- common_nighthawk_grid_y[cell[i]] + 1
  }
}
common_nighthawk_grid_n <- common_nighthawk_grid_n[common_nighthawk_grid != 2]
common_nighthawk_grid_y <- common_nighthawk_grid_y[common_nighthawk_grid != 2]

western_bluebird_grid <- rep(2, nrow(sg))
western_bluebird_grid_n <- rep(0, nrow(sg))
western_bluebird_grid_y <- rep(0, nrow(sg))
western_bluebird_grid[unique(cell[!western_bluebird])] <- 0
western_bluebird_grid[unique(cell[western_bluebird])]  <- 1
for (i in 1:length(western_bluebird)) {
  western_bluebird_grid_n[cell[i]] <- western_bluebird_grid_n[cell[i]] + 1
  if (western_bluebird[i]) {
    western_bluebird_grid_y[cell[i]] <- western_bluebird_grid_y[cell[i]] + 1
  }
}
western_bluebird_grid_n <- western_bluebird_grid_n[western_bluebird_grid != 2]
western_bluebird_grid_y <- western_bluebird_grid_y[western_bluebird_grid != 2]

vesper_sparrow_grid <- rep(2, nrow(sg))
vesper_sparrow_grid[unique(cell[!vesper_sparrow])] <- 0
vesper_sparrow_grid[unique(cell[vesper_sparrow])]  <- 1
vesper_sparrow_grid_n <- rep(0, nrow(sg))
vesper_sparrow_grid_y <- rep(0, nrow(sg))
for (i in 1:length(vesper_sparrow)) {
  vesper_sparrow_grid_n[cell[i]] <- vesper_sparrow_grid_n[cell[i]] + 1
  if (vesper_sparrow[i]) {
    vesper_sparrow_grid_y[cell[i]] <- vesper_sparrow_grid_y[cell[i]] + 1
  }
}
vesper_sparrow_grid_n <- vesper_sparrow_grid_n[vesper_sparrow_grid != 2]
vesper_sparrow_grid_y <- vesper_sparrow_grid_y[vesper_sparrow_grid != 2]

bluewinged_teal_grid <- rep(2, nrow(sg))
bluewinged_teal_grid_n <- rep(0, nrow(sg))
bluewinged_teal_grid_y <- rep(0, nrow(sg))
bluewinged_teal_grid[unique(cell[!bluewinged_teal])] <- 0
bluewinged_teal_grid[unique(cell[bluewinged_teal])]  <- 1
for (i in 1:length(bluewinged_teal)) {
  bluewinged_teal_grid_n[cell[i]] <- bluewinged_teal_grid_n[cell[i]] + 1
  if (bluewinged_teal[i]) {
    bluewinged_teal_grid_y[cell[i]] <- bluewinged_teal_grid_y[cell[i]] + 1
  }
}
bluewinged_teal_grid_n <- bluewinged_teal_grid_n[bluewinged_teal_grid != 2]
bluewinged_teal_grid_y <- bluewinged_teal_grid_y[bluewinged_teal_grid != 2]

whiteeyed_vireo_grid <- rep(2, nrow(sg))
whiteeyed_vireo_grid_n <- rep(0, nrow(sg))
whiteeyed_vireo_grid_y <- rep(0, nrow(sg))
whiteeyed_vireo_grid[unique(cell[!whiteeyed_vireo])] <- 0
whiteeyed_vireo_grid[unique(cell[whiteeyed_vireo])]  <- 1
for (i in 1:length(whiteeyed_vireo)) {
  whiteeyed_vireo_grid_n[cell[i]] <- whiteeyed_vireo_grid_n[cell[i]] + 1
  if (whiteeyed_vireo[i]) {
    whiteeyed_vireo_grid_y[cell[i]] <- whiteeyed_vireo_grid_y[cell[i]] + 1
  }
}
whiteeyed_vireo_grid_n <- whiteeyed_vireo_grid_n[whiteeyed_vireo_grid != 2]
whiteeyed_vireo_grid_y <- whiteeyed_vireo_grid_y[whiteeyed_vireo_grid != 2]

sharpshinned_hawk_grid <- rep(2, nrow(sg))
sharpshinned_hawk_grid_n <- rep(0, nrow(sg))
sharpshinned_hawk_grid_y <- rep(0, nrow(sg))
sharpshinned_hawk_grid[unique(cell[!sharpshinned_hawk])] <- 0
sharpshinned_hawk_grid[unique(cell[sharpshinned_hawk])]  <- 1
for (i in 1:length(sharpshinned_hawk)) {
  sharpshinned_hawk_grid_n[cell[i]] <- sharpshinned_hawk_grid_n[cell[i]] + 1
  if (sharpshinned_hawk[i]) {
    sharpshinned_hawk_grid_y[cell[i]] <- sharpshinned_hawk_grid_y[cell[i]] + 1
  }
}
sharpshinned_hawk_grid_n <- sharpshinned_hawk_grid_n[sharpshinned_hawk_grid != 2]
sharpshinned_hawk_grid_y <- sharpshinned_hawk_grid_y[sharpshinned_hawk_grid != 2]

lesser_goldfinch_grid <- rep(2, nrow(sg))
lesser_goldfinch_grid_n <- rep(0, nrow(sg))
lesser_goldfinch_grid_y <- rep(0, nrow(sg))
lesser_goldfinch_grid[unique(cell[!lesser_goldfinch])] <- 0
lesser_goldfinch_grid[unique(cell[lesser_goldfinch])]  <- 1
for (i in 1:length(lesser_goldfinch)) {
  lesser_goldfinch_grid_n[cell[i]] <- lesser_goldfinch_grid_n[cell[i]] + 1
  if (lesser_goldfinch[i]) {
    lesser_goldfinch_grid_y[cell[i]] <- lesser_goldfinch_grid_y[cell[i]] + 1
  }
}
lesser_goldfinch_grid_n <- lesser_goldfinch_grid_n[lesser_goldfinch_grid != 2]
lesser_goldfinch_grid_y <- lesser_goldfinch_grid_y[lesser_goldfinch_grid != 2]

longbilled_curlew_grid <- rep(2, nrow(sg))
longbilled_curlew_grid_n <- rep(0, nrow(sg))
longbilled_curlew_grid_y <- rep(0, nrow(sg))
longbilled_curlew_grid[unique(cell[!longbilled_curlew])] <- 0
longbilled_curlew_grid[unique(cell[longbilled_curlew])]  <- 1
for (i in 1:length(longbilled_curlew)) {
  longbilled_curlew_grid_n[cell[i]] <- longbilled_curlew_grid_n[cell[i]] + 1
  if (longbilled_curlew[i]) {
    longbilled_curlew_grid_y[cell[i]] <- longbilled_curlew_grid_y[cell[i]] + 1
  }
}
longbilled_curlew_grid_n <- longbilled_curlew_grid_n[longbilled_curlew_grid != 2]
longbilled_curlew_grid_y <- longbilled_curlew_grid_y[longbilled_curlew_grid != 2]

common_grounddove_grid <- rep(2, nrow(sg))
common_grounddove_grid_n <- rep(0, nrow(sg))
common_grounddove_grid_y <- rep(0, nrow(sg))
common_grounddove_grid[unique(cell[!common_grounddove])] <- 0
common_grounddove_grid[unique(cell[common_grounddove])]  <- 1
for (i in 1:length(common_grounddove)) {
  common_grounddove_grid_n[cell[i]] <- common_grounddove_grid_n[cell[i]] + 1
  if (common_grounddove[i]) {
    common_grounddove_grid_y[cell[i]] <- common_grounddove_grid_y[cell[i]] + 1
  }
}
common_grounddove_grid_n <- common_grounddove_grid_n[common_grounddove_grid != 2]
common_grounddove_grid_y <- common_grounddove_grid_y[common_grounddove_grid != 2]

mountain_bluebird_grid <- rep(2, nrow(sg))
mountain_bluebird_grid_n <- rep(0, nrow(sg))
mountain_bluebird_grid_y <- rep(0, nrow(sg))
mountain_bluebird_grid[unique(cell[!mountain_bluebird])] <- 0
mountain_bluebird_grid[unique(cell[mountain_bluebird])]  <- 1
for (i in 1:length(mountain_bluebird)) {
  mountain_bluebird_grid_n[cell[i]] <- mountain_bluebird_grid_n[cell[i]] + 1
  if (mountain_bluebird[i]) {
    mountain_bluebird_grid_y[cell[i]] <- mountain_bluebird_grid_y[cell[i]] + 1
  }
}
mountain_bluebird_grid_n <- mountain_bluebird_grid_n[mountain_bluebird_grid != 2]
mountain_bluebird_grid_y <- mountain_bluebird_grid_y[mountain_bluebird_grid != 2]

greater_white_goose_grid <- rep(2, nrow(sg))
greater_white_goose_grid_n <- rep(0, nrow(sg))
greater_white_goose_grid_y <- rep(0, nrow(sg))
greater_white_goose_grid[unique(cell[!greater_white_goose])] <- 0
greater_white_goose_grid[unique(cell[greater_white_goose])]  <- 1
for (i in 1:length(greater_white_goose)) {
  greater_white_goose_grid_n[cell[i]] <- greater_white_goose_grid_n[cell[i]] + 1
  if (greater_white_goose[i]) {
    greater_white_goose_grid_y[cell[i]] <- greater_white_goose_grid_y[cell[i]] + 1
  }
}
greater_white_goose_grid_n <- greater_white_goose_grid_n[greater_white_goose_grid != 2]
greater_white_goose_grid_y <- greater_white_goose_grid_y[greater_white_goose_grid != 2]

snowy_plover_grid <- rep(2, nrow(sg))
snowy_plover_grid_n <- rep(0, nrow(sg))
snowy_plover_grid_y <- rep(0, nrow(sg))
snowy_plover_grid[unique(cell[!snowy_plover])] <- 0
snowy_plover_grid[unique(cell[snowy_plover])]  <- 1
for (i in 1:length(snowy_plover)) {
  snowy_plover_grid_n[cell[i]] <- snowy_plover_grid_n[cell[i]] + 1
  if (snowy_plover[i]) {
    snowy_plover_grid_y[cell[i]] <- snowy_plover_grid_y[cell[i]] + 1
  }
}
snowy_plover_grid_n <- snowy_plover_grid_n[snowy_plover_grid != 2]
snowy_plover_grid_y <- snowy_plover_grid_y[snowy_plover_grid != 2]

longeared_owl_grid <- rep(2, nrow(sg))
longeared_owl_grid_n <- rep(0, nrow(sg))
longeared_owl_grid_y <- rep(0, nrow(sg))
longeared_owl_grid[unique(cell[!longeared_owl])] <- 0
longeared_owl_grid[unique(cell[longeared_owl])]  <- 1
for (i in 1:length(longeared_owl)) {
  longeared_owl_grid_n[cell[i]] <- longeared_owl_grid_n[cell[i]] + 1
  if (longeared_owl[i]) {
    longeared_owl_grid_y[cell[i]] <- longeared_owl_grid_y[cell[i]] + 1
  }
}
longeared_owl_grid_n <- longeared_owl_grid_n[longeared_owl_grid != 2]
longeared_owl_grid_y <- longeared_owl_grid_y[longeared_owl_grid != 2]

piping_plover_grid <- rep(2, nrow(sg))
piping_plover_grid_n <- rep(0, nrow(sg))
piping_plover_grid_y <- rep(0, nrow(sg))
piping_plover_grid[unique(cell[!piping_plover])] <- 0
piping_plover_grid[unique(cell[piping_plover])]  <- 1
for (i in 1:length(piping_plover)) {
  piping_plover_grid_n[cell[i]] <- piping_plover_grid_n[cell[i]] + 1
  if (piping_plover[i]) {
    piping_plover_grid_y[cell[i]] <- piping_plover_grid_y[cell[i]] + 1
  }
}
piping_plover_grid_n <- piping_plover_grid_n[piping_plover_grid != 2]
piping_plover_grid_y <- piping_plover_grid_y[piping_plover_grid != 2]

hooded_oriole_grid <- rep(2, nrow(sg))
hooded_oriole_grid_n <- rep(0, nrow(sg))
hooded_oriole_grid_y <- rep(0, nrow(sg))
hooded_oriole_grid[unique(cell[!hooded_oriole])] <- 0
hooded_oriole_grid[unique(cell[hooded_oriole])]  <- 1
for (i in 1:length(hooded_oriole)) {
  hooded_oriole_grid_n[cell[i]] <- hooded_oriole_grid_n[cell[i]] + 1
  if (hooded_oriole[i]) {
    hooded_oriole_grid_y[cell[i]] <- hooded_oriole_grid_y[cell[i]] + 1
  }
}
hooded_oriole_grid_n <- hooded_oriole_grid_n[hooded_oriole_grid != 2]
hooded_oriole_grid_y <- hooded_oriole_grid_y[hooded_oriole_grid != 2]

us_map <- map("state")

this.species <- cattle_egret
this.species.grid <- cattle_egret_grid
this.species.name <- "Cattle Egret"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/cattle_egret_actual.pdf"
filename.grid <- "plots/cattle_egret_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- common_nighthawk
this.species.grid <- common_nighthawk_grid
this.species.name <- "Common Nighthawk"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/common_nighthawk_actual.pdf"
filename.grid <- "plots/common_nighthawk_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- western_bluebird
this.species.grid <- western_bluebird_grid
this.species.name <- "Western Bluebird"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/western_bluebird_actual.pdf"
filename.grid <- "plots/western_bluebird_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- vesper_sparrow
this.species.grid <- vesper_sparrow_grid
this.species.name <- "Vesper Sparrow"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/vesper_sparrow_actual.pdf"
filename.grid <- "plots/vesper_sparrow_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- bluewinged_teal
this.species.grid <- bluewinged_teal_grid
this.species.name <- "Blue-winged Teal"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/bluewinged_teal_actual.pdf"
filename.grid <- "plots/bluewinged_teal_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- whiteeyed_vireo
this.species.grid <- whiteeyed_vireo_grid
this.species.name <- "White-eyed Vireo"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/whiteeyed_vireo_actual.pdf"
filename.grid <- "plots/whiteeyed_vireo_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- sharpshinned_hawk
this.species.grid <- sharpshinned_hawk_grid
this.species.name <- "Sharp-shinned Hawk"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/sharpshinned_hawk_actual.pdf"
filename.grid <- "plots/sharpshinned_hawk_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- sharpshinned_hawk
this.species.grid <- sharpshinned_hawk_grid
this.species.name <- "Sharp-shinned Hawk"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/sharpshinned_hawk_actual.pdf"
filename.grid <- "plots/sharpshinned_hawk_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- lesser_goldfinch
this.species.grid <- lesser_goldfinch_grid
this.species.name <- "Lesser Goldfinch"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/lesser_goldfinch_actual.pdf"
filename.grid <- "plots/lesser_goldfinch_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- longbilled_curlew
this.species.grid <- longbilled_curlew_grid
this.species.name <- "Long-billed Curlew"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/longbilled_curlew_actual.pdf"
filename.grid <- "plots/longbilled_curlew_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- common_grounddove
this.species.grid <- common_grounddove_grid
this.species.name <- "Common Ground-Dove"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/common_grounddove_actual.pdf"
filename.grid <- "plots/common_grounddove_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- mountain_bluebird
this.species.grid <- mountain_bluebird_grid
this.species.name <- "Mountain Bluebird"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/mountain_bluebird_actual.pdf"
filename.grid <- "plots/mountain_bluebird_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- greater_white_goose
this.species.grid <- greater_white_goose_grid
this.species.name <- "Greater White-fronted Goose"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/greater_white_goose_actual.pdf"
filename.grid <- "plots/greater_white_goose_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- snowy_plover
this.species.grid <- snowy_plover_grid
this.species.name <- "Snowy plover"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/snowy_plover_actual.pdf"
filename.grid <- "plots/snowy_plover_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- longeared_owl
this.species.grid <- longeared_owl_grid
this.species.name <- "Long-eared owl"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/longeared_owl_actual.pdf"
filename.grid <- "plots/longeared_owl_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- piping_plover
this.species.grid <- piping_plover_grid
this.species.name <- "Piping plover"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/piping_plover_actual.pdf"
filename.grid <- "plots/piping_plover_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
dev.off()

this.species <- hooded_oriole
this.species.grid <- hooded_oriole_grid
this.species.name <- "Hooded oriole"
this.title.main.act <- paste("Actual", this.species.name, "sightings")
this.title.main.grid <- paste("Gridded", this.species.name, "sightings")
filename.act <- "plots/hooded_oriole_actual.pdf"
filename.grid <- "plots/hooded_oriole_grid.pdf"

dev.new()
map("state",
    xlim = range(c(s[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(s[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.act, cex.main = 2,
      sub = paste(round(mean(this.species) * 100, 2), "%", sep = ""))
points(s[!this.species, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s[this.species, ], pch = 21, col = "firebrick4", bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.act)
dev.off()

dev.new()
map("state",
    xlim = range(c(sg[, 1], us_map$x), na.rm = TRUE),
    ylim = range(c(sg[, 2], us_map$y), na.rm = TRUE))
title(main = this.title.main.grid, cex.main = 2,
      sub = paste(round(mean(this.species.grid == 1) * 100, 2), "%", sep = ""))
points(sg[this.species.grid == 0, ], pch = 21, col = "dodgerblue4",
       bg = "dodgerblue1")
points(sg[this.species.grid == 1, ], pch = 21, col = "firebrick4",
       bg = "firebrick1")
map("state", add = TRUE)
dev.print(device = pdf, file = filename.grid)
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
knots.1 <- knots.2 <- knots.5 <- vector(mode = "list", length = 2)
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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(cattle_egret, cattle_egret_grid_n, cattle_egret_grid_y,
     s, knots.1, knots.2, knots.5, knots.10, knots.15, knots.20,
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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(common_nighthawk, common_nighthawk_grid_n, common_nighthawk_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "common_nighthawk.RData")

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(western_bluebird, western_bluebird_grid_n, western_bluebird_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "western_bluebird.RData")

#### vesper_sparrow: 10.32% ####
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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(vesper_sparrow, vesper_sparrow_grid_n, vesper_sparrow_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15, knots.20,
     cv.idx, file = "vesper_sparrow.RData")

#### Long-billed Curlew ####
longbilled_curlew <- longbilled_curlew_grid[longbilled_curlew_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(longbilled_curlew == 1) / ns
these.ones  <- which(longbilled_curlew == 1)
these.zeros <- which(longbilled_curlew == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(longbilled_curlew, longbilled_curlew_grid_n, longbilled_curlew_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "longbilled_curlew.RData")

#### Common Ground Dove ####
common_grounddove <- common_grounddove_grid[common_grounddove_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(common_grounddove == 1) / ns
these.ones  <- which(common_grounddove == 1)
these.zeros <- which(common_grounddove == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(common_grounddove, common_grounddove_grid_n, common_grounddove_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "common_grounddove.RData")

#### Mountain bluebird ####
mountain_bluebird <- mountain_bluebird_grid[mountain_bluebird_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(mountain_bluebird == 1) / ns
these.ones  <- which(mountain_bluebird == 1)
these.zeros <- which(mountain_bluebird == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(mountain_bluebird, mountain_bluebird_grid_n, mountain_bluebird_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "mountain_bluebird.RData")

#### Greater White-fronted Goose ####
greater_white_goose <- greater_white_goose_grid[greater_white_goose_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(greater_white_goose == 1) / ns
these.ones  <- which(greater_white_goose == 1)
these.zeros <- which(greater_white_goose == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(greater_white_goose, greater_white_goose_grid_n, greater_white_goose_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "greater_white_goose.RData")

#### Blue-winged Teal ####
bluewinged_teal <- bluewinged_teal_grid[bluewinged_teal_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(bluewinged_teal == 1) / ns
these.ones  <- which(bluewinged_teal == 1)
these.zeros <- which(bluewinged_teal == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(bluewinged_teal, bluewinged_teal_grid_n, bluewinged_teal_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "bluewinged_teal.RData")

#### White-eyed Vireo ####
whiteeyed_vireo <- whiteeyed_vireo_grid[whiteeyed_vireo_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(whiteeyed_vireo == 1) / ns
these.ones  <- which(whiteeyed_vireo == 1)
these.zeros <- which(whiteeyed_vireo == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(whiteeyed_vireo, whiteeyed_vireo_grid_n, whiteeyed_vireo_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "whiteeyed_vireo.RData")

#### Sharp-shinned Hawk ####
sharpshinned_hawk <- sharpshinned_hawk_grid[sharpshinned_hawk_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(sharpshinned_hawk == 1) / ns
these.ones  <- which(sharpshinned_hawk == 1)
these.zeros <- which(sharpshinned_hawk == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(sharpshinned_hawk, sharpshinned_hawk_grid_n, sharpshinned_hawk_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "sharpshinned_hawk.RData")

#### Lesser Goldfinch ####
lesser_goldfinch <- lesser_goldfinch_grid[lesser_goldfinch_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(lesser_goldfinch == 1) / ns
these.ones  <- which(lesser_goldfinch == 1)
these.zeros <- which(lesser_goldfinch == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(lesser_goldfinch, lesser_goldfinch_grid_n, lesser_goldfinch_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15,
     knots.20, cv.idx, file = "lesser_goldfinch.RData")

#### Snowy plover ####
snowy_plover <- snowy_plover_grid[snowy_plover_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(snowy_plover == 1) / ns
these.ones  <- which(snowy_plover == 1)
these.zeros <- which(snowy_plover == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(snowy_plover, snowy_plover_grid_n, snowy_plover_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15, knots.20,
     cv.idx, file = "snowy_plover.RData")

#### Long-eared owl ####
longeared_owl <- longeared_owl_grid[longeared_owl_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(longeared_owl == 1) / ns
these.ones  <- which(longeared_owl == 1)
these.zeros <- which(longeared_owl == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(longeared_owl, longeared_owl_grid_n, longeared_owl_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15, knots.20,
     cv.idx, file = "longeared_owl.RData")

#### Piping plover ####
piping_plover <- piping_plover_grid[piping_plover_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(piping_plover == 1) / ns
these.ones  <- which(piping_plover == 1)
these.zeros <- which(piping_plover == 0)

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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(piping_plover, piping_plover_grid_n, piping_plover_grid_y, 
     s, knots.1, knots.2, knots.5, knots.10, knots.15, knots.20,
     cv.idx, file = "piping_plover.RData")

#### Hooded oriole ####
hooded_oriole <- hooded_oriole_grid[hooded_oriole_grid != 2]
cv.idx  <- vector(mode = "list", length = 2)
prop.ones <- sum(hooded_oriole == 1) / ns
these.ones  <- which(hooded_oriole == 1)
these.zeros <- which(hooded_oriole == 0)

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
knots.1 <- knots.2 <- knots.5 <- vector(mode = "list", length = 2)
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

# 1, 2, and 5 are after 10, 15, 20 because original cv only had 10, 15, 20 in it
# Knots at 1% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.01)
  knots.1[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 2% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.02)
  knots.2[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

# Knots at 5% of the sites
for (i in 1:2) {
  nknots <- floor(length(cv.idx[[i]]) * 0.05)
  knots.5[[i]] <- cover.design(R = s[cv.idx[[i]], ], nd = nknots)$design
}

save(hooded_oriole, hooded_oriole_grid_n, hooded_oriole_grid_y,
     s, knots.1, knots.2, knots.5, knots.10, knots.15, knots.20,
     cv.idx, file = "hooded_oriole.RData")