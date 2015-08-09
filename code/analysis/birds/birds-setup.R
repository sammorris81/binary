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
these <- sample(x = length(cattle_egret), size = 10000)
s_sub <- s[these, ]

cattle_egret_sub <- cattle_egret[these]
plot(s_sub, type = "n")
points(s_sub[!cattle_egret_sub, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s_sub[cattle_egret_sub, ], pch = 21, col = "firebrick4", bg = "firebrick1")

sanderling_sub <- sanderling[these]
plot(s_sub, type = "n")
points(s_sub[!sanderling_sub, ], pch = 21, col = "dodgerblue4", bg = "dodgerblue1")
points(s_sub[sanderling_sub, ], pch = 21, col = "firebrick4", bg = "firebrick1")

theme_clean <- function(base_size = 12) {
  require(grid)
  theme_grey(base_size) %+replace%
    theme(
      axis.title      =   element_blank(),
      axis.text       =   element_blank(),
      panel.background    =   element_blank(),
      panel.grid      =   element_blank(),
      axis.ticks.length   =   unit(0,"cm"),
      axis.ticks.margin   =   unit(0,"cm"),
      panel.margin    =   unit(0,"lines"),
      plot.margin     =   unit(c(0,0,0,0),"lines"),
      complete = TRUE
    )
}

us <- map_data(map = "state")
p <- ggplot()
p <- p + geom_polygon(data = us, aes(x = long, y = lat, group = group), 
                      colour = "grey35", fill = "white")
p <- p + theme_clean()
p
