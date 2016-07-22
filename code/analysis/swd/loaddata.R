#Krishna Pacifici
#July 15, 2013
#get data formatted
library(fields)
dat<-read.csv("SWD.csv")

get_a_var<-function(dat,index){
  D<-dat[,index]
  D<-matrix(D,200,200,byrow=T)
  D<-t(D[200:1,])
}

Y1      <- get_a_var(dat,10)
Y2      <- get_a_var(dat,9)
tci     <- get_a_var(dat,3)
solar_s <- get_a_var(dat,4)   
slope   <- get_a_var(dat,5)
road_di <- get_a_var(dat,6)
elev    <- get_a_var(dat,7)


X1 <-(tci-mean(tci))/sd(as.vector(tci))
X2 <-(solar_s-mean(solar_s))/sd(as.vector(solar_s))
X3 <-(slope-mean(slope))/sd(as.vector(slope))
X4 <-(road_di-mean(road_di))/sd(as.vector(road_di))
X5 <-(elev-mean(elev))/sd(as.vector(elev))


rm(get_a_var,dat,tci,solar_s,slope,road_di,elev)

Y1 <- ifelse(Y1 > 0, 1, 0)  # 2.84% - d_count
Y2 <- ifelse(Y2 > 0, 1, 0)  # 0.54% - e_count

# Y1[1, ] left column
# Y1[200, ] right column
# Y1[, 1] bottom row
# Y1[, 200] top row

#### look to make sure that locations are assigned correctly ####
s1 <- seq(0, 1, length = 200)
s2 <- seq(0, 1, length = 200)
s <- as.matrix(expand.grid(s1, s2))
par(mfrow = c(1, 2))
quilt.plot(x = s[, 1], y = s[, 2], z = as.vector(Y1), nx = 200, ny = 200)
image.plot(Y1)

quilt.plot(x = s[, 1], y = s[, 2], z = as.vector(Y2), nx = 200, ny = 200)

new.s1 <- seq(0, 1, length = 100)
new.s2 <- seq(0, 1, length = 100)
new.s  <- as.matrix(expand.grid(new.s1, new.s2))
newY1 <- rep(0, nrow(new.s))
newY2 <- rep(0, nrow(new.s))
for (i in 1:nrow(s)) {
  this.cell <- which.min(rdist(s[i, , drop = F], new.s))
  newY1[this.cell] <- newY1[this.cell] + Y1[i]
  newY2[this.cell] <- newY2[this.cell] + Y2[i]
}
newY1 <- ifelse(newY1 > 0, 1, 0)
newY2 <- ifelse(newY2 > 0, 1, 0)
quilt.plot(new.s[, 1], new.s[, 2], newY1, nx = 100, ny = 100)

save(Y1, Y2, s, file = "plant_inventory.RData")

# generate the samples ahead of time
source("../../../../usefulR/usefulfunctions.R", chdir = TRUE)
nsets <- 200
n     <- 100

clu.lst.Y1 <- clu.lst.Y2 <- vector(mode = "list", length = nsets)
srs.lst.Y1 <- srs.lst.Y2 <- vector(mode = "list", length = nsets)
for (set in 1:nsets) {
  set.seed(1000 + set)
  nobs <- 0
  while(nobs < 3) {  
    # keep repeating the sampling until there are at least 3 observations
    these.train <- sort(sample(length(Y1), n))
    y.o <- Y1[these.train]
    these.cluster <- y.o == 1
    these.cluster.ids <- these.train[these.cluster]
    for (i in 1:length(these.cluster.ids)) {
      # get the location of the cell where y == 1
      this.cell <- get.arr.idx(these.cluster.ids[i], nrows = 200)
      this.row  <- this.cell[1]
      this.col  <- this.cell[2]
      
      # account for the boundary
      neighbors.row <- c(this.row + 1, this.row - 1)
      neighbors.col <- c(this.col + 1, this.col - 1)
      neighbors.row <- neighbors.row[neighbors.row > 0 & neighbors.row < 200]
      neighbors.col <- neighbors.col[neighbors.col > 0 & neighbors.col < 200]
      for (j in 1:length(neighbors.row)) {
        these.train <- c(these.train,
                         get.idx(row = neighbors.row[j], col = this.col, 
                                 nrows = 200))
      }
      for (j in 1:length(neighbors.col)) {
        these.train <- c(these.train, 
                         get.idx(row = this.row, col = neighbors.col[j], 
                                 nrows = 200))
      }
      
    }
    these.train <- sort(unique(these.train))
    y.o <- Y1[these.train]
    nobs <- sum(y.o)
  }
  clu.lst.Y1[[set]] <- these.train
  
  set.seed(2000 + set)
  nobs <- 0
  nsamp <- length(these.train)  # make sure that we get the same sample size
  while (nobs < 3) {
    these.train <- sort(sample(length(Y1), nsamp))
    y.o <- Y1[these.train]
    nobs <- sum(y.o)
  }
  srs.lst.Y1[[set]] <- these.train
}

for (set in 1:nsets) {
  set.seed(1000 + set)
  nobs <- 0
  while(nobs < 3) {  
    # keep repeating the sampling until there are at least 3 observations
    these.train <- sort(sample(length(Y2), n))
    y.o <- Y2[these.train]
    these.cluster <- y.o == 1
    these.cluster.ids <- these.train[these.cluster]
    for (i in 1:length(these.cluster.ids)) {
      # get the location of the cell where y == 1
      this.cell <- get.arr.idx(these.cluster.ids[i], nrows = 200)
      this.row  <- this.cell[1]
      this.col  <- this.cell[2]
      
      # account for the boundary
      neighbors.row <- c(this.row + 1, this.row - 1)
      neighbors.col <- c(this.col + 1, this.col - 1)
      neighbors.row <- neighbors.row[neighbors.row > 0 & neighbors.row < 200]
      neighbors.col <- neighbors.col[neighbors.col > 0 & neighbors.col < 200]
      for (j in 1:length(neighbors.row)) {
        these.train <- c(these.train,
                         get.idx(row = neighbors.row[j], col = this.col, 
                                 nrows = 200))
      }
      for (j in 1:length(neighbors.col)) {
        these.train <- c(these.train, 
                         get.idx(row = this.row, col = neighbors.col[j], 
                                 nrows = 200))
      }
      
    }
    these.train <- sort(unique(these.train))
    y.o <- Y2[these.train]
    nobs <- sum(y.o)
  }
  clu.lst.Y2[[set]] <- these.train
  
  set.seed(2000 + set)
  nobs <- 0
  nsamp <- length(these.train)  # make sure that we get the same sample size
  while (nobs < 3) {
    these.train <- sort(sample(length(Y2), nsamp))
    y.o <- Y2[these.train]
    nobs <- sum(y.o)
  }
  srs.lst.Y2[[set]] <- these.train
}

save(Y1, Y2, clu.lst.Y1, clu.lst.Y2, srs.lst.Y1, srs.lst.Y2, s, 
     file = "plant_inventory.RData")

load("plant_inventory.RData")
source("../../../../usefulR/usefulfunctions.R", chdir = TRUE)
clu.lst.Y1.100 <- clu.lst.Y1
clu.lst.Y2.100 <- clu.lst.Y2
srs.lst.Y1.100 <- srs.lst.Y1
srs.lst.Y2.100 <- srs.lst.Y2
nsets <- 200
n     <- 250

clu.lst.Y1.250 <- clu.lst.Y2.250 <- vector(mode = "list", length = nsets)
srs.lst.Y1.250 <- srs.lst.Y2.250 <- vector(mode = "list", length = nsets)
for (set in 1:nsets) {
  set.seed(1000 + set)
  nobs <- 0
  while(nobs < 3) {  
    # keep repeating the sampling until there are at least 3 observations
    these.train <- sort(sample(length(Y1), n))
    y.o <- Y1[these.train]
    these.cluster <- y.o == 1
    these.cluster.ids <- these.train[these.cluster]
    for (i in 1:length(these.cluster.ids)) {
      # get the location of the cell where y == 1
      this.cell <- get.arr.idx(these.cluster.ids[i], nrows = 200)
      this.row  <- this.cell[1]
      this.col  <- this.cell[2]
      
      # account for the boundary
      neighbors.row <- c(this.row + 1, this.row - 1)
      neighbors.col <- c(this.col + 1, this.col - 1)
      neighbors.row <- neighbors.row[neighbors.row > 0 & neighbors.row < 200]
      neighbors.col <- neighbors.col[neighbors.col > 0 & neighbors.col < 200]
      for (j in 1:length(neighbors.row)) {
        these.train <- c(these.train,
                         get.idx(row = neighbors.row[j], col = this.col, 
                                 nrows = 200))
      }
      for (j in 1:length(neighbors.col)) {
        these.train <- c(these.train, 
                         get.idx(row = this.row, col = neighbors.col[j], 
                                 nrows = 200))
      }
      
    }
    these.train <- sort(unique(these.train))
    y.o <- Y1[these.train]
    nobs <- sum(y.o)
  }
  clu.lst.Y1.250[[set]] <- these.train
  
  set.seed(2000 + set)
  nobs <- 0
  nsamp <- length(these.train)  # make sure that we get the same sample size
  while (nobs < 3) {
    these.train <- sort(sample(length(Y1), nsamp))
    y.o <- Y1[these.train]
    nobs <- sum(y.o)
  }
  srs.lst.Y1.250[[set]] <- these.train
}

for (set in 1:nsets) {
  set.seed(1000 + set)
  nobs <- 0
  while(nobs < 3) {  
    # keep repeating the sampling until there are at least 3 observations
    these.train <- sort(sample(length(Y2), n))
    y.o <- Y2[these.train]
    these.cluster <- y.o == 1
    these.cluster.ids <- these.train[these.cluster]
    for (i in 1:length(these.cluster.ids)) {
      # get the location of the cell where y == 1
      this.cell <- get.arr.idx(these.cluster.ids[i], nrows = 200)
      this.row  <- this.cell[1]
      this.col  <- this.cell[2]
      
      # account for the boundary
      neighbors.row <- c(this.row + 1, this.row - 1)
      neighbors.col <- c(this.col + 1, this.col - 1)
      neighbors.row <- neighbors.row[neighbors.row > 0 & neighbors.row < 200]
      neighbors.col <- neighbors.col[neighbors.col > 0 & neighbors.col < 200]
      for (j in 1:length(neighbors.row)) {
        these.train <- c(these.train,
                         get.idx(row = neighbors.row[j], col = this.col, 
                                 nrows = 200))
      }
      for (j in 1:length(neighbors.col)) {
        these.train <- c(these.train, 
                         get.idx(row = this.row, col = neighbors.col[j], 
                                 nrows = 200))
      }
      
    }
    these.train <- sort(unique(these.train))
    y.o <- Y2[these.train]
    nobs <- sum(y.o)
  }
  clu.lst.Y2.250[[set]] <- these.train
  
  set.seed(2000 + set)
  nobs <- 0
  nsamp <- length(these.train)  # make sure that we get the same sample size
  while (nobs < 3) {
    these.train <- sort(sample(length(Y2), nsamp))
    y.o <- Y2[these.train]
    nobs <- sum(y.o)
  }
  srs.lst.Y2.250[[set]] <- these.train
}

nsets <- 200
n     <- 500

clu.lst.Y1.500 <- clu.lst.Y2.500 <- vector(mode = "list", length = nsets)
srs.lst.Y1.500 <- srs.lst.Y2.500 <- vector(mode = "list", length = nsets)
for (set in 1:nsets) {
  set.seed(1000 + set)
  nobs <- 0
  while(nobs < 3) {  
    # keep repeating the sampling until there are at least 3 observations
    these.train <- sort(sample(length(Y1), n))
    y.o <- Y1[these.train]
    these.cluster <- y.o == 1
    these.cluster.ids <- these.train[these.cluster]
    for (i in 1:length(these.cluster.ids)) {
      # get the location of the cell where y == 1
      this.cell <- get.arr.idx(these.cluster.ids[i], nrows = 200)
      this.row  <- this.cell[1]
      this.col  <- this.cell[2]
      
      # account for the boundary
      neighbors.row <- c(this.row + 1, this.row - 1)
      neighbors.col <- c(this.col + 1, this.col - 1)
      neighbors.row <- neighbors.row[neighbors.row > 0 & neighbors.row < 200]
      neighbors.col <- neighbors.col[neighbors.col > 0 & neighbors.col < 200]
      for (j in 1:length(neighbors.row)) {
        these.train <- c(these.train,
                         get.idx(row = neighbors.row[j], col = this.col, 
                                 nrows = 200))
      }
      for (j in 1:length(neighbors.col)) {
        these.train <- c(these.train, 
                         get.idx(row = this.row, col = neighbors.col[j], 
                                 nrows = 200))
      }
      
    }
    these.train <- sort(unique(these.train))
    y.o <- Y1[these.train]
    nobs <- sum(y.o)
  }
  clu.lst.Y1.500[[set]] <- these.train
  
  set.seed(2000 + set)
  nobs <- 0
  nsamp <- length(these.train)  # make sure that we get the same sample size
  while (nobs < 3) {
    these.train <- sort(sample(length(Y1), nsamp))
    y.o <- Y1[these.train]
    nobs <- sum(y.o)
  }
  srs.lst.Y1.500[[set]] <- these.train
}

for (set in 1:nsets) {
  set.seed(1000 + set)
  nobs <- 0
  while(nobs < 3) {  
    # keep repeating the sampling until there are at least 3 observations
    these.train <- sort(sample(length(Y2), n))
    y.o <- Y2[these.train]
    these.cluster <- y.o == 1
    these.cluster.ids <- these.train[these.cluster]
    for (i in 1:length(these.cluster.ids)) {
      # get the location of the cell where y == 1
      this.cell <- get.arr.idx(these.cluster.ids[i], nrows = 200)
      this.row  <- this.cell[1]
      this.col  <- this.cell[2]
      
      # account for the boundary
      neighbors.row <- c(this.row + 1, this.row - 1)
      neighbors.col <- c(this.col + 1, this.col - 1)
      neighbors.row <- neighbors.row[neighbors.row > 0 & neighbors.row < 200]
      neighbors.col <- neighbors.col[neighbors.col > 0 & neighbors.col < 200]
      for (j in 1:length(neighbors.row)) {
        these.train <- c(these.train,
                         get.idx(row = neighbors.row[j], col = this.col, 
                                 nrows = 200))
      }
      for (j in 1:length(neighbors.col)) {
        these.train <- c(these.train, 
                         get.idx(row = this.row, col = neighbors.col[j], 
                                 nrows = 200))
      }
      
    }
    these.train <- sort(unique(these.train))
    y.o <- Y2[these.train]
    nobs <- sum(y.o)
  }
  clu.lst.Y2.500[[set]] <- these.train
  
  set.seed(2000 + set)
  nobs <- 0
  nsamp <- length(these.train)  # make sure that we get the same sample size
  while (nobs < 3) {
    these.train <- sort(sample(length(Y2), nsamp))
    y.o <- Y2[these.train]
    nobs <- sum(y.o)
  }
  srs.lst.Y2.500[[set]] <- these.train
}

save(Y1, Y2, s, 
     clu.lst.Y1.100, clu.lst.Y2.100, srs.lst.Y1.100, srs.lst.Y2.100, 
     clu.lst.Y1.250, clu.lst.Y2.250, srs.lst.Y1.250, srs.lst.Y2.250,
     clu.lst.Y1.500, clu.lst.Y2.500, srs.lst.Y1.500, srs.lst.Y2.500, 
     file = "plant_inventory.RData")