rm(list = ls())

library(ggplot2)
library(maps)
library(Rcpp)
source("../../code/R/aux_gev.R", chdir = TRUE)
source("../../../usefulR/usefulfunctions.R")

# Include functions to generate latitude/longitude pairs to cover the US
# These functions were written by Neal Grantham
generate_grid <- function(x = 100, y = 100){
  # purpose: generate grid of points over the US at which to make predictions
  # returns: two column matrix, longitude and latitude of all grid points
  library(maps)  # for map.where
  corners <- enclose_USA()
  # create grid
  grid <- expand.grid(seq(corners[1], corners[2], length = x), 
                      seq(corners[3], corners[4], length = y))
  # retain only points that fall over US soil
  inUSA <- !is.na(map.where("usa", x = grid[, 1], y = grid[, 2]))
  grid <- as.matrix(grid[inUSA, ])
  # name the columns and rows
  colnames(grid) <- c("lon", "lat")
  rownames(grid) <- 1:nrow(grid)
  return(grid)
}

enclose_USA <- function() {
  # purpose: geographic coordinate limits within the continental US
  # returns: vector of corners of USA
  max.lat <- 49.384472    # most northern
  min.lat <- 24.520833    # most southern
  max.lon <- -66.947028   # most eastern
  min.lon <- -124.733056  # most western
  corners <- c(min.lon, max.lon, min.lat, max.lat)
}

# Mapping functions
map.points <- function (lat, lon, data, bwFlag=F,
                        xlim=NULL, ylim=NULL, zlim=NULL,
                        mainTitle=NULL, legendTitle=NULL) { 
  # map.points creates a point map of a data vector given a set of locations
  #
  # Args:
  #   lat - vector (n x 1) of latitude values
  #   lon - vector (n x 1) of longitude values
  #   data - vector (n x 1) of data values to be plotted
  #   bwFlag - boolean indicator to force the color scheme to black/white
  #   xlim - vector (2 x 1) of longitude values to serve as (min, max)
  #   ylim - vector (2 x 1) of latitude values to serve as (min, max)
  #   zlim - vector (2 x 1) of data values to serve as cutoff values (min, max)
  #   mainTitle - character string to be used as the map title
  #   legendTitle - character string to be used as the legend heading
  #
  # Returns:
  #   ggplot2 object
  
  # Store the base data of the underlying map
  baseData <- map_data("state")
  
  # Combine the data into a dataframe
  dfMap <- as.data.frame(cbind(lon, lat, data))
  colnames(dfMap) <- c("lon", "lat", "Value")
  
  # Set limits for x, y, z if not specified as parameters
  if (is.null(xlim)) { xlim <- c(min(baseData$long)-1, max(baseData$long)+1) }
  if (is.null(ylim)) { ylim <- c(min(baseData$lat)-1, max(baseData$lat)+1) }
  if (is.null(zlim)) { zlim <- range(data) }
  
  # Remove observations that fall outside of the range set in zlim
  rmWhich <- (dfMap$Value < min(zlim)) + (dfMap$Value > max(zlim))
  dfMapRm <- dfMap[rmWhich==0, ]
  
  # Create the plot
  p <- ggplot(dfMapRm, aes(x=lon, y=lat)) + theme_bw()
  p <- p + theme(panel.border=element_blank(), panel.grid=element_blank(),
                 axis.ticks=element_blank(), axis.text=element_blank())
  p <- p + geom_point(aes(colour = Value)) + scale_colour_gradient(limits=zlim)
  p <- p + geom_polygon(data=baseData, aes(x=long, y=lat, group=group), 
                        colour="black", fill="white", alpha=0) 
  p <- p + labs(title=mainTitle, x="", y="")
  p <- p + coord_fixed(ratio=1.1, xlim=xlim, ylim=ylim)
  
  # Change the plot to a black/white scale if necessary
  if (bwFlag) { p <- p + scale_colour_gradient(low="white", high="black",
                                               limits=zlim) }
  
  # Change the legend title if necessary
  if (!is.null(legendTitle)) { p <- p + labs(colour=legendTitle) }
  
  return(p)  
}

map.heatmap <- function (lat, lon, data, bwFlag=F,
                         xlim=NULL, ylim=NULL, zlim=NULL,
                         mainTitle=NULL, legendTitle=NULL) {
  # map.heatmap creates a heatmap of a data vector given a set of locations
  #
  # Args:
  #   lat - vector (n x 1) of latitude values
  #   lon - vector (n x 1) of longitude values
  #   data - vector (n x 1) of data values to be plotted
  #   bwFlag - boolean indicator to force the color scheme to black/white
  #   xlim - vector (2 x 1) of longitude values to serve as (min, max)
  #   ylim - vector (2 x 1) of latitude values to serve as (min, max)
  #   zlim - vector (2 x 1) of data values to serve as cutoff values (min, max)
  #   mainTitle - character string to be used as the map title
  #   legendTitle - character string to be used as the legend heading
  #
  # Returns:
  #   ggplot2 object
  
  # Store the base data of the underlying map
  baseData <- map_data("state")
  
  # Combine the data into a dataframe
  dfMap <- as.data.frame(cbind(lon, lat, data))
  colnames(dfMap) <- c("lon", "lat", "Value")
  
  # Set limits for x, y, z if not specified as parameters
  if (is.null(xlim)) { xlim <- c(min(baseData$long)-1, max(baseData$long)+1) }
  if (is.null(ylim)) { ylim <- c(min(baseData$lat)-1, max(baseData$lat)+1) }
  if (is.null(zlim)) { zlim <- range(data) }
  
  # Remove observations that fall outside of the range set in zlim
  rmWhich <- (dfMap$Value < min(zlim)) + (dfMap$Value > max(zlim))
  dfMapRm <- dfMap[rmWhich==0, ]
  
  # Create the plot
  p <- ggplot(dfMapRm, aes(x=lon, y=lat, fill=Value)) + theme_bw()
  p <- p + theme(panel.border=element_blank(), panel.grid=element_blank(),
                 axis.ticks=element_blank(), axis.text=element_blank())
  p <- p + geom_tile()
  p <- p + geom_polygon(data=baseData, aes(x=long, y=lat, group=group), 
                        colour="black", fill="white", alpha=0) 
  p <- p + labs(title=mainTitle, x="", y="")
  p <- p + coord_fixed(ratio=1.1, xlim=xlim, ylim=ylim)
  
  # Change the plot to a black/white scale if necessary
  if (bwFlag) { p <- p + scale_fill_gradient(low="white", high="black",
                                             limits=zlim) 
  } else {
    p <- p + scale_fill_gradient(limits=zlim)
  }
  
  # Change the legend title if necessary
  if (!is.null(legendTitle)) { p <- p + labs(colour=legendTitle) }
  
  return(p)  
}

rRareBinarySpat <- function(x, s, knots, beta, xi, alpha, rho, nt = 1,
                            prob.success = 0.05, dw2 = NULL, a = NULL) {
  
  p <- length(beta)
  if (nt == 1) {
    ns <- nrow(x)
  } else {
    if (nrow(x) %% nt != 0) {
      stop("The number of rows of x must be a multiple of the number of days")
    }
    ns <- nrow(x) / nt
  }
  
  y      <- matrix(NA, ns, nt)
  nknots <- nrow(knots)
  
  x.beta <- matrix(0, ns, nt)
  
  # get weights
  if (is.null(dw2)) {  # for predictions, already have dw2
    dw2 <- as.matrix(rdist(s, knots))^2  # dw2 is ns x nknots
  }
  
  w <- getW(rho = rho, dw2 = dw2, a.cutoff = NULL) # w is ns x nknots
  
  # get random effects and theta.star
  if (is.null(a)) {
    a <- matrix(rPS(n = nknots * nt, alpha = alpha), nknots, nt)
  }
  
  w.star <- getWStar(alpha = alpha, w = w)
  theta  <- (w.star %*% a)^alpha
  
  # get underlying latent variable
  u <- matrix(rgev(n = ns * nt, 1, alpha, alpha), ns, nt)
  z <- u * theta
  
  # z = u * theta ~ GEV(1, 1, 1)
  # h = x.beta + (z^xi - 1) / xi ~ GEV(x.beta, 1, xi)
  if (xi == 0) {
    h <- x.beta + log(z)
  } else {
    h <- x.beta + (z^xi - 1) / xi
  }
  
  # set the threshold for success at whatever value will give us
  # our desired percentage of 1s.
  thresh <- quantile(h, probs = (1 - prob.success))
  
  h <- h - thresh
  prob <- 1 - exp(-exp(h))
  
  return(list(h = h, prob = prob))
}

library(fields)
library(SpatialTools)
set.seed(100)
s <- as.matrix(expand.grid(x = seq(0, 1, length = 75), 
                           y = seq(0, 1, length = 75)))
d <- rdist(s)
diag(d) <- 0
nknots <- 441
ns <- nrow(s)
rho <- 0.05
var <- 7

knots.grid <- as.matrix(expand.grid(x = seq(0, 1, length = 21), 
                                    y = seq(0, 1, length = 21)))
these.knots <- sample(x = ns, size = nknots, replace = FALSE)
knots.rand <- s[these.knots, ]

# generate a MVN surface using evenly spaced knots GP


# generate a MVN surface using randomly selected knots GP


# genearte a MVN surface using all the sites as knots
cov <- simple.cov.sp(D=d, sp.type="exponential", sp.par=c(var, rho), 
                     error.var=0, finescale.var=0)
z.log.3 <- -4 + t(chol(cov)) %*% rnorm(ns)
p.log.3 <- transform$inv.logit(z.log.3)


# generate a mGEV surface using evenly spaced knots
gev.1 <- rRareBinarySpat(x = matrix(1, ns, 1), s = s, knots = knots.grid, 
                           beta = 0, xi = 0, alpha = 0.3, rho = 0.05, nt = 1, 
                           prob.success = 0.05)
h.gev.1 <- gev.1$h
p.gev.1 <- gev.1$prob

# generate a mGEV surface using randomly selected knots
gev.2 <- rRareBinarySpat(x = matrix(1, ns, 1), s = s, knots = knots.rand, 
                           beta = 0, xi = 0, alpha = 0.3, rho = 0.05, nt = 1, 
                           prob.success = 0.05)
h.gev.2 <- gev.2$h
p.gev.2 <- gev.2$prob

# generate a mGEV surface using all the sites as knots
gev.3 <- rRareBinarySpat(x = matrix(1, ns, 1), s = s, knots = s, beta = 0, 
                           xi = 0, alpha = 0.3, rho = 0.05, nt = 1, 
                           prob.success = 0.05)
h.gev.3 <- gev.3$h
p.gev.3 <- gev.3$prob

# generate a mGEV surface using evenly spaced knots
knots.coarse <- as.matrix(expand.grid(x = seq(0, 1, length = 6),
                                      y = seq(0, 1, length = 6)))
gev.4 <- rRareBinarySpat(x = matrix(1, ns, 1), s = s, knots = knots.coarse, 
                         beta = 0, xi = 0, alpha = 0.3, rho = 0.08, nt = 1, 
                         prob.success = 0.05)
h.gev.4 <- gev.4$h
p.gev.4 <- gev.4$prob

par(mfrow = c(2, 2))
# quilt.plot(x = s[, 1], y = s[, 2], z = p.log.3, nx = 75, ny = 75, main = "logit")
quilt.plot(x = s[, 1], y = s[, 2], z = h.gev.1, nx = 75, ny = 75, main = "gev-1")
# points(knots.grid, pch = 21, col = "firebrick4", bg = "firebrick1")
quilt.plot(x = s[, 1], y = s[, 2], z = h.gev.2, nx = 75, ny = 75, main = "gev-2")
# points(knots.rand, pch = 21, col = "firebrick4", bg = "firebrick1")
quilt.plot(x = s[, 1], y = s[, 2], z = h.gev.3, nx = 75, ny = 75, main = "gev-3")
quilt.plot(x = s[, 1], y = s[, 2], z = h.gev.4, nx = 75, ny = 75, main = "gev-4")

# quilt.plot(x = s[, 1], y = s[, 2], z = p.log.3, nx = 75, ny = 75, main = "logit")
quilt.plot(x = s[, 1], y = s[, 2], z = p.gev.1, nx = 75, ny = 75, main = "gev-1")
# points(knots.grid, pch = 21, col = "firebrick4", bg = "firebrick1")
quilt.plot(x = s[, 1], y = s[, 2], z = p.gev.2, nx = 75, ny = 75, main = "gev-2")
# points(knots.rand, pch = 21, col = "firebrick4", bg = "firebrick1")
quilt.plot(x = s[, 1], y = s[, 2], z = p.gev.3, nx = 75, ny = 75, main = "gev-3")
quilt.plot(x = s[, 1], y = s[, 2], z = p.gev.4, nx = 75, ny = 75, main = "gev-4")