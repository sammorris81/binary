library(ggplot2)
library(colorspace)
library(scales)
plot.post.heatmap <- function(df, main, zlim, midpoint = NULL, 
                              legend.title = NULL) {
  if (is.null(legend.title)) {
    legend.title <- as.expression(names(df)[1])
  }
  if (is.null(midpoint)) {
    midpoint <- min((zlim[2] - zlim[1]) / 2)
  }
  p <- ggplot(df, aes(x = s1, y = s2))
  p <- p + geom_raster(aes(fill = Y))
  # p <- p + scale_fill_gradient(legend.title, low = "white",
  #                              high = "firebrick4", limits = zlim)
  p <- p + scale_fill_gradient2(legend.title, low = "dodgerblue4",
                                mid = "white", high = "firebrick4",
                                limits = zlim, midpoint = midpoint, 
                                na.value = "grey7")
  p <- p + labs(x = NULL, y = NULL, title = main)
  p <- p + theme_bw()
  p <- p + theme(axis.ticks = element_blank(), 
                 axis.text = element_blank(), 
                 panel.grid.major = element_blank(), 
                 panel.background = element_blank(), 
                 panel.border = element_blank())
  return(p)
}

plot.species <- function(df, main, legend.title = NULL) {
  if (is.null(legend.title)) {
    legend.title <- as.expression(names(df)[1])
  }
  p <- ggplot(df, aes(x = s1, y = s2))
  p <- p + geom_raster(aes(fill = Y))
  p <- p + scale_fill_manual(legend.title, 
                             values = c("dodgerblue4", "firebrick4"), 
                             labels = c("Not obs", "Obs"))
  p <- p + labs(x = NULL, y = NULL, title = main)
  p <- p + theme_bw()
  p <- p + theme(axis.ticks = element_blank(), 
                 axis.text = element_blank(), 
                 panel.grid.major = element_blank(), 
                 panel.background = element_blank(), 
                 panel.border = element_blank())
  return(p)
}