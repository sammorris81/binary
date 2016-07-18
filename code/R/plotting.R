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

plot.roc.prc <- function(pred.gev, pred.pro, pred.log) {
  roc.gev  <- performance(pred.gev, "tpr", "fpr")
  roc.pro  <- performance(pred.pro, "tpr", "fpr")
  roc.log  <- performance(pred.log, "tpr", "fpr")
  prc.gev  <- performance(pred.gev, "prec", "rec")
  prc.pro  <- performance(pred.pro, "prec", "rec")
  prc.log  <- performance(pred.log, "prec", "rec")
  
  par(mfrow = c(1, 2))
  plot(roc.gev@x.values[[1]], roc.gev@y.values[[1]], 
       xlim = c(0, 1), ylim = c(0, 1),
       type = "l", col = "grey20", 
       main = "ROC Curve", lwd = 2, 
       ylab = "True positive rate", xlab = "False positive rate")
  lines(roc.pro@x.values[[1]], roc.pro@y.values[[1]],
        col = "firebrick2", lwd = 2)
  lines(roc.log@x.values[[1]], roc.log@y.values[[1]],
        col = "dodgerblue2", lwd = 2)
  legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"), 
         legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)

  plot(prc.gev@x.values[[1]], prc.gev@y.values[[1]],
       xlim = c(0, 1), ylim = c(0, 1),
       type = "l", col = "grey20", 
       main = "Precision Recall Curve", lwd = 2, 
       ylab = "Precision", xlab = "Recall")
  lines(roc.pro@x.values[[1]], roc.pro@y.values[[1]],
        col = "firebrick2", lwd = 2)
  lines(roc.log@x.values[[1]], roc.log@y.values[[1]],
        col = "dodgerblue2", lwd = 2)
  # legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"), 
  #        legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)
}