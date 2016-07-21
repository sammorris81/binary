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
  p <- p + scale_fill_gradient2(legend.title, low = "skyblue1",
                                mid = "plum3", high = "firebrick3",
                                limits = zlim, midpoint = midpoint, 
                                na.value = "grey7")
  p <- p + labs(x = NULL, y = NULL, title = main)
  p <- p + coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))
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
                             values = c("skyblue1", "firebrick3"), 
                             labels = c("Not observed", "Observed"))
  p <- p + labs(x = NULL, y = NULL, title = main)
  p <- p + coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))
  p <- p + theme_bw()
  p <- p + theme(axis.ticks = element_blank(), 
                 axis.text = element_blank(), 
                 panel.grid.major = element_blank(), 
                 panel.background = element_blank(), 
                 panel.border = element_blank(),
                 legend.title = element_text(size = 20, face = "bold.italic", 
                                             vjust = 1),
                 legend.text = element_text(size = 16),
                 legend.key.height = unit(24, 'pt'),
                 legend.key.width  = unit(24, 'pt'))
  return(p)
}

plot.roc.prc <- function(pred.gev, pred.pro, pred.log, avg = "vertical", 
                         main = NULL) {
  roc.gev  <- performance(pred.gev, "tpr", "fpr")
  roc.pro  <- performance(pred.pro, "tpr", "fpr")
  roc.log  <- performance(pred.log, "tpr", "fpr")
  prc.gev  <- performance(pred.gev, "prec", "rec")
  prc.pro  <- performance(pred.pro, "prec", "rec")
  prc.log  <- performance(pred.log, "prec", "rec")
  
  roc.main <- paste("ROC Curve: ", main, sep = "")
  prc.main <- paste("Precision Recall Curve: ", main, sep = "")
  
  par(mfrow = c(1, 2))
  plot(roc.gev, avg = avg, col = "grey20", lwd = 2, main = roc.main)
  plot(roc.pro, avg = avg, col = "firebrick2", lwd = 2, add = TRUE)
  plot(roc.log, avg = avg, col = "dodgerblue2", lwd = 2, add = TRUE)
  legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"), 
         legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)

  plot(prc.gev, avg = avg, col = "grey20", lwd = 2, main = prc.main)
  plot(prc.pro, avg = avg, col = "firebrick2", lwd = 2, add = TRUE)
  plot(prc.log, avg = avg, col = "dodgerblue2", lwd = 2, add = TRUE)
  # legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"), 
  #        legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)
}