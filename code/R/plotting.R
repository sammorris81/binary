library(ggplot2)
library(gridExtra)
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
                 legend.title = element_text(size = rel(1.5), 
                                             face = "bold.italic",
                                             vjust = 1),
                 legend.text = element_text(size = rel(1.5)),
                 # legend.key.height = unit(24, 'pt'),
                 # legend.key.width  = unit(24, 'pt'),
                 plot.title = element_text(size = rel(1.5)))
  return(p)
}

plot.simdata <- function(y, s) {
  plot(s[y == 1, 1], s[y == 1, 2], col = "firebrick4", bg = "firebrick1",
       pch = 21, xlim = c(0, 1), ylim = c(0, 1),
       main = paste(round(100 * mean(y), 2), "%"))
  points(s[y == 0, 1], s[y == 0, 2], col = "dodgerblue4", bg = "dodgerblue1",
         pch = 21)
}

plot.roc.gg <- function(pred.gev.100, pred.pro.100, pred.log.100,
                     pred.gev.250, pred.pro.250, pred.log.250,
                     main = NULL, span = 1) {
  roc.gev.clu <- performance(pred.gev.clu, "tpr", "fpr")
  roc.pro.clu <- performance(pred.pro.clu, "tpr", "fpr")
  roc.log.clu <- performance(pred.log.clu, "tpr", "fpr")
  roc.gev.srs <- performance(pred.gev.srs, "tpr", "fpr")
  roc.pro.srs <- performance(pred.pro.srs, "tpr", "fpr")
  roc.log.srs <- performance(pred.log.srs, "tpr", "fpr")

  roc.main <- paste("ROC Curve: ", main, sep = "")
  # prc.main <- paste("Precision Recall Curve: ", main, sep = "")
  
  ngev.clu <- length(unlist(roc.gev.clu@x.values))
  npro.clu <- length(unlist(roc.pro.clu@x.values))
  nlog.clu <- length(unlist(roc.log.clu@x.values))
  ngev.srs <- length(unlist(roc.gev.srs@x.values))
  npro.srs <- length(unlist(roc.pro.srs@x.values))
  nlog.srs <- length(unlist(roc.log.srs@x.values))
  x.vals <- c(unlist(roc.gev.clu@x.values), unlist(roc.pro.clu@x.values), 
              unlist(roc.log.clu@x.values), unlist(roc.gev.srs@x.values),
              unlist(roc.pro.srs@x.values), unlist(roc.log.srs@x.values))
  y.vals <- c(unlist(roc.gev.clu@y.values), unlist(roc.pro.clu@y.values), 
              unlist(roc.log.clu@y.values), unlist(roc.gev.srs@y.values),
              unlist(roc.pro.srs@y.values), unlist(roc.log.srs@y.values))
  method <- as.factor(c(rep("GEV", ngev.clu), rep("Probit", npro.clu), 
                        rep("Logistic", nlog.clu), rep("GEV", ngev.srs),
                        rep("Probit", npro.srs), rep("Logistic", nlog.srs)))
  sample <- as.factor(c(rep("clu", ngev.clu + npro.clu + nlog.clu), 
                        rep("srs", ngev.srs + npro.srs + nlog.srs)))
  df <- data.frame(x.vals, y.vals, method, sample)
  
  p1 <- ggplot(df, aes(x = x.vals, y = y.vals, color = method, linetype = sample))
  p1 <- p1 + geom_smooth()
  p1 <- p1 + labs(x = "False positive rate", y = "Average true positive rate",
                  title = roc.main)
  return(p1)
  
  # plot(roc.gev, avg = avg, col = "grey20", lwd = 2, main = roc.main)
  # plot(roc.pro, avg = avg, col = "firebrick2", lwd = 2, add = TRUE)
  # plot(roc.log, avg = avg, col = "dodgerblue2", lwd = 2, add = TRUE)
  # legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"),
  #        legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)

  # plot(prc.gev, avg = avg, col = "grey20", lwd = 2, main = prc.main,
  #      ylim = c(0, 1))
  # plot(prc.pro, avg = avg, col = "firebrick2", lwd = 2, add = TRUE)
  # plot(prc.log, avg = avg, col = "dodgerblue2", lwd = 2, add = TRUE)
  # legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"),
  #        legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)
}

plot.roc <- function(pred.gev, pred.pro, pred.log,
                     avg = "vertical", main = NULL) {
  roc.gev <- performance(pred.gev, "tpr", "fpr")
  roc.pro <- performance(pred.pro, "tpr", "fpr")
  roc.log <- performance(pred.log, "tpr", "fpr")
  
  roc.main <- paste("ROC Curve: ", main, sep = "")
  # prc.main <- paste("Precision Recall Curve: ", main, sep = "")
 
  # plot(roc.gev, avg = avg, col = "grey20", main = roc.main)
  plot(roc.gev, avg = avg, col = "grey20", spread.estimate = "stderror", 
       show.spread.at = seq(0.1, 0.9, 0.1), main = roc.main, 
       cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  # plot(roc.pro, avg = avg, col = "firebrick2", add = TRUE)
  plot(roc.pro, avg = avg, col = "firebrick2", plotCI.col = "firebrick2", 
       spread.estimate = "stderror", show.spread.at = seq(0.1, 0.9, 0.1), 
       add = TRUE)
  # plot(roc.log, avg = avg, col = "dodgerblue2", add = TRUE)
  plot(roc.log, avg = avg, plotCI.col = "dodgerblue2", col = "dodgerblue2", 
       spread.estimate = "stderror", show.spread.at = seq(0.1, 0.9, 0.1), 
       add = TRUE)

  legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"),
         title = "Link function:",
         legend = c("GEV", "Probit", "Logistic"), 
         lty = 1, lwd = 1.5, cex = 1.5)

  # plot(prc.gev, avg = avg, col = "grey20", lwd = 2, main = prc.main,
  #      ylim = c(0, 1))
  # plot(prc.pro, avg = avg, col = "firebrick2", lwd = 2, add = TRUE)
  # plot(prc.log, avg = avg, col = "dodgerblue2", lwd = 2, add = TRUE)
  # legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"),
  #        legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)
}

plot.prc <- function(pred.gev, pred.pro, pred.log,
                     avg = "vertical", main = NULL) {
  roc.gev <- performance(pred.gev, "prec", "rec")
  roc.pro <- performance(pred.pro, "prec", "rec")
  roc.log <- performance(pred.log, "prec", "rec")
  
  roc.main <- paste("PRC Curve: ", main, sep = "")
  # prc.main <- paste("Precision Recall Curve: ", main, sep = "")
  
  # plot(roc.gev, avg = avg, col = "grey20", main = roc.main)
  plot(roc.gev, avg = avg, col = "grey20", spread.estimate = "stderror", 
       show.spread.at = seq(0.1, 0.9, 0.1), main = roc.main, 
       cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  # plot(roc.pro, avg = avg, col = "firebrick2", add = TRUE)
  plot(roc.pro, avg = avg, col = "firebrick2", plotCI.col = "firebrick2", 
       spread.estimate = "stderror", show.spread.at = seq(0.1, 0.9, 0.1), 
       add = TRUE)
  # plot(roc.log, avg = avg, col = "dodgerblue2", add = TRUE)
  plot(roc.log, avg = avg, plotCI.col = "dodgerblue2", col = "dodgerblue2", 
       spread.estimate = "stderror", show.spread.at = seq(0.1, 0.9, 0.1), 
       add = TRUE)
  
  legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"),
         title = "Link function:",
         legend = c("GEV", "Probit", "Logistic"), 
         lty = 1, lwd = 1.5, cex = 1.5)
  
  # plot(prc.gev, avg = avg, col = "grey20", lwd = 2, main = prc.main,
  #      ylim = c(0, 1))
  # plot(prc.pro, avg = avg, col = "firebrick2", lwd = 2, add = TRUE)
  # plot(prc.log, avg = avg, col = "dodgerblue2", lwd = 2, add = TRUE)
  # legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"),
  #        legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)
}

plot.smooth <- function(df) {
  p1 <- ggplot(df, aes(x = rareness, y = bs, color = method))
  p1 <- p1 + geom_point(alpha = 0.8)
  p1 <- p1 + geom_smooth(alpha = 0.4)
  # p1 <- p1 + scale_x_continuous(breaks = seq(0.05, 0.4, by = 0.02))
  p1 <- p1 + scale_color_manual(values = c("grey20", "firebrick2", "dodgerblue2"))
  p1 <- p1 + labs(x = "Rareness", y = "Brier Score", 
                  title = "Rareness vs Brier Score")
  p1 <- p1 + theme_bw()
  
  p2 <- ggplot(df, aes(x = rareness, y = auc, color = method))
  p2 <- p2 + geom_point(alpha = 0.8)
  p2 <- p2 + geom_smooth(alpha = 0.4)
  # p2 <- p2 + scale_x_continuous(breaks = seq(0.05, 0.4, by = 0.02))
  p2 <- p2 + scale_color_manual(values = c("grey20", "firebrick2", "dodgerblue2"))
  p2 <- p2 + labs(x = "Rareness", y = "AUROC", 
                  title = "Rareness vs AUROC")
  p2 <- p2 + theme_bw()
  
  
  layout.mtx = matrix(1:2, nrow = 1, ncol = 2)
  panel <- arrangeGrob(p1, p2, ncol = 2, layout_matrix = layout.mtx)
  return(panel)
}