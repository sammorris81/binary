rm(list = ls())

library(ggplot2)
library(gridExtra)
library(pROC)
library(ROCR)
source("./package_load.R")
load("plant_inventory.RData")
source("../../R/plotting.R", chdir = TRUE)
source("../../../../usefulR/usefulfunctions.R", chdir = TRUE)
prefix <- "ss-tables/"
files <- list.files(path = prefix)
ns <- c(100, 250)
samp.types <- c("clu", "srs")

nmethods <- 3
nsets <- 100
nfits <- nsets * length(ns) * length(samp.types)

bs.results <- auc.results <- vector(length = 2,
                                    mode = "list")
for (i in 1:2) {
  bs.results[[i]]  <- matrix(NA, nfits, nmethods)
  auc.results[[i]] <- matrix(NA, nfits, nmethods)
  these.rownames <- c(paste("clu-100-", 1:nsets, sep = ""),
                      paste("srs-100-", 1:nsets, sep = ""),
                      paste("clu-250-", 1:nsets, sep = ""),
                      paste("srs-250-", 1:nsets, sep = ""))
  these.colnames <- c("gev", "probit", "logit")
  rownames(bs.results[[i]])  <- these.rownames
  rownames(auc.results[[i]]) <- these.rownames
  colnames(bs.results[[i]])  <- these.colnames
  colnames(auc.results[[i]]) <- these.colnames
}

for (i in 1:length(files)) {
  split         <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  species.idx   <- as.numeric(split[2])
  samp.type.idx <- which(samp.types == split[1])
  n.idx         <- which(ns == as.numeric(split[3]))
  set           <- as.numeric(split[4])
  this.row      <- (n.idx - 1) * nsets * 2 + (samp.type.idx - 1) * nsets + set
  table.set     <- read.table(paste(prefix, files[i], sep = ""))
  if (all(!is.na(table.set))) {
    bs.results[[species.idx]][this.row, 1]  <- table.set[1, 1]
    auc.results[[species.idx]][this.row, 1] <- table.set[1, 2]
    bs.results[[species.idx]][this.row, 2]  <- table.set[2, 1]
    auc.results[[species.idx]][this.row, 2] <- table.set[2, 2]
    bs.results[[species.idx]][this.row, 3]  <- table.set[3, 1]
    auc.results[[species.idx]][this.row, 3] <- table.set[3, 2]
  }
}

bs.results.combined  <- vector(mode = "list", length = 2)
auc.results.combined <- vector(mode = "list", length = 2)
species.rate <- c(mean(Y1), mean(Y2))
for (i in 1:2) {
  species.idx <- i
  these.nrows <- length(ns) * length(samp.types)
  bs.results.combined[[species.idx]] <- matrix(NA, these.nrows, nmethods)
  auc.results.combined[[species.idx]] <- matrix(NA, these.nrows, nmethods)
  these.rownames <- paste(rep(samp.types, 2), "-", rep(ns, each = 2), sep = "")
  rownames(bs.results.combined[[species.idx]])  <- these.rownames
  rownames(auc.results.combined[[species.idx]]) <- these.rownames
  colnames(bs.results.combined[[species.idx]])  <- c("gev", "probit", "logit")
  colnames(auc.results.combined[[species.idx]]) <- c("gev", "probit", "logit")

  this.row <- apply(bs.results[[species.idx]][1:nsets, ],
                    2, mean, na.rm = TRUE)
  bs.results.combined[[species.idx]][1, ] <- this.row
  this.row <- apply(bs.results[[species.idx]][(nsets + 1):(2 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.combined[[species.idx]][2, ] <- this.row
  this.row <- apply(bs.results[[species.idx]][(2 * nsets +1):(3 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.combined[[species.idx]][3, ] <- this.row
  this.row <- apply(bs.results[[species.idx]][(3 * nsets + 1):(4 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.combined[[species.idx]][4, ] <- this.row

  this.row <- apply(auc.results[[species.idx]][1:nsets, ],
                    2, mean, na.rm = TRUE)
  auc.results.combined[[species.idx]][1, ] <- this.row
  this.row <- apply(auc.results[[species.idx]][(nsets + 1):(2 * nsets), ],
                    2, mean, na.rm = TRUE)
  auc.results.combined[[species.idx]][2, ] <- this.row
  this.row <- apply(auc.results[[species.idx]][(2 * nsets + 1):(3 * nsets), ],
                    2, mean, na.rm = TRUE)
  auc.results.combined[[species.idx]][3, ] <- this.row
  this.row <- apply(auc.results[[species.idx]][(3 * nsets + 1):(4 * nsets), ],
                    2, mean, na.rm = TRUE)
  auc.results.combined[[species.idx]][3, ] <- this.row
}

#### look at over ROC curves and PRC curves
#### TODO: Add in BS for Y = 1 and Y = 0
for (species.idx in 1:2) { for (sample.idx in 1:2) { for (n.idx in 1:2) {
  for (set.idx in 1:nsets) {
    which.y <- paste("Y", species.idx, sep = "")
    samp.type <- samp.types[sample.idx]
    n <- ns[n.idx]
    results.file <- paste("./ss-results/", samp.type, "-", species.idx, "-", n, 
                          "-", set.idx, ".RData", sep = "")
    ss.list <- paste(samp.type, ".lst.", which.y, sep = "")
    post.plot.file <- paste("./plots/post-prob-", samp.type, "-", species.idx, 
                            "-", n, "-", set.idx, ".pdf", sep = "")
    roc.plot.file <- paste("./plots/roc-", samp.type, "-", species.idx, 
                           "-", n, ".pdf", sep = "")
    prc.plot.file <- paste("./plots/prc-", samp.type, "-", species.idx, 
                           "-", n, ".pdf", sep = "")
    
    if (file.exists(results.file)) {
      load(results.file)
      
      # Looping over all the datasets to get all the points for the ROC and PRC
      # curves. Then running a smoother to come up with a single curve for the 
      # method, sampling type, and sample size
      pred.gev <- prediction(post.prob.gev, y.p)
      pred.pro <- prediction(post.prob.pro, y.p)
      pred.log <- prediction(post.prob.log, y.p)
      
      roc.gev  <- performance(pred.gev, "tpr", "fpr")
      roc.pro  <- performance(pred.pro, "tpr", "fpr")
      roc.log  <- performance(pred.log, "tpr", "fpr")
      prc.gev  <- performance(pred.gev, "prec", "rec")
      prc.pro  <- performance(pred.pro, "prec", "rec")
      prc.log  <- performance(pred.log, "prec", "rec")
      
      if (set.idx == 1) {
        roc.x <- roc.gev@x.values[[1]]
        roc.y <- roc.gev@y.values[[1]]
        prc.x <- prc.gev@x.values[[1]]
        prc.y <- prc.gev@y.values[[1]]
        
        roc.x <- c(roc.x, roc.pro@x.values[[1]], roc.log@x.values[[1]])
        roc.y <- c(roc.y, roc.pro@y.values[[1]], roc.log@y.values[[1]])
        prc.x <- c(prc.x, prc.pro@x.values[[1]], prc.log@x.values[[1]])
        prc.y <- c(prc.y, prc.pro@y.values[[1]], prc.log@y.values[[1]])
        
        roc.method <- c(rep("GEV", length(roc.gev@x.values[[1]])),
                        rep("Probit", length(roc.pro@x.values[[1]])),
                        rep("Logit", length(roc.log@x.values[[1]])))
        prc.method <- c(rep("GEV", length(prc.gev@x.values[[1]])),
                        rep("Probit", length(prc.pro@x.values[[1]])),
                        rep("Logit", length(prc.log@x.values[[1]])))
      } else {
        roc.x <- c(roc.x, roc.gev@x.values[[1]], 
                   roc.pro@x.values[[1]], roc.log@x.values[[1]])
        roc.y <- c(roc.y, roc.gev@y.values[[1]], 
                   roc.pro@y.values[[1]], roc.log@y.values[[1]])
        prc.x <- c(prc.x, prc.gev@x.values[[1]], 
                   prc.pro@x.values[[1]], prc.log@x.values[[1]])
        prc.y <- c(prc.y, prc.gev@y.values[[1]], 
                   prc.pro@y.values[[1]], prc.log@y.values[[1]])
        
        roc.method <- c(roc.method, 
                        rep("GEV", length(roc.gev@x.values[[1]])),
                        rep("Probit", length(roc.pro@x.values[[1]])),
                        rep("Logit", length(roc.log@x.values[[1]])))
        prc.method <- c(prc.method,
                        rep("GEV", length(prc.gev@x.values[[1]])),
                        rep("Probit", length(prc.pro@x.values[[1]])),
                        rep("Logit", length(prc.log@x.values[[1]])))
      
        if (set.idx %% 10 == 0) {
          obs <- as.vector(as.factor(get(which.y)))
          df <- data.frame(Y = obs, s1 = s[, 1], s2 = s[, 2])
          main <- "Census of species 1"
          legend.title <- "Species 1"
          p1 <- plot.species(df = df, main = main, legend.title = legend.title)
          
          zlim <- range(c(post.prob.gev, post.prob.pro, post.prob.log))
          post.s       <- rbind(s.p, s.o)
          legend.title <- "P(Y = 1)"
          
          post.prob <- c(post.prob.gev, rep(NA, length = nrow(s.o)))
          df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
          main <- "GEV"
          p2 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                                  legend.title = legend.title)
          
          post.prob <- c(post.prob.pro, rep(NA, length = nrow(s.o)))
          df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
          main <- "Probit"
          p3 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                                  legend.title = legend.title)
          
          post.prob <- c(post.prob.log, rep(NA, length = nrow(s.o)))
          df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
          main <- "Logit"
          p4 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                                  legend.title = legend.title)
          
          layout.mtx <- matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE)
          panel <- arrangeGrob(p1, p2, p3, p4, ncol = 2, layout_matrix = layout.mtx)
          ggsave(post.plot.file, plot = panel, width = 13, height = 8)
          
          print(paste("Species: ", species.idx, ", Sampling: ", samp.type, 
                      ", n: ", n, ", Set", set.idx, sep = ""))
        }
        
      }
    }
  }
  
  df.roc <- data.frame(x = roc.x, y = roc.y, method = roc.method) 
  df.roc$method <- factor(df.roc$method, levels = c("GEV", "Probit", "Logit"))
  
  if (sample.idx == 1) {
    title.sample <- "Cluster"
  } else {
    title.sample <- "Simple Random Sample"
  }
  title.n <- paste("n = ", n, sep = "")
  
  title.roc <- paste("ROC Curve: ", title.sample, " with ", title.n, sep = "")
  p <- ggplot(df.roc, aes(x = x, y = y, color = method))
  p <- p + geom_smooth()
  p <- p + labs(title = title.roc, x = "False positive rate", 
                y = "True positive rate", color = "Method")
  p <- p + theme_bw()
  ggsave(roc.plot.file, plot = p, width = 8, height = 8)
  
  prc.method <- prc.method[!is.nan(prc.y)]
  prc.x      <- prc.x[!is.nan(prc.y)]
  prc.y      <- prc.y[!is.nan(prc.y)]
  df.prc <- data.frame(x = prc.x, y = prc.y, method = prc.method)
  df.prc$method <- factor(df.prc$method, levels = c("GEV", "Probit", "Logit"))
  
  title.prc <- paste("Precision Recall Curve: ", title.sample, " with ", 
                     title.n, sep = "")
  p <- ggplot(df.prc, aes(x = x, y = y, color = method))
  p <- p + geom_smooth()
  p <- p + labs(title = title.prc, x = "Recall",
                y = "Precision", color = "Method")
  p <- p + theme_bw()
  ggsave(prc.plot.file, plot = p, width = 8, height = 8)
}}}




# ROC Curve


plot.roc.prc(pred.gev, pred.pro, pred.log)


roc <- performance(pred.gev, "tpr", "fpr")
x.plot <- seq(0, 1, 0.001)
bin.n <- rep(0, length(x.plot))
bin.sum <- rep(0, length(x.plot))
for (i in 1:length(roc@x.values[[1]])) {
  this.bin <- which.min(roc@x.values[[1]][i] > x.plot)
  bin.n[this.bin] <- bin.n[this.bin] + 1
  bin.sum[this.bin] <-bin.sum[this.bin] + roc@y.values[[1]][i]
}


pred.gev.2 <- prediction(post.prob.gev, y.p)
roc.2 <- performance(pred.gev.2, "tpr", "fpr")

roc.gev <- roc(y.p ~ post.prob.gev)
roc.pro <- roc(y.p ~ post.prob.pro)
roc.log <- roc(y.p ~ post.prob.log)
quartz(width = 8, height = 8)

dev.print(device = pdf, file = paste("plots/roc-clu-1-", this.set, ".pdf",
                                     sep = ""))
dev.off()


#### MSE for predicting Y
# lower is better for each of these
mean((y.p[y.p == 1] - post.prob.gev[y.p == 1])^2)  # 0.5764
mean((y.p[y.p == 1] - post.prob.pro[y.p == 1])^2)  # 0.6281
mean((y.p[y.p == 1] - post.prob.log[y.p == 1])^2)  # 0.9010

mean((y.p[y.p == 0] - post.prob.gev[y.p == 0])^2)  # 0.0113
mean((y.p[y.p == 0] - post.prob.pro[y.p == 0])^2)  # 0.0084
mean((y.p[y.p == 0] - post.prob.log[y.p == 0])^2)  # 0.0003

#### look at posterior probability species 2
this.set <- 7
load(paste("./ss-results/clu-2-100-", this.set, ".RData", sep = ""))
plot.prob.gev <- Y2
plot.prob.gev[-clu.lst.Y2[[this.set]]] <- post.prob.gev
plot.prob.gev[Y2 == 1] <- 1

plot.prob.pro <- Y2
plot.prob.pro[-clu.lst.Y2[[this.set]]] <- post.prob.pro
plot.prob.pro[Y2 == 1] <- 1

obs <- as.vector(as.factor(Y2))
df <- data.frame(Y = obs, s1 = s[, 1], s2 = s[, 2])
main <- "Actual species map"
legend.title <- "Species 2"
p1 <- plot.species(df = df, main = main, legend.title = legend.title)

zlim         <- range(c(post.prob.gev, post.prob.pro, post.prob.log))
post.s       <- rbind(s.p, s.o)
legend.title <- "P(Y = 1)"

post.prob <- c(post.prob.gev, rep(NA, length = nrow(s.o)))
df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
main <- "GEV"
p2 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                        legend.title = legend.title)

post.prob <- c(post.prob.pro, rep(NA, length = nrow(s.o)))
df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
main <- "Probit"
p3 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                        legend.title = legend.title)

post.prob <- c(post.prob.log, rep(NA, length = nrow(s.o)))
df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
main <- "Logit"
p4 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                        legend.title = legend.title)

layout.mtx <- matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE)
panel <- arrangeGrob(p1, p2, p3, p4, ncol = 2, layout_matrix = layout.mtx)
ggsave(paste("plots/post-prob-2-", this.set, ".pdf", sep = ""), plot = panel)

# ROC Curve
roc.gev <- roc(y.p ~ post.prob.gev)
roc.pro <- roc(y.p ~ post.prob.pro)
roc.log <- roc(y.p ~ post.prob.log)
quartz(width = 8, height = 8)
plot(roc.gev, col = "grey20", main = "Cluster sample: Species 2")
plot(roc.pro, add = TRUE, col = "firebrick2")
plot(roc.log, add = TRUE, col = "dodgerblue2")
legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"),
       legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)
dev.print(device = pdf, file = paste("plots/roc-clu-2-", this.set, ".pdf",
                                     sep = ""))
dev.off()

#### MSE for predicting Y
# lower is better for each of these
mean((y.p[y.p == 1] - post.prob.gev[y.p == 1])^2)  # 0.9120
mean((y.p[y.p == 1] - post.prob.pro[y.p == 1])^2)  # 0.9222
mean((y.p[y.p == 1] - post.prob.log[y.p == 1])^2)  # 0.9712

mean((y.p[y.p == 0] - post.prob.gev[y.p == 0])^2)  # 0.0010
mean((y.p[y.p == 0] - post.prob.pro[y.p == 0])^2)  # 0.0013
mean((y.p[y.p == 0] - post.prob.log[y.p == 0])^2)  # 0.0002

#### look at posterior probability species 1
this.set <- 2
load(paste("./ss-results/srs-1-100-", this.set, ".RData", sep = ""))
plot.prob.gev <- Y1
plot.prob.gev[-clu.lst.Y1[[this.set]]] <- post.prob.gev
plot.prob.gev[Y1 == 1] <- 1

plot.prob.pro <- Y1
plot.prob.pro[-clu.lst.Y1[[this.set]]] <- post.prob.pro
plot.prob.pro[Y1 == 1] <- 1

obs <- as.vector(as.factor(Y1))
df <- data.frame(Y = obs, s1 = s[, 1], s2 = s[, 2])
main <- "Actual species map"
legend.title <- "Species 1"
p1 <- plot.species(df = df, main = main, legend.title = legend.title)

zlim <- range(c(post.prob.gev, post.prob.pro, post.prob.log))
post.s       <- rbind(s.p, s.o)
legend.title <- "P(Y = 1)"

post.prob <- c(post.prob.gev, rep(NA, length = nrow(s.o)))
df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
main <- "GEV"
p2 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                        legend.title = legend.title)

post.prob <- c(post.prob.pro, rep(NA, length = nrow(s.o)))
df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
main <- "Probit"
p3 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                        legend.title = legend.title)

post.prob <- c(post.prob.log, rep(NA, length = nrow(s.o)))
df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
main <- "Logit"
p4 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                        legend.title = legend.title)

layout.mtx <- matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE)
panel <- arrangeGrob(p1, p2, p3, p4, ncol = 2, layout_matrix = layout.mtx)
ggsave(paste("plots/post-prob-srs-1-", this.set, ".pdf", sep = ""),
       plot = panel)

# ROC Curve
roc.gev <- roc(y.p ~ post.prob.gev)
roc.pro <- roc(y.p ~ post.prob.pro)
roc.log <- roc(y.p ~ post.prob.log)
quartz(width = 8, height = 8)
plot(roc.gev, col = "grey20", main = "Simple random sample: Species 1")
plot(roc.pro, add = TRUE, col = "firebrick2")
plot(roc.log, add = TRUE, col = "dodgerblue2")
legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"),
       legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)
dev.print(device = pdf, file = paste("plots/roc-srs-1-", this.set, ".pdf",
                                     sep = ""))
dev.off()

#### MSE for predicting Y
# lower is better for each of these
mean((y.p[y.p == 1] - post.prob.gev[y.p == 1])^2)  # 0.0586
mean((y.p[y.p == 1] - post.prob.pro[y.p == 1])^2)  # 0.9032
mean((y.p[y.p == 1] - post.prob.log[y.p == 1])^2)  # 0.9787

mean((y.p[y.p == 0] - post.prob.gev[y.p == 0])^2)  # 0.0258
mean((y.p[y.p == 0] - post.prob.pro[y.p == 0])^2)  # 0.0013
mean((y.p[y.p == 0] - post.prob.log[y.p == 0])^2)  # 0.0001

#### look at posterior probability species 2
this.set <- 2
load(paste("./ss-results/srs-2-100-", this.set, ".RData", sep = ""))
plot.prob.gev <- Y2
plot.prob.gev[-srs.lst.Y2[[this.set]]] <- post.prob.gev
plot.prob.gev[Y1 == 1] <- 1

plot.prob.pro <- Y2
plot.prob.pro[-srs.lst.Y2[[this.set]]] <- post.prob.pro
plot.prob.pro[Y2 == 1] <- 1

obs <- as.vector(as.factor(Y2))
df <- data.frame(Y = obs, s1 = s[, 1], s2 = s[, 2])
main <- "Actual species map"
legend.title <- "Species 2"
p1 <- plot.species(df = df, main = main, legend.title = legend.title)

zlim <- range(c(post.prob.gev, post.prob.pro, post.prob.log))
post.s       <- rbind(s.p, s.o)
legend.title <- "P(Y = 1)"

post.prob <- c(post.prob.gev, rep(NA, length = nrow(s.o)))
df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
main <- "GEV"
p2 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                        legend.title = legend.title)

post.prob <- c(post.prob.pro, rep(NA, length = nrow(s.o)))
df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
main <- "Probit"
p3 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                        legend.title = legend.title)

post.prob <- c(post.prob.log, rep(NA, length = nrow(s.o)))
df <- data.frame(Y = post.prob, s1 = post.s[, 1], s2 = post.s[, 2])
main <- "Logit"
p4 <- plot.post.heatmap(df = df, main = main, zlim = zlim,
                        legend.title = legend.title)

layout.mtx <- matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE)
panel <- arrangeGrob(p1, p2, p3, p4, ncol = 2, layout_matrix = layout.mtx)
ggsave(paste("plots/post-prob-srs-2-", this.set, ".pdf", sep = ""),
       plot = panel)

# ROC Curve
roc.gev <- roc(y.p ~ post.prob.gev)
roc.pro <- roc(y.p ~ post.prob.pro)
roc.log <- roc(y.p ~ post.prob.log)
quartz(width = 8, height = 8)
plot(roc.gev, col = "grey20", main = "Simple random sample: Species 2")
plot(roc.pro, add = TRUE, col = "firebrick2")
plot(roc.log, add = TRUE, col = "dodgerblue2")
legend("bottomright", col = c("grey20", "firebrick2", "dodgerblue2"),
       legend = c("GEV", "Probit", "Logit"), lty = 1, lwd = 2)
dev.print(device = pdf, file = paste("plots/roc-srs-2-", this.set, ".pdf",
                                     sep = ""))
dev.off()

#### MSE for predicting Y
# lower is better for each of these
mean((y.p[y.p == 1] - post.prob.gev[y.p == 1])^2)  # 0.0586
mean((y.p[y.p == 1] - post.prob.pro[y.p == 1])^2)  # 0.9032
mean((y.p[y.p == 1] - post.prob.log[y.p == 1])^2)  # 0.9787

mean((y.p[y.p == 0] - post.prob.gev[y.p == 0])^2)  # 0.0258
mean((y.p[y.p == 0] - post.prob.pro[y.p == 0])^2)  # 0.0013
mean((y.p[y.p == 0] - post.prob.log[y.p == 0])^2)  # 0.0001
