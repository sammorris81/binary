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
method.types <- c("gev", "pro", "log")

nmethods <- length(method.types)
nsets <- 100
nfits <- nsets * length(ns) * length(samp.types)

bs.results <- auc.results <- vector(length = 2,
                                    mode = "list")
bs.results.1 <- bs.results.0 <- vector(length = 2,
                                       mode = "list")

for (i in 1:2) {
  bs.results[[i]]   <- matrix(NA, nfits, nmethods)
  bs.results.1[[i]] <- matrix(NA, nfits, nmethods)
  bs.results.0[[i]] <- matrix(NA, nfits, nmethods)
  auc.results[[i]]  <- matrix(NA, nfits, nmethods)
  these.rownames <- c(paste("clu-100-", 1:nsets, sep = ""),
                      paste("srs-100-", 1:nsets, sep = ""),
                      paste("clu-250-", 1:nsets, sep = ""),
                      paste("srs-250-", 1:nsets, sep = ""))
  these.colnames <- method.types
  rownames(bs.results[[i]])   <- these.rownames
  rownames(auc.results[[i]])  <- these.rownames
  colnames(bs.results[[i]])   <- these.colnames
  colnames(auc.results[[i]])  <- these.colnames
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
  colnames(bs.results.combined[[species.idx]])  <- method.types
  colnames(auc.results.combined[[species.idx]]) <- method.types

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
  auc.results.combined[[species.idx]][4, ] <- this.row
}

prefix <- "ss-results/"
files <- list.files(path = prefix)

bs.results.1 <- bs.results.0 <- vector(length = 2, mode = "list")
for (i in 1:2) {
  bs.results.1[[i]] <- matrix(NA, nfits, nmethods)
  bs.results.0[[i]] <- matrix(NA, nfits, nmethods)
  these.rownames <- c(paste("clu-100-", 1:nsets, sep = ""),
                      paste("srs-100-", 1:nsets, sep = ""),
                      paste("clu-250-", 1:nsets, sep = ""),
                      paste("srs-250-", 1:nsets, sep = ""))
  these.colnames <- c("gev", "probit", "logit")
  rownames(bs.results.1[[i]]) <- these.rownames
  rownames(bs.results.0[[i]]) <- these.rownames
  colnames(bs.results.1[[i]]) <- these.colnames
  colnames(bs.results.0[[i]]) <- these.colnames
}

for (i in 1:length(files)) {
  split         <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  species.idx   <- as.numeric(split[2])
  samp.type.idx <- which(samp.types == split[1])
  n.idx         <- which(ns == as.numeric(split[3]))
  set           <- as.numeric(split[4])
  this.row      <- (n.idx - 1) * nsets * 2 + (samp.type.idx - 1) * nsets + set
  results.file  <- paste("ss-results/", files[i], sep = "")
  
  load(results.file)
  bs.results.1[[species.idx]][this.row, 1] <- mean((y.p[y.p == 1] - post.prob.gev[y.p == 1])^2)
  bs.results.1[[species.idx]][this.row, 2] <- mean((y.p[y.p == 1] - post.prob.pro[y.p == 1])^2)
  bs.results.1[[species.idx]][this.row, 3] <- mean((y.p[y.p == 1] - post.prob.log[y.p == 1])^2)
  
  bs.results.0[[species.idx]][this.row, 1] <- mean((y.p[y.p == 0] - post.prob.gev[y.p == 0])^2)
  bs.results.0[[species.idx]][this.row, 2] <- mean((y.p[y.p == 0] - post.prob.pro[y.p == 0])^2)
  bs.results.0[[species.idx]][this.row, 3] <- mean((y.p[y.p == 0] - post.prob.log[y.p == 0])^2)
}


bs.results.1.combined <- vector(mode = "list", length = 2)
bs.results.0.combined <- vector(mode = "list", length = 2)
for (i in 1:2) {
  species.idx <- i
  these.nrows <- length(ns) * length(samp.types)
  bs.results.1.combined[[species.idx]] <- matrix(NA, these.nrows, nmethods)
  bs.results.0.combined[[species.idx]] <- matrix(NA, these.nrows, nmethods)
  these.rownames <- paste(rep(samp.types, 2), "-", rep(ns, each = 2), sep = "")
  rownames(bs.results.1.combined[[species.idx]]) <- these.rownames
  rownames(bs.results.0.combined[[species.idx]]) <- these.rownames
  colnames(bs.results.1.combined[[species.idx]]) <- c("gev", "probit", "logit")
  colnames(bs.results.0.combined[[species.idx]]) <- c("gev", "probit", "logit")
  
  this.row <- apply(bs.results.1[[species.idx]][1:nsets, ],
                    2, mean, na.rm = TRUE)
  bs.results.1.combined[[species.idx]][1, ] <- this.row
  this.row <- apply(bs.results.1[[species.idx]][(nsets + 1):(2 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.1.combined[[species.idx]][2, ] <- this.row
  this.row <- apply(bs.results.1[[species.idx]][(2 * nsets +1):(3 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.1.combined[[species.idx]][3, ] <- this.row
  this.row <- apply(bs.results.1[[species.idx]][(3 * nsets + 1):(4 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.1.combined[[species.idx]][4, ] <- this.row
  
  this.row <- apply(bs.results.0[[species.idx]][1:nsets, ],
                    2, mean, na.rm = TRUE)
  bs.results.0.combined[[species.idx]][1, ] <- this.row
  this.row <- apply(bs.results.0[[species.idx]][(nsets + 1):(2 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.0.combined[[species.idx]][2, ] <- this.row
  this.row <- apply(bs.results.0[[species.idx]][(2 * nsets + 1):(3 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.0.combined[[species.idx]][3, ] <- this.row
  this.row <- apply(bs.results.0[[species.idx]][(3 * nsets + 1):(4 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.0.combined[[species.idx]][4, ] <- this.row
}

# Apparently ROCR can take in a list over all 100 datasets to 
# come up with an averaged cross-validation curve

for (species.idx in 1:2) { for (sample.idx in 1:2) { for (n.idx in 1:2) {
  this.samp    <- samp.types[sample.idx]
  this.species <- species.idx
  this.n       <- ns[n.idx]
  
  this.gev.pred <- vector(mode = "list", length = nsets)
  this.pro.pred <- vector(mode = "list", length = nsets)
  this.log.pred <- vector(mode = "list", length = nsets)
  this.yp       <- vector(mode = "list", length = nsets)
  sets.done <- rep(FALSE, nsets)
  
  for (set in 1:nsets) {
    results.file <- paste("./ss-results/", this.samp, "-", this.species, "-",
                          this.n, "-", set, ".RData", sep = "")
    if (file.exists(results.file)) {
      sets.done[set] <- TRUE
      load(results.file)
      this.yp[[set]] <- y.p
      this.gev.pred[[set]] <- post.prob.gev
      this.pro.pred[[set]] <- post.prob.pro
      this.log.pred[[set]] <- post.prob.log
    }
  }
  
  pred.gev.name <- paste("pred.gev.", this.samp, ".",
                         this.species, ".", this.n, sep = "")
  pred.pro.name <- paste("pred.pro.", this.samp, ".",
                         this.species, ".", this.n, sep = "")
  pred.log.name <- paste("pred.log.", this.samp, ".",
                         this.species, ".", this.n, sep = "")
  yp.name   <- paste("yp.", this.samp, ".", this.species, ".", this.n, sep = "")
  # the syntax here is single bracket to access multiple elements of the list
  assign(pred.gev.name, this.gev.pred[sets.done])
  assign(pred.pro.name, this.pro.pred[sets.done])
  assign(pred.log.name, this.log.pred[sets.done])
  assign(yp.name, this.yp[sets.done])
  print(paste(this.samp, this.species, this.n))
}}}

#### look at over ROC curves and PRC curves
pred.gev <- prediction(pred.gev.clu.1.100, yp.clu.1.100)
pred.pro <- prediction(pred.pro.clu.1.100, yp.clu.1.100)
pred.log <- prediction(pred.log.clu.1.100, yp.clu.1.100)
quartz(width = 16, height = 8)
main <- "Species 1, Cluster sample, n = 100"
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main)
dev.print(device = pdf, "./plots/perf-clu-1-100-thresh.pdf")
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main, avg = "vertical")
dev.print(device = pdf, "./plots/perf-clu-1-100-vert.pdf")
dev.off()

pred.gev <- prediction(pred.gev.clu.2.100, yp.clu.2.100)
pred.pro <- prediction(pred.pro.clu.2.100, yp.clu.2.100)
pred.log <- prediction(pred.log.clu.2.100, yp.clu.2.100)
quartz(width = 16, height = 8)
main <- "Species 2, Cluster sample, n = 100"
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main)
dev.print(device = pdf, "./plots/perf-clu-2-100-thresh.pdf")
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main, avg = "vertical")
dev.print(device = pdf, "./plots/perf-clu-2-100-vert.pdf")
dev.off()

pred.gev <- prediction(pred.gev.srs.1.100, yp.srs.1.100)
pred.pro <- prediction(pred.pro.srs.1.100, yp.srs.1.100)
pred.log <- prediction(pred.log.srs.1.100, yp.srs.1.100)
quartz(width = 16, height = 8)
main <- "Species 1, Simple random sample, n = 100"
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main)
dev.print(device = pdf, "./plots/perf-srs-1-100-thresh.pdf")
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main, avg = "vertical")
dev.print(device = pdf, "./plots/perf-srs-1-100-vert.pdf")
dev.off()

pred.gev <- prediction(pred.gev.srs.2.100, yp.srs.2.100)
pred.pro <- prediction(pred.pro.srs.2.100, yp.srs.2.100)
pred.log <- prediction(pred.log.srs.2.100, yp.srs.2.100)
quartz(width = 16, height = 8)
main <- "Species 2, Simple random sample, n = 100"
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main)
dev.print(device = pdf, "./plots/perf-srs-2-100-thresh.pdf")
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main, avg = "vertical")
dev.print(device = pdf, "./plots/perf-srs-2-100-vert.pdf")
dev.off()

pred.gev <- prediction(pred.gev.clu.1.250, yp.clu.1.250)
pred.pro <- prediction(pred.pro.clu.1.250, yp.clu.1.250)
pred.log <- prediction(pred.log.clu.1.250, yp.clu.1.250)
quartz(width = 16, height = 8)
main <- "Species 1, Cluster sample, n = 250"
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main)
dev.print(device = pdf, "./plots/perf-clu-1-250-thresh.pdf")
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main, avg = "vertical")
dev.print(device = pdf, "./plots/perf-clu-1-250-vert.pdf")
dev.off()

pred.gev <- prediction(pred.gev.clu.2.250, yp.clu.2.250)
pred.pro <- prediction(pred.pro.clu.2.250, yp.clu.2.250)
pred.log <- prediction(pred.log.clu.2.250, yp.clu.2.250)
quartz(width = 16, height = 8)
main <- "Species 2, Cluster sample, n = 250"
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main)
dev.print(device = pdf, "./plots/perf-clu-2-250-thresh.pdf")
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main, avg = "vertical")
dev.print(device = pdf, "./plots/perf-clu-2-250-vert.pdf")
dev.off()

pred.gev <- prediction(pred.gev.srs.1.250, yp.srs.1.250)
pred.pro <- prediction(pred.pro.srs.1.250, yp.srs.1.250)
pred.log <- prediction(pred.log.srs.1.250, yp.srs.1.250)
quartz(width = 16, height = 8)
main <- "Species 1, Simple random sample, n = 250"
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main)
dev.print(device = pdf, "./plots/perf-srs-1-250-thresh.pdf")
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main, avg = "vertical")
dev.print(device = pdf, "./plots/perf-srs-1-250-vert.pdf")
dev.off()

pred.gev <- prediction(pred.gev.srs.2.250, yp.srs.2.250)
pred.pro <- prediction(pred.pro.srs.2.250, yp.srs.2.250)
pred.log <- prediction(pred.log.srs.2.250, yp.srs.2.250)
quartz(width = 16, height = 8)
main <- "Species 2, Simple random sample, n = 250"
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main)
dev.print(device = pdf, "./plots/perf-srs-2-250-thresh.pdf")
plot.roc.prc(pred.gev, pred.pro, pred.log, main = main, avg = "vertical")
dev.print(device = pdf, "./plots/perf-srs-2-250-vert.pdf")
dev.off()


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
                      ", n: ", n, ", Set: ", set.idx, sep = ""))
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
  
  title.roc <- paste("ROC Curve: ", title.sample, " of Species ", species.idx, 
                     " with ", title.n, sep = "")
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
  
  title.prc <- paste("Precision Recall Curve: ", title.sample, " of Species ", 
                     species.idx, " with ", title.n, sep = "")
  p <- ggplot(df.prc, aes(x = x, y = y, color = method))
  p <- p + geom_smooth()
  p <- p + labs(title = title.prc, x = "Recall",
                y = "Precision", color = "Method")
  p <- p + theme_bw()
  ggsave(prc.plot.file, plot = p, width = 8, height = 8)
}}}