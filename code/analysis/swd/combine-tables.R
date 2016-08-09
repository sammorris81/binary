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
nsets <- 50
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
  if (set <= nsets) {
    this.row      <- (n.idx - 1) * nsets * 2 + (samp.type.idx - 1) * nsets + set
    table.set     <- read.table(paste(prefix, files[i], sep = ""))
    if (all(!is.na(table.set))) {
      bs.results[[species.idx]][this.row, 1]    <- table.set[1, 1]
      auc.results[[species.idx]][this.row, 1]   <- table.set[1, 2]
      bs.results.1[[species.idx]][this.row, 1]  <- table.set[1, 3]
      bs.results.0[[species.idx]][this.row, 1]  <- table.set[1, 4]
      bs.results[[species.idx]][this.row, 2]    <- table.set[2, 1]
      auc.results[[species.idx]][this.row, 2]   <- table.set[2, 2]
      bs.results.1[[species.idx]][this.row, 2]  <- table.set[2, 3]
      bs.results.0[[species.idx]][this.row, 2]  <- table.set[2, 4]
      bs.results[[species.idx]][this.row, 3]    <- table.set[3, 1]
      auc.results[[species.idx]][this.row, 3]   <- table.set[3, 2]
      bs.results.1[[species.idx]][this.row, 3]  <- table.set[3, 3]
      bs.results.0[[species.idx]][this.row, 3]  <- table.set[3, 4]
    }
  }
}

bs.results.combined  <- vector(mode = "list", length = 2)
bs.results.comb.se   <- vector(mode = "list", length = 2)
auc.results.combined <- vector(mode = "list", length = 2)
auc.results.comb.se  <- vector(mode = "list", length = 2)
bs.results.1.combined <- vector(mode = "list", length = 2)
bs.results.0.combined <- vector(mode = "list", length = 2)
species.rate <- c(mean(Y1), mean(Y2))
for (i in 1:2) {
  species.idx <- i
  these.nrows <- length(ns) * length(samp.types)
  bs.results.combined[[species.idx]]  <- matrix(NA, these.nrows, nmethods)
  bs.results.comb.se[[species.idx]]   <- matrix(NA, these.nrows, nmethods)
  auc.results.combined[[species.idx]] <- matrix(NA, these.nrows, nmethods)
  auc.results.comb.se[[species.idx]]  <- matrix(NA, these.nrows, nmethods)
  bs.results.1.combined[[species.idx]] <- matrix(NA, these.nrows, nmethods)
  bs.results.0.combined[[species.idx]] <- matrix(NA, these.nrows, nmethods)

  these.rownames <- paste(rep(samp.types, 2), "-", rep(ns, each = 2), sep = "")
  rownames(bs.results.combined[[species.idx]])   <- these.rownames
  rownames(bs.results.comb.se[[species.idx]])    <- these.rownames
  rownames(auc.results.combined[[species.idx]])  <- these.rownames
  rownames(auc.results.comb.se[[species.idx]])   <- these.rownames
  rownames(bs.results.1.combined[[species.idx]]) <- these.rownames
  rownames(bs.results.0.combined[[species.idx]]) <- these.rownames
  colnames(bs.results.combined[[species.idx]])   <- method.types
  colnames(auc.results.combined[[species.idx]])  <- method.types
  colnames(bs.results.1.combined[[species.idx]]) <- method.types
  colnames(bs.results.0.combined[[species.idx]]) <- method.types

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

  this.row <- apply(bs.results[[species.idx]][1:nsets, ],
                    2, sd, na.rm = TRUE) / sqrt(nsets)
  bs.results.comb.se[[species.idx]][1, ] <- this.row
  this.row <- apply(bs.results[[species.idx]][(nsets + 1):(2 * nsets), ],
                    2, sd, na.rm = TRUE) / sqrt(nsets)
  bs.results.comb.se[[species.idx]][2, ] <- this.row
  this.row <- apply(bs.results[[species.idx]][(2 * nsets +1):(3 * nsets), ],
                    2, sd, na.rm = TRUE) / sqrt(nsets)
  bs.results.comb.se[[species.idx]][3, ] <- this.row
  this.row <- apply(bs.results[[species.idx]][(3 * nsets + 1):(4 * nsets), ],
                    2, sd, na.rm = TRUE) / sqrt(nsets)
  bs.results.comb.se[[species.idx]][4, ] <- this.row

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

  this.row <- apply(auc.results[[species.idx]][1:nsets, ],
                    2, sd, na.rm = TRUE) / sqrt(nsets)
  auc.results.comb.se[[species.idx]][1, ] <- this.row
  this.row <- apply(auc.results[[species.idx]][(nsets + 1):(2 * nsets), ],
                    2, sd, na.rm = TRUE) / sqrt(nsets)
  auc.results.comb.se[[species.idx]][2, ] <- this.row
  this.row <- apply(auc.results[[species.idx]][(2 * nsets + 1):(3 * nsets), ],
                    2, sd, na.rm = TRUE) / sqrt(nsets)
  auc.results.comb.se[[species.idx]][3, ] <- this.row
  this.row <- apply(auc.results[[species.idx]][(3 * nsets + 1):(4 * nsets), ],
                    2, sd, na.rm = TRUE) / sqrt(nsets)
  auc.results.comb.se[[species.idx]][4, ] <- this.row

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
  this.row <- apply(bs.results.0[[species.idx]][(2 * nsets +1):(3 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.0.combined[[species.idx]][3, ] <- this.row
  this.row <- apply(bs.results.0[[species.idx]][(3 * nsets + 1):(4 * nsets), ],
                    2, mean, na.rm = TRUE)
  bs.results.0.combined[[species.idx]][4, ] <- this.row
}

prefix <- "ss-results/"
files <- list.files(path = prefix)

round(bs.results.combined[[1]] * 100, 3)
round(bs.results.comb.se[[1]] * 100, 3)
round(auc.results.combined[[1]], 3)
round(auc.results.comb.se[[1]], 3)

round(bs.results.combined[[2]] * 100, 3)
round(bs.results.comb.se[[2]] * 100, 3)
round(auc.results.combined[[2]], 3)
round(auc.results.comb.se[[2]], 3)

rareness <- matrix(0, nsets, 4)
y <- Y1
these.cols <- rep(NA, 4)

for (sample.idx in 1:2) {
  if (sample.idx == 1) {
    samp.type <- "clu"
  } else {
    samp.type <- "srs"
  }

  for (n.idx in 1:2) {
    if (n.idx == 1) {
      n <- 100
    } else {
      n <- 250
    }

    sample.list <- paste(samp.type, ".lst.Y1.", n, sep = "")
    this.col <- (sample.idx - 1) * 2 + n.idx
    these.cols[this.col] <- paste(samp.type, n, sep = "-")

    for (set in 1:nsets) {
      these.train <- get(sample.list)[[set]]
      y.o <- y[these.train]
      rareness[set, this.col] <- mean(y.o)
    }
  }
}

colnames(rareness) <- these.cols
#### smooth of rareness by brier score
for (setting in 1:4) {
  # only want a single sample type in each plot
  start <- (setting - 1) * nsets + 1
  end   <- setting * nsets

  plot.file <- paste("./plots/byrareness-", setting, ".pdf", sep = "")
  df <- data.frame(rareness = rep(rareness[, setting], 3),
                   bs = c(bs.results[[1]][start:end, 1],
                          bs.results[[1]][start:end, 2],
                          bs.results[[1]][start:end, 3]),
                   auc = c(auc.results[[1]][start:end, 1],
                           auc.results[[1]][start:end, 2],
                           auc.results[[1]][start:end, 3]),
                   method = as.factor(rep(c("GEV", "Probit", "Logistic"),
                                          each = length(start:end))))
  df$method <- factor(df$method, levels = c("GEV", "Probit", "Logistic"))

  panel <- plot.smooth(df)
  ggsave(plot.file, plot = panel, width = 9, height = 4.5)
}

# Apparently ROCR can take in a list over all 100 datasets to
# come up with an averaged cross-validation curve

for (species.idx in 1:2) {
  for (sample.idx in 1:2) {
    for (n.idx in 1:2) {
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

quartz(width = 16, height = 16)
par(mfrow = c(2, 2))
#### look at over ROC curves and PRC curves
pred.gev <- prediction(pred.gev.clu.1.100, yp.clu.1.100)
pred.pro <- prediction(pred.pro.clu.1.100, yp.clu.1.100)
pred.log <- prediction(pred.log.clu.1.100, yp.clu.1.100)
main <- expression(paste(italic("Tamarix ramosissima"), 
                         ": Cluster sample, n = 100", sep = ""))
plot.roc(pred.gev, pred.pro, pred.log, main = main)

pred.gev <- prediction(pred.gev.clu.1.250, yp.clu.1.250)
pred.pro <- prediction(pred.pro.clu.1.250, yp.clu.1.250)
pred.log <- prediction(pred.log.clu.1.250, yp.clu.1.250)
main <- expression(paste(italic("Tamarix ramosissima"), 
                         ": Cluster sample, n = 250", sep = ""))
plot.roc(pred.gev, pred.pro, pred.log, main = main)

pred.gev <- prediction(pred.gev.srs.1.100, yp.srs.1.100)
pred.pro <- prediction(pred.pro.srs.1.100, yp.srs.1.100)
pred.log <- prediction(pred.log.srs.1.100, yp.srs.1.100)
main <- expression(paste(italic("Tamarix ramosissima"), 
                         ": Simple random sample, n = 100", sep = ""))
plot.roc(pred.gev, pred.pro, pred.log, main = main)

pred.gev <- prediction(pred.gev.srs.1.250, yp.srs.1.250)
pred.pro <- prediction(pred.pro.srs.1.250, yp.srs.1.250)
pred.log <- prediction(pred.log.srs.1.250, yp.srs.1.250)
main <- expression(paste(italic("Tamarix ramosissima"), 
                         ": Simple random sample, n = 250", sep = ""))
plot.roc(pred.gev, pred.pro, pred.log, main = main)

dev.print(device = pdf, "./plots/data-perf-species1.pdf")
dev.off()

quartz(width = 16, height = 16)
par(mfrow = c(2, 2))
pred.gev <- prediction(pred.gev.clu.2.100, yp.clu.2.100)
pred.pro <- prediction(pred.pro.clu.2.100, yp.clu.2.100)
pred.log <- prediction(pred.log.clu.2.100, yp.clu.2.100)
main <- expression(paste(italic("Hedysarum scoparium"), 
                         ": Cluster sample, n = 100", sep = ""))
plot.roc(pred.gev, pred.pro, pred.log, main = main)

pred.gev <- prediction(pred.gev.clu.2.250, yp.clu.2.250)
pred.pro <- prediction(pred.pro.clu.2.250, yp.clu.2.250)
pred.log <- prediction(pred.log.clu.2.250, yp.clu.2.250)
main <- expression(paste(italic("Hedysarum scoparium"), 
                         ": Cluster sample, n = 250", sep = ""))
plot.roc(pred.gev, pred.pro, pred.log, main = main)

pred.gev <- prediction(pred.gev.srs.2.100, yp.srs.2.100)
pred.pro <- prediction(pred.pro.srs.2.100, yp.srs.2.100)
pred.log <- prediction(pred.log.srs.2.100, yp.srs.2.100)
main <- expression(paste(italic("Hedysarum scoparium"), 
                         ": Simple random sample, n = 100", sep = ""))
plot.roc(pred.gev, pred.pro, pred.log, main = main)

pred.gev <- prediction(pred.gev.srs.2.250, yp.srs.2.250)
pred.pro <- prediction(pred.pro.srs.2.250, yp.srs.2.250)
pred.log <- prediction(pred.log.srs.2.250, yp.srs.2.250)
main <- expression(paste(italic("Hedysarum scoparium"), 
                         ": Simple random sample, n = 250", sep = ""))
plot.roc(pred.gev, pred.pro, pred.log, main = main)

dev.print(device = pdf, "./plots/data-perf-species2.pdf")
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


    if (set.idx %% 10 == 0) {
      if (file.exists(results.file)) {
        load(results.file)

        obs <- as.vector(as.factor(get(which.y)))
        df <- data.frame(Y = obs, s1 = s[, 1], s2 = s[, 2])
        main <- paste("Census of species ", species.idx, sep = "")
        legend.title <- paste("Species ", species.idx, sep = "")
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
        ggsave(post.plot.file, plot = panel, width = 10, height = 8)

        print(paste("Species: ", species.idx, ", Sampling: ", samp.type,
                    ", n: ", n, ", Set: ", set.idx, sep = ""))

      }
    }
  }
}}}


df <- data.frame(Y = as.vector(as.factor(Y1)), s1 = s[, 1], s2 = s[, 2])
p1 <- plot.species(df = df, main = NULL, legend.title = "Tamarix\nramosissima")
ggsave("./plots/tamarix-census.pdf", plot = p1, width = 8, height = 8)

df <- data.frame(Y = as.vector(as.factor(Y2)), s1 = s[, 1], s2 = s[, 2])
p2 <- plot.species(df = df, main = NULL, legend.title = "Hedysarum\nscoparium")

layout.mtx <- matrix(1:2, nrow = 1, ncol = 2)
panel <- arrangeGrob(p1, p2, ncol = 2, layout_matrix = layout.mtx)
ggsave("./plots/plant-census.pdf", plot = panel, width = 16, height = 8)


timing.gev <- matrix(NA, nsets, 8)
timing.pro <- matrix(NA, nsets, 8)
timing.log <- matrix(NA, nsets, 8)
files <- list.files("./ss-fit")
for (i in 1:length(files)) {
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  species.idx <- as.numeric(split[2])
  if (split[1] == "clu") {
    samp.idx <- 1
  } else {
    samp.idx <- 2
  }
  if (as.numeric(split[3] == 100)) {
    n.idx <- 1
  } else {
    n.idx <- 2
  }
  set <- as.numeric(split[4])
  if (set <= 50) {
    setting <- (species.idx - 1) * 4 + (n.idx - 1) * 2 + samp.idx
    
    set       <- as.numeric(split[4])
    load(paste("./ss-fit/", files[i], sep = ""))
    timing.gev[set, setting] <- timings[1]
    timing.pro[set, setting] <- timings[2]
    timing.log[set, setting] <- timings[3]
  }
  if (i %% 10 == 0) {
    print(paste("Finished: ", i, " of ", length(files), sep = ""))
  }
}

round(apply(timing.gev * 60 / 25, 2, mean, na.rm = TRUE), 1)
apply(timing.pro, 2, mean, na.rm = TRUE)
apply(timing.log, 2, mean, na.rm = TRUE)
