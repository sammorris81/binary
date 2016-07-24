rm(list = ls())
source("../../R/plotting.R", chdir = TRUE)
load("simdata.RData")
tbl.dir <- "./sim-tables/"
res.dir <- "./sim-results/"
files <- list.files(path = tbl.dir)
library(ROCR)

nsettings <- 6
nsets <- 100
method.types <- c("gev", "pro", "log")
nmethods <- length(method.types)
rareness <- matrix(NA, nsets, nsettings)
for (setting in 1:nsettings) {
  for (set in 1:nsets) {
    rareness[set, setting] <- mean(simdata[[setting]]$y[, set])
  }
}

bs.results <- auc.results <- vector(length = nsettings, mode = "list")
bs.results.1 <- bs.results.0 <- vector(length = nsettings, mode = "list")
for (setting in 1:nsettings) {
  bs.results[[setting]]   <- matrix(NA, nsets, nmethods)
  bs.results.1[[setting]] <- matrix(NA, nsets, nmethods)
  bs.results.0[[setting]] <- matrix(NA, nsets, nmethods)
  auc.results[[setting]]  <- matrix(NA, nsets, nmethods)
  rownames(bs.results[[setting]])   <- paste("set", 1:nsets)
  rownames(bs.results.1[[setting]]) <- paste("set", 1:nsets)
  rownames(bs.results.1[[setting]]) <- paste("set", 1:nsets)
  rownames(auc.results[[setting]])  <- paste("set", 1:nsets)
  colnames(bs.results[[setting]])   <- method.types
  colnames(bs.results.1[[setting]]) <- method.types
  colnames(bs.results.0[[setting]]) <- method.types
  colnames(auc.results[[setting]])  <- method.types
}

these.sets <- which(rowSums(!is.na(auc.results[[1]])) == 3)
this.order <- these.sets[order(rareness[these.sets])]
plot(rareness[this.order, 1], bs.results[[1]][this.order, 1], type = "l")
lines(rareness[this.order, 1], bs.results[[1]][this.order, 2], lty = 2)
lines(rareness[this.order, 1], bs.results[[1]][this.order, 3], lty = 3)

for (i in 1:length(files)) {
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  setting   <- as.numeric(split[1])
  set       <- as.numeric(split[2])
  table.set <- read.table(paste(tbl.dir, files[i], sep = ""))
  bs.results[[setting]][set, 1]   <- table.set[1, 1]
  auc.results[[setting]][set, 1]  <- table.set[1, 2]
  bs.results.1[[setting]][set, 1] <- table.set[1, 3]
  bs.results.0[[setting]][set, 1] <- table.set[1, 4]
  bs.results[[setting]][set, 2]   <- table.set[2, 1]
  auc.results[[setting]][set, 2]  <- table.set[2, 2]
  bs.results.1[[setting]][set, 2] <- table.set[2, 3]
  bs.results.0[[setting]][set, 2] <- table.set[2, 4]
  bs.results[[setting]][set, 3]   <- table.set[3, 1]
  auc.results[[setting]][set, 3]  <- table.set[3, 2]
  bs.results.1[[setting]][set, 3] <- table.set[3, 3]
  bs.results.0[[setting]][set, 3] <- table.set[3, 4]
}

bs.results.combined   <- matrix(NA, nsettings, nmethods)
bs.results.1.combined <- matrix(NA, nsettings, nmethods)
bs.results.0.combined <- matrix(NA, nsettings, nmethods)
auc.results.combined  <- matrix(NA, nsettings, nmethods)
# vector of sets numbers that are complete
finished.sets <- vector(mode = "list", length = nsettings)
for (setting in 1:nsettings) {
  these.sets <- which(rowSums(is.na(bs.results[[setting]])) == 0)
  finished.sets[[setting]] <- these.sets
  bs.results.combined[setting, ]   <- apply(bs.results[[setting]][these.sets, ],
                                           2, mean, na.rm = TRUE)
  bs.results.1.combined[setting, ] <- apply(bs.results.1[[setting]][these.sets, ],
                                           2, mean, na.rm = TRUE)
  bs.results.0.combined[setting, ] <- apply(bs.results.0[[setting]][these.sets, ],
                                            2, mean, na.rm = TRUE)
  auc.results.combined[setting, ]  <- apply(auc.results[[setting]][these.sets, ],
                                           2, mean, na.rm = TRUE)
}

for (setting in 1:nsettings) {
  print(length(finished.sets[[setting]]))
}

# Apparently ROCR can take in a list over all 100 datasets to
# come up with an averaged cross-validation curve
files <- list.files(path = res.dir)
for (setting.idx in 1:6) {
  this.setting <- setting.idx

  this.gev.pred <- vector(mode = "list", length = nsets)
  this.pro.pred <- vector(mode = "list", length = nsets)
  this.log.pred <- vector(mode = "list", length = nsets)
  this.yp       <- vector(mode = "list", length = nsets)
  sets.done <- rep(FALSE, nsets)

  for (set in 1:nsets) {
    results.file <- paste(res.dir, this.setting, "-", set, ".RData", sep = "")
    if (file.exists(results.file)) {
      sets.done[set] <- TRUE
      load(results.file)
      this.yp[[set]] <- y.i.p
      this.gev.pred[[set]] <- post.prob.gev
      this.pro.pred[[set]] <- post.prob.pro
      this.log.pred[[set]] <- post.prob.log
    }
  }

  pred.gev.name <- paste("pred.gev.", this.setting, sep = "")
  pred.pro.name <- paste("pred.pro.", this.setting, sep = "")
  pred.log.name <- paste("pred.log.", this.setting, sep = "")
  yp.name       <- paste("yp.", this.setting, sep = "")

  # the syntax here is single bracket to access multiple elements of the list
  assign(pred.gev.name, this.gev.pred[sets.done])
  assign(pred.pro.name, this.pro.pred[sets.done])
  assign(pred.log.name, this.log.pred[sets.done])
  assign(yp.name, this.yp[sets.done])
  print(paste(this.setting))
}

#### look at over ROC curves and PRC curves
for (i in 1:nsettings) {
  this.gev <- paste("pred.gev.", i, sep = "")
  this.pro <- paste("pred.pro.", i, sep = "")
  this.log <- paste("pred.log.", i, sep = "")
  this.yp  <- paste("yp.", i, sep = "")

  pred.gev <- prediction(get(this.gev), get(this.yp))
  pred.pro <- prediction(get(this.pro), get(this.yp))
  pred.log <- prediction(get(this.log), get(this.yp))

  quartz(width = 16, height = 8)
  main <- paste("Setting ", i, sep = "")
  plot.roc.prc(pred.gev, pred.pro, pred.log, main = main)
  dev.print(device = pdf, paste("./plots/perf-setting-", i, ".pdf", sep = ""))
  dev.off()
}

#### smooth of rareness by brier score
for (setting in 1:nsettings) {
  these.include <- which(rowSums(!is.na(bs.results[[setting]])) == 3)
  plot.file <- paste("./plots/byrareness-", setting, ".pdf", sep = "")
  df <- data.frame(rareness = rep(rareness[these.include, 1], 3),
                   bs = c(bs.results[[setting]][these.include, 1], 
                          bs.results[[setting]][these.include, 2],
                          bs.results[[setting]][these.include, 3]),
                   auc = c(auc.results[[setting]][these.include, 1],
                           auc.results[[setting]][these.include, 2], 
                           auc.results[[setting]][these.include, 3]),
                   method = as.factor(rep(c("GEV", "Probit", "Logistic"), 
                                          each = length(these.include))))
  df$method <- factor(df$method, levels = c("GEV", "Probit", "Logistic"))
  
  panel <- plot.smooth(df)
  ggsave(plot.file, plot = panel, width = 8, height = 4)
}




# how many have finished
colSums(!is.na(bs.results[[1]]))
colSums(!is.na(bs.results[[2]]))
colSums(!is.na(bs.results[[3]]))
colSums(!is.na(bs.results[[4]]))
colSums(!is.na(bs.results[[5]]))
colSums(!is.na(bs.results[[6]]))

round(bs.results.combined[, 1] / bs.results.combined[, 3], 4)
round(bs.results.combined[, 2] / bs.results.combined[, 3], 4)
round(auc.results.combined[, 1], 4)
round(auc.results.combined[, 2], 4)
round(auc.results.combined[, 3], 4)

# Check for differences
# First do Friedman test (one-way repeated measures)
#   friedman.test(y ~ trt | block, data)
# Then follow up with the Wilcoxon, Nemenyi, McDonald-Thompson test
# pWNMT(x, b, trt, method, n.mc)
#     x: list of values
#     b: vector of blocks (only needed if x is a vector)
#     trt: vector of treatments
#     method: "Exact", "Monte Carlo" or "Asymptotic"

library(NSM3)
set.seed(6727)  #npar
groups <- as.factor(rep(1:nmethods, each=50))
dataset <- as.factor(rep(1:50, times=nmethods))
results.friedman <- matrix(0, nsettings, 2)
colnames(results.friedman) <- c("bs", "auc")

for (setting in 1:nsettings) {
  scores <- as.vector(bs.results[[setting]])
  combine <- data.frame(scores, groups, dataset)
  results.friedman[setting, 1] <- friedman.test(scores ~ groups | dataset,
                                                data=combine)$p.value

  scores <- as.vector(auc.results[[setting]])
  combine <- data.frame(scores, groups, dataset)
  results.friedman[setting, 2] <- friedman.test(scores ~ groups | dataset,
                                                data=combine)$p.value
}

# posthoc is  Wilcoxon, Nemenyi, McDonald-Thompson test
bs.results.wnmt  <- matrix(0, choose(nmethods, 2), nsettings)
auc.results.wnmt <- matrix(0, choose(nmethods, 2), nsettings)
for (setting in 1:nsettings) {
  scores <- as.vector(bs.results[[setting]])
  combine <- data.frame(scores, groups, dataset)
  bs.results.wnmt[, setting] <- pWNMT(x=combine$scores, b=combine$dataset,
                                         trt=combine$groups, n.mc=20000)$p.val

  scores <- as.vector(auc.results[[setting]])
  combine <- data.frame(scores, groups, dataset)
  auc.results.wnmt[, setting] <- pWNMT(x=combine$scores, b=combine$dataset,
                                          trt=combine$groups, n.mc=20000)$p.val

  print(paste("setting:", setting))
}

save(bs.results, bs.results.combined, bs.results.wnmt,
     auc.results, auc.results.combined, auc.results.wnmt,
     results.friedman,
     file = "results.RData")

# # look at a few iteration plots
# set <- 1
# setting <- 1
# dataset <- paste("sim-results/", setting, "-", set, ".RData", sep = "")
# load(dataset)
#
# par(mfrow = c(4, 5))
# for (i in 1:4) {
#   plot(log(fit.gev$a[, i, ]), type = "l",
#        main = paste("log(a[", i, "])", sep = ""))
# }
# plot(fit.gev$beta, type = "l", main = bquote(beta[0]))
#
# for (i in 8:11) {
#   plot(log(fit.gev$a[, i, ]), type = "l",
#        main = paste("log(a[", i, "])", sep = ""))
# }
# plot(fit.gev$alpha, type = "l", main = bquote(alpha))
#
# for (i in 1:4) {
#   plot(fit.gev$b[, i, ], type = "l", main = paste("b[", i, "]", sep = ""))
# }
# plot(fit.gev$rho, type = "l", main = bquote(rho))
#
# for (i in 6:10) {
#   plot(fit.gev$b[, i, ], type = "l", main = paste("b[", i, "]", sep = ""))
# }