rm(list = ls())

files <- list.files(path = "sim-tables/")

nsettings <- 6
nsets <- 100
nmethods <- 3
bs.results <- auc.results <- vector(length = nsettings, mode = "list")
for (setting in 1:nsettings) {
  bs.results[[setting]]  <- matrix(NA, nsets, nmethods)
  auc.results[[setting]] <- matrix(NA, nsets, nmethods)
  rownames(bs.results[[setting]])  <- paste("set", 1:nsets)
  rownames(auc.results[[setting]]) <- paste("set", 1:nsets)
  colnames(bs.results[[setting]])  <- c("gev", "probit", "logit")
  colnames(auc.results[[setting]]) <- c("gev", "probit", "logit")
}

for (i in 1:length(files)) {
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  setting   <- as.numeric(split[1])
  set       <- as.numeric(split[2])
  table.set <- read.table(paste("sim-tables/", files[i], sep = ""))
  bs.results[[setting]][set, 1]  <- table.set[1, 1]
  auc.results[[setting]][set, 1] <- table.set[1, 2]
  bs.results[[setting]][set, 2]  <- table.set[2, 1]
  auc.results[[setting]][set, 2] <- table.set[2, 2]
  bs.results[[setting]][set, 3]  <- table.set[3, 1]
  auc.results[[setting]][set, 3] <- table.set[3, 2]
}

bs.results.combined  <- matrix(NA, nsettings, nmethods)
auc.results.combined <- matrix(NA, nsettings, nmethods)
# vector of sets numbers that are complete
finished.sets <- vector(mode = "list", length = nsettings)
for (setting in 1:nsettings) {
  these.sets <- which(rowSums(is.na(bs.results[[setting]])) == 0)
  finished.sets[[setting]] <- these.sets
  bs.results.combined[setting, ]  <- apply(bs.results[[setting]][these.sets, ],
                                           2, mean, na.rm = TRUE)
  auc.results.combined[setting, ] <- apply(auc.results[[setting]][these.sets, ],
                                           2, mean, na.rm = TRUE)
}

for (setting in 1:nsettings) {
  print(length(finished.sets[[setting]]))
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