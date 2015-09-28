rm(list = ls())

files <- list.files(path = "sim-tables/")

nsettings <- 4
nsets <- 100
nmethods <- 3
results <- vector(length = nsettings, mode = "list")
for (setting in 1:nsettings) {
  results[[setting]] <- matrix(NA, nsets, nmethods)
  rownames(results[[setting]]) <- paste("set", 1:100)
  colnames(results[[setting]]) <- c("gev", "probit", "logit")
}

for (i in 1:length(files)) {
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  setting   <- as.numeric(split[1])
  set       <- as.numeric(split[2])
  table.set <- read.table(paste("sim-tables/", files[i], sep = ""))
  if (nrow(table.set) == 2) {
    results[[setting]][set, 1] <- table.set[1, ]
    results[[setting]][set, 2] <- table.set[2, ]
  }
}

results.combined <- matrix(NA, nsettings, nmethods)
results.sd       <- matrix(NA, nsettings, nmethods)
for (setting in 1:nsettings) {
  results.combined[setting, ] <- apply(results[[setting]], 2, mean, 
                                       na.rm = TRUE)
  results.sd[setting, ] <- apply(results[[setting]], 2, sd, 
                                 na.rm = TRUE)
}

# how many have finished
colSums(!is.na(results[[1]]))
colSums(!is.na(results[[2]]))
colSums(!is.na(results[[3]]))
colSums(!is.na(results[[4]]))

setting.1 <- which(rowSums(is.na(results[[1]])) == 3)
setting.2 <- which(rowSums(is.na(results[[2]])) == 3)
setting.3 <- which(rowSums(is.na(results[[3]])) == 3)
setting.4 <- which(rowSums(is.na(results[[4]])) == 3)

# Check for differences
# First do Friedman test (one-way repeated measures)
#   friedman.test(y ~ trt | block, data)
# Then follow up with the Wilcoxon, Nemenyi, McDonald-Thompson test
# pWNMT(x, b, trt, method, n.mc)
#     x: list of values
#     b: vector of blocks (only needed if x is a vector)
#     trt: vector of treatments
#     method: "Exact", "Monte Carlo" or "Asymptotic"

groups <- rep(1:(nmethods - 1), each=100)
dataset <- rep(1:100, times=(nmethods - 1))

library(NSM3)
results.friedman <- rep(0, nsettings)
for (setting in 1:nsettings) {
  scores <- as.vector(results[[setting]][, c(1, 2)])
  combine <- data.frame(scores, groups, dataset)
  results.friedman[setting] <- friedman.test(scores ~ groups | dataset,
                                             data=combine)$p.value
}

# posthoc is  Wilcoxon, Nemenyi, McDonald-Thompson test
results.wnmt <- matrix(0, nrow = choose((nmethods - 1), 2), ncol = nsettings)
for (setting in 1:nsettings) {
  scores <- as.vector(results[[setting]][, c(1, 2)])
  combine <- data.frame(scores, groups, dataset)
  results.wnmt[, setting] <- pWNMT(x=combine$scores, b=combine$dataset,
                                   trt=combine$groups, n.mc=20000)$p.val
  print(paste("setting:", setting))
}
