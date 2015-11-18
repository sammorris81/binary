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
  results[[setting]][set, 1] <- table.set[1, ]
  results[[setting]][set, 2] <- table.set[2, ]
  results[[setting]][set, 3] <- table.set[3, ]
}

results.combined <- matrix(NA, nsettings, nmethods)
for (setting in 1:nsettings) {
  results.combined[setting, ] <- apply(results[[setting]], 2, mean, 
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

round(results.combined[, 1] / results.combined[, 3], 4)
round(results.combined[, 2] / results.combined[, 3], 4)

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