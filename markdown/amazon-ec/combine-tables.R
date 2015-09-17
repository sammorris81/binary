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
