rm(list = ls())

files <- list.files(path = "sim-tables/")

nsettings <- 3
nsets <- 20
nmethods <- 9
bs.results <- auc.results <- vector(length = nsettings, mode = "list")
for (setting in 1:nsettings) {
  bs.results[[setting]]  <- matrix(NA, nsets, nmethods)
  auc.results[[setting]] <- matrix(NA, nsets, nmethods)
  rownames(bs.results[[setting]])  <- paste("set", 1:nsets)
  rownames(auc.results[[setting]]) <- paste("set", 1:nsets)
  colnames(bs.results[[setting]])  <- c("gev-1", "gev-2", "gev-3",
                                        "probit-1", "probit-2", "probit-3",
                                        "logit-1", "logit-2", "logit-3")
  colnames(auc.results[[setting]]) <- c("gev-1", "gev-2", "gev-3",
                                        "probit-1", "probit-2", "probit-3",
                                        "logit-1", "logit-2", "logit-3")
}

for (i in 1:length(files)) {
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  setting   <- as.numeric(split[1])
  set       <- as.numeric(split[2])
  knots     <- as.numeric(split[3])
  table.set <- read.table(paste("sim-tables/", files[i], sep = ""))
  if (knots == 1) {
    bs.results[[setting]][set, 1]  <- table.set[1, 1]
    auc.results[[setting]][set, 1] <- table.set[1, 2]
    bs.results[[setting]][set, 4]  <- table.set[4, 1]
    auc.results[[setting]][set, 4] <- table.set[4, 2]
    bs.results[[setting]][set, 7]  <- table.set[7, 1]
    auc.results[[setting]][set, 7] <- table.set[7, 2]
  } else if (knots == 2) {
    bs.results[[setting]][set, 2]  <- table.set[2, 1]
    auc.results[[setting]][set, 2] <- table.set[2, 2]
    bs.results[[setting]][set, 5]  <- table.set[5, 1]
    auc.results[[setting]][set, 5] <- table.set[5, 2]
    bs.results[[setting]][set, 8]  <- table.set[8, 1]
    auc.results[[setting]][set, 8] <- table.set[8, 2]
  } else {
    bs.results[[setting]][set, 3]  <- table.set[3, 1]
    auc.results[[setting]][set, 3] <- table.set[3, 2]
    bs.results[[setting]][set, 6]  <- table.set[6, 1]
    auc.results[[setting]][set, 6] <- table.set[6, 2]
    bs.results[[setting]][set, 9]  <- table.set[9, 1]
    auc.results[[setting]][set, 9] <- table.set[9, 2]
  }
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

settings <- c("GEV Link - Knots on 21 x 21 grid",
              "GEV Link - Knots at 441 randomly selected sites",
              "Logit Link")

for (setting in 1:nsettings) {
  these.sets <- finished.sets[[setting]]
  boxplot(bs.results[[setting]][these.sets, ], 
          main = paste("Brier scores for ", settings[setting], ", ", 
                       length(these.sets), " sets finished", sep = ""))
  points(1:9, bs.results.combined[setting, ], pch = 21, 
         bg = "firebrick1", col = "firebrick4")
  dev.print(device = pdf, paste("plots/bs-", setting, ".pdf", sep = ""))
  
  boxplot(auc.results[[setting]][finished.sets[[setting]], ], 
          main = paste("AUC for ", settings[setting], ", ", 
                       length(these.sets), " sets finished", sep = ""))
  points(1:9, auc.results.combined[setting, ], pch = 21, 
         bg = "firebrick1", col = "firebrick4")
  dev.print(device = pdf, paste("plots/auc-", setting, ".pdf", sep = ""))
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