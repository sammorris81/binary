rm(list = ls())

prefix <- "cv-tables/"
files <- list.files(path = prefix)
knot.percents <- c(10, 15, 20)
cvs <- c(1, 2)
species.list <- c("bluewinged_teal", "cattle_egret", "common_grounddove",
                  "common_nighthawk", "greater_white_goose", "hooded_oriole",
                  "lesser_goldfinch", "longbilled_curlew", "longeared_owl",
                  "mountain_bluebird", "piping_plover", "sharpshinned_hawk",
                  "snowy_plover", "vesper_sparrow", "western_bluebird",
                  "whiteeyed_vireo"
                  )

nmethods <- 3
nsets <- length(knot.percents) * length(cvs)

bs.results <- auc.results <- vector(length = length(species.list),
                                    mode = "list")
for (species in species.list) {
  species.idx <- which(species.list == species)
  bs.results[[species.idx]]  <- matrix(NA, nsets, nmethods)
  auc.results[[species.idx]] <- matrix(NA, nsets, nmethods)
  rownames(bs.results[[species.idx]])  <- c(paste("knots-10-", 1:2, sep = ""),
                                            paste("knots-15-", 1:2, sep = ""),
                                            paste("knots-20-", 1:2, sep = ""))
  rownames(auc.results[[species.idx]]) <- c(paste("knots-10-", 1:2, sep = ""),
                                            paste("knots-15-", 1:2, sep = ""),
                                            paste("knots-20-", 1:2, sep = ""))
  colnames(bs.results[[species.idx]])  <- c("gev", "probit", "logit")
  colnames(auc.results[[species.idx]]) <- c("gev", "probit", "logit")
}

for (i in 1:length(files)) {
  split        <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  species.idx  <- which(species.list == split[1])
  knot.percent <- as.numeric(split[2])
  cv           <- as.numeric(split[3])
  this.row     <- (which(knot.percents == knot.percent) - 1) * 2 + cv
  table.set    <- read.table(paste(prefix, files[i], sep = ""))
  bs.results[[species.idx]][this.row, 1]  <- table.set[1, 1]
  auc.results[[species.idx]][this.row, 1] <- table.set[1, 2]
  bs.results[[species.idx]][this.row, 2]  <- table.set[2, 1]
  auc.results[[species.idx]][this.row, 2] <- table.set[2, 2]
  bs.results[[species.idx]][this.row, 3]  <- table.set[3, 1]
  auc.results[[species.idx]][this.row, 3] <- table.set[3, 2]
}

bs.results.combined  <- vector(mode = "list", length = length(species))
auc.results.combined <- vector(mode = "list", length = length(species))
species.rate <- rep(NA, length(species))
for (species in species.list) {
  species.idx <- which(species.list == species)
  load(paste(species, ".RData", sep = ""))
  species.rate[species.idx] <- sum(get(species))
  bs.results.combined[[species.idx]] <- matrix(NA, length(knot.percents),
                                               nmethods)
  auc.results.combined[[species.idx]] <- matrix(NA, length(knot.percents),
                                                nmethods)
  rownames(bs.results.combined[[species.idx]])  <- paste("knots-",
                                                         knot.percents,
                                                         sep = "")
  rownames(auc.results.combined[[species.idx]]) <- paste("knots-",
                                                         knot.percents,
                                                         sep = "")
  colnames(bs.results.combined[[species.idx]])  <- c("gev", "probit", "logit")
  colnames(auc.results.combined[[species.idx]]) <- c("gev", "probit", "logit")

  this.row <- apply(bs.results[[species.idx]][1:2, ], 2, mean, na.rm = TRUE)
  bs.results.combined[[species.idx]][1, ] <- this.row
  this.row <- apply(bs.results[[species.idx]][3:4, ], 2, mean, na.rm = TRUE)
  bs.results.combined[[species.idx]][2, ] <- this.row
  this.row <- apply(bs.results[[species.idx]][5:6, ], 2, mean, na.rm = TRUE)
  bs.results.combined[[species.idx]][3, ] <- this.row

  this.row <- apply(auc.results[[species.idx]][1:2, ], 2, mean, na.rm = TRUE)
  auc.results.combined[[species.idx]][1, ] <- this.row
  this.row <- apply(auc.results[[species.idx]][3:4, ], 2, mean, na.rm = TRUE)
  auc.results.combined[[species.idx]][2, ] <- this.row
  this.row <- apply(auc.results[[species.idx]][5:6, ], 2, mean, na.rm = TRUE)
  auc.results.combined[[species.idx]][3, ] <- this.row
}

specied.list[5]
