rm(list = ls())

prefix <- "cv-tables-samp/"
files <- list.files(path = prefix)
ns <- c(100, 200)
samp.types <- c("clu", "srs")
species.list <- c("broadbill_hummingbird", "cattle_egret", "common_grounddove",
                  "eurasian_wigeon", "greater_white_goose", "hooded_oriole",
                  "lapland_longspur", "longbilled_curlew", "longeared_owl",
                  "mountain_bluebird", "northern_sawwhet_owl", "piping_plover",
                  "snowy_plover", "tricolored_blackbird", "vesper_sparrow")

nmethods <- 3
nsets <- 50
nfits <- nsets * length(ns) * length(samp.types)

bs.results <- auc.results <- vector(length = length(species.list),
                                    mode = "list")
for (species in species.list) {
  species.idx <- which(species.list == species)
  bs.results[[species.idx]]  <- matrix(NA, nfits, nmethods)
  auc.results[[species.idx]] <- matrix(NA, nfits, nmethods)
  these.rownames <- c(paste("clu-100-", 1:nsets, sep = ""),
                      paste("srs-100-", 1:nsets, sep = ""),
                      paste("clu-200-", 1:nsets, sep = ""),
                      paste("srs-200-", 1:nsets, sep = ""))
  these.colnames <- c("gev", "probit", "logit")
  rownames(bs.results[[species.idx]])  <- these.rownames
  rownames(auc.results[[species.idx]]) <- these.rownames
  colnames(bs.results[[species.idx]])  <- these.colnames
  colnames(auc.results[[species.idx]]) <- these.colnames
}

for (i in 1:length(files)) {
  split         <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  species.idx   <- which(species.list == split[1])
  samp.type.idx <- which(samp.types == split[2])
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

bs.results.combined  <- vector(mode = "list", length = length(species))
auc.results.combined <- vector(mode = "list", length = length(species))
species.rate <- rep(NA, length(species))
for (species in species.list) {
  species.idx <- which(species.list == species)
  load(paste(species, ".RData", sep = ""))
  species.rate[species.idx] <- mean(get(species))
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

specied.list[5]
