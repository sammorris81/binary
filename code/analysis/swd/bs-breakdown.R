rm(list = ls())

ns <- c(100, 250)
samp.types <- c("clu", "srs")

nmethods <- 3
nsets <- 100
nfits <- nsets * length(ns) * length(samp.types)

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

