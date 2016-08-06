rm(list = ls())
source("../../R/plotting.R", chdir = TRUE)
load("simdata-grid.RData")
tbl.dir <- "./sim-tables-2/"
res.dir <- "./sim-results-2/"
files <- list.files(path = tbl.dir)
library(ROCR)

ngen.methods <- 3
nsets <- 100
method.types <- c("gev", "pro", "log")
genmethod.names <- c("gev", "log", "hot")
nmethods <- length(method.types)
rareness <- matrix(NA, nsets, ngen.methods)
for (setting in 1:ngen.methods) {
  for (set in 1:nsets) {
    rareness[set, setting] <- mean(simdata[[setting]]$y.grid[, set])
  }
}

nsettings <- 12
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

setting.names <- rep(NA, nsettings)
for (i in 1:length(files)) {
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  gen.method <- as.numeric(split[3])
  if (split[1] == "clu") {
    samp.idx <- 1
  } else {
    samp.idx <- 2
  }
  if (as.numeric(split[2] == 100)) {
    n.idx <- 1
  } else {
    n.idx <- 2
  }
  set <- as.numeric(split[4])
  setting <- (gen.method - 1) * 4 + (n.idx - 1) * 2 + samp.idx
  setting.names[setting] <- paste(genmethod.names[gen.method], "-", split[1], 
                                  "-", split[2], sep = "")
  
  set       <- as.numeric(split[4])
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
bs.results.comb.se    <- matrix(NA, nsettings, nmethods)
bs.results.1.combined <- matrix(NA, nsettings, nmethods)
bs.results.0.combined <- matrix(NA, nsettings, nmethods)
auc.results.combined  <- matrix(NA, nsettings, nmethods)
auc.results.comb.se   <- matrix(NA, nsettings, nmethods)
# vector of sets numbers that are complete
finished.sets <- vector(mode = "list", length = nsettings)
for (setting in 1:nsettings) {
  these.sets <- which(rowSums(is.na(bs.results[[setting]])) == 0)
  finished.sets[[setting]] <- these.sets
  bs.results.combined[setting, ]   <- apply(bs.results[[setting]][these.sets, ],
                                            2, mean, na.rm = TRUE)
  bs.results.comb.se[setting, ]    <- apply(bs.results[[setting]][these.sets, ],
                                            2, sd, na.rm = TRUE) / sqrt(length(these.sets))
  bs.results.1.combined[setting, ] <- apply(bs.results.1[[setting]][these.sets, ],
                                           2, mean, na.rm = TRUE)
  bs.results.0.combined[setting, ] <- apply(bs.results.0[[setting]][these.sets, ],
                                            2, mean, na.rm = TRUE)
  auc.results.combined[setting, ]  <- apply(auc.results[[setting]][these.sets, ],
                                           2, mean, na.rm = TRUE)
  auc.results.comb.se[setting, ]   <- apply(auc.results[[setting]][these.sets, ],
                                            2, sd, na.rm = TRUE) / sqrt(length(these.sets))
}

rownames(bs.results.combined)  <- setting.names
rownames(auc.results.combined) <- setting.names
rownames(bs.results.comb.se)   <- setting.names
rownames(auc.results.comb.se)  <- setting.names
colnames(bs.results.combined)  <- method.types
colnames(auc.results.combined) <- method.types
colnames(bs.results.comb.se)   <- method.types
colnames(auc.results.comb.se)  <- method.types

round(bs.results.combined * 100, 2)
round(bs.results.comb.se * 100, 2)
round(auc.results.combined, 3)
round(auc.results.comb.se, 3)

df1 <- data.frame(Y = as.factor(simdata[[1]]$y.grid[, 16]), s1 = s.grid[, 1], 
                  s2 = s.grid[, 2])
p1 <- plot.species(df = df1, main = "Simulated GEV dataset")

df2 <- data.frame(Y = as.factor(simdata[[2]]$y.grid[, 23]), s1 = s.grid[, 1], 
                  s2 = s.grid[, 2])
p2 <- plot.species(df = df2, main = "Simulated probit dataset")

df3 <- data.frame(Y = as.factor(simdata[[3]]$y.grid[, 18]), s1 = s.grid[, 1], 
                  s2 = s.grid[, 2])
p3 <- plot.species(df = df3, main = "Simulated hotspot dataset")

panel <- arrangeGrob(p1, p2, p3, ncol = 3)
ggsave("./plots/simulateddata.pdf", panel, width = 17, height = 4.5)


for (setting in 1:nsettings) {
  print(length(finished.sets[[setting]]))
}

# Apparently ROCR can take in a list over all 100 datasets to
# come up with an averaged cross-validation curve
files <- list.files(path = res.dir)
samp.types <- c("clu", "srs")
ns <- c(100, 250)
for (gen.method in 1:3) {
  for (n.idx in 1:2) {
    for (samp.idx in 1:2) {
      this.setting <- (gen.method - 1) * 4 + (n.idx - 1) * 2 + samp.idx
      
      this.gev.pred <- vector(mode = "list", length = nsets)
      this.pro.pred <- vector(mode = "list", length = nsets)
      this.log.pred <- vector(mode = "list", length = nsets)
      this.yp       <- vector(mode = "list", length = nsets)
      sets.done <- rep(FALSE, nsets)
      
      for (set in 1:nsets) {
        results.file <- paste(res.dir, samp.types[samp.idx], "-", ns[n.idx], 
                              "-", gen.method, "-", set, ".RData", sep = "")
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
  }
}

#### look at over ROC curves and PRC curves
pred.gev.clu.100.1 <- prediction(pred.gev.1, yp.1)
pred.gev.srs.100.1 <- prediction(pred.gev.2, yp.2)
pred.pro.clu.100.1 <- prediction(pred.pro.1, yp.1)
pred.pro.srs.100.1 <- prediction(pred.pro.2, yp.2)
pred.log.clu.100.1 <- prediction(pred.log.1, yp.1)
pred.log.srs.100.1 <- prediction(pred.log.2, yp.2)

pred.gev.clu.250.1 <- prediction(pred.gev.3, yp.3)
pred.gev.srs.250.1 <- prediction(pred.gev.4, yp.4)
pred.pro.clu.250.1 <- prediction(pred.pro.3, yp.3)
pred.pro.srs.250.1 <- prediction(pred.pro.4, yp.4)
pred.log.clu.250.1 <- prediction(pred.log.3, yp.3)
pred.log.srs.250.1 <- prediction(pred.log.4, yp.4)

quartz(width = 16, height = 16)
par(mfrow = c(2, 2))
plot.roc(pred.gev.clu.100.1, pred.pro.clu.100.1, pred.log.clu.100.1,
         main = "Setting: GEV, Sample: CLU-100")
plot.roc(pred.gev.clu.250.1, pred.pro.clu.250.1, pred.log.clu.250.1,
         main = "Setting: GEV, Sample: CLU-250")
plot.roc(pred.gev.srs.100.1, pred.pro.srs.100.1, pred.log.srs.100.1,
         main = "Setting: GEV, Sample: SRS-100")
plot.roc(pred.gev.srs.250.1, pred.pro.srs.250.1, pred.log.srs.250.1,
         main = "Setting: GEV, Sample: SRS-250")
dev.print(device = pdf, "./plots/sim-perf-gev.pdf")
dev.off()

pred.gev.clu.100.2 <- prediction(pred.gev.5, yp.5)
pred.gev.srs.100.2 <- prediction(pred.gev.6, yp.6)
pred.pro.clu.100.2 <- prediction(pred.pro.5, yp.5)
pred.pro.srs.100.2 <- prediction(pred.pro.6, yp.6)
pred.log.clu.100.2 <- prediction(pred.log.5, yp.5)
pred.log.srs.100.2 <- prediction(pred.log.6, yp.6)

pred.gev.clu.250.2 <- prediction(pred.gev.7, yp.7)
pred.gev.srs.250.2 <- prediction(pred.gev.8, yp.8)
pred.pro.clu.250.2 <- prediction(pred.pro.7, yp.7)
pred.pro.srs.250.2 <- prediction(pred.pro.8, yp.8)
pred.log.clu.250.2 <- prediction(pred.log.7, yp.7)
pred.log.srs.250.2 <- prediction(pred.log.8, yp.8)

quartz(width = 16, height = 16)
par(mfrow = c(2, 2))
plot.roc(pred.gev.clu.100.2, pred.pro.clu.100.2, pred.log.clu.100.2,
         main = "Setting: Probit, Sample: CLU-100")
plot.roc(pred.gev.clu.250.2, pred.pro.clu.250.2, pred.log.clu.250.2,
         main = "Setting: Probit, Sample: CLU-250")
plot.roc(pred.gev.srs.100.2, pred.pro.srs.100.2, pred.log.srs.100.2,
         main = "Setting: Probit, Sample: SRS-100")
plot.roc(pred.gev.srs.250.2, pred.pro.srs.250.2, pred.log.srs.250.2,
         main = "Setting: Probit, Sample: SRS-250")
dev.print(device = pdf, "./plots/sim-perf-probit.pdf")
dev.off()

pred.gev.clu.100.3 <- prediction(pred.gev.9, yp.9)
pred.gev.srs.100.3 <- prediction(pred.gev.10, yp.10)
pred.pro.clu.100.3 <- prediction(pred.pro.9, yp.9)
pred.pro.srs.100.3 <- prediction(pred.pro.10, yp.10)
pred.log.clu.100.3 <- prediction(pred.log.9, yp.9)
pred.log.srs.100.3 <- prediction(pred.log.10, yp.10)

pred.gev.clu.250.3 <- prediction(pred.gev.11, yp.11)
pred.gev.srs.250.3 <- prediction(pred.gev.12, yp.12)
pred.pro.clu.250.3 <- prediction(pred.pro.11, yp.11)
pred.pro.srs.250.3 <- prediction(pred.pro.12, yp.12)
pred.log.clu.250.3 <- prediction(pred.log.11, yp.11)
pred.log.srs.250.3 <- prediction(pred.log.12, yp.12)

quartz(width = 16, height = 16)
par(mfrow = c(2, 2))
plot.roc(pred.gev.clu.100.3, pred.pro.clu.100.3, pred.log.clu.100.3,
         main = "Setting: Hotspot, Sample: CLU-100")
plot.roc(pred.gev.clu.250.3, pred.pro.clu.250.3, pred.log.clu.250.3,
         main = "Setting: Hotspot, Sample: CLU-250")
plot.roc(pred.gev.srs.100.3, pred.pro.srs.100.3, pred.log.srs.100.3,
         main = "Setting: Hotspot, Sample: SRS-100")
plot.roc(pred.gev.srs.250.3, pred.pro.srs.250.3, pred.log.srs.250.3,
         main = "Setting: Hotspot, Sample: SRS-250")
dev.print(device = pdf, "./plots/sim-perf-hotspot.pdf")
dev.off()

#### smooth of rareness by brier score
for (setting in 1:nsettings) {
  this.grid <- ceiling(setting / 4)
  these.include <- which(rowSums(!is.na(bs.results[[setting]])) == 3)
  plot.file <- paste("./plots/byrareness-", setting, ".pdf", sep = "")
  df <- data.frame(rareness = rep(rareness[these.include, this.grid], 3),
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
  ggsave(plot.file, plot = panel, width = 9, height = 4.5)
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
nsets <- 5

factor.1 <- as.factor(c(100, 250))
factor.2 <- as.factor(c("clu", "srs"))
factor.3 <- as.factor(c("gev", "probit", "logistic"))

a <- length(factor.1)
b <- length(factor.2)
c <- length(factor.3)

sample.size <- rep(factor.1, each = b * c * nsets)
sample.type <- rep(factor.2, each = c * nsets, times = a)
fit.method  <- rep(factor.3, each = nsets, times = a * b)
dataset     <- as.factor(rep(1:nsets, times = a * b * c))
results.friedman <- matrix(0, 3, 2)  # one row for each of gev, pro, log
colnames(results.friedman) <- c("bs", "auc")

scores.bs <- c(as.vector(bs.results[[1]][1:nsets, ]),
               as.vector(bs.results[[2]][1:nsets, ]),
               as.vector(bs.results[[3]][1:nsets, ]),
               as.vector(bs.results[[4]][1:nsets, ]))
scores.auc <- c(as.vector(auc.results[[1]][1:nsets, ]),
                as.vector(auc.results[[2]][1:nsets, ]),
                as.vector(auc.results[[3]][1:nsets, ]),
                as.vector(auc.results[[4]][1:nsets, ]))

boxplot(scores.bs ~ sample.size * sample.type * fit.method)
boxplot(scores.auc ~ sample.size * sample.type * fit.method)

combine.bs <- data.frame(scores.bs, sample.size, sample.type, fit.method, dataset)
lme.1.bs.full <- lmer(scores.bs ~ sample.size*sample.type*fit.method + 
                        (1 | dataset), data = combine.bs, REML = FALSE)
lme.1.bs.red1 <- lmer(scores.bs ~ (sample.size + sample.type + fit.method)^2 + 
                        (1 | dataset), data = combine.bs, REML = FALSE)
lme.1.bs.red2<- lmer(scores.bs ~ sample.size*sample.type + fit.method + 
                       (1 | dataset), data = combine.bs, REML = FALSE)
lme.1.bs.red3 <- lmer(scores.bs ~ sample.size + sample.type + fit.method + 
                        (1 | dataset), data = combine.bs, REML = FALSE)
lme.1.bs.red4 <- lmer(scores.bs ~ sample.size + sample.type + 
                       (1 | dataset), data = combine.bs, REML = FALSE)
lme.1.bs.red5 <- lmer(scores.bs ~ sample.size + (1 | dataset), data = combine.bs, REML = FALSE)
BIC(lme.1.bs.full, lme.1.bs.red1, lme.1.bs.red2, lme.1.bs.red3, lme.1.bs.red4)

combine.auc <- data.frame(scores.auc, sample.size, sample.type, fit.method, dataset)
lme.1.auc <- lmer(scores.auc ~ sample.size + sample.type + fit.method + 
                    (1 | dataset), data = combine.auc)

results.friedman[1, 1] <- friedman.test(scores ~ sample.size + sample.type + 
                                          fit.method | dataset, 
                                        data = combine)

for (setting in 1:nsettings) {

  
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