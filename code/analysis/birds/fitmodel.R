# get the datasets
load(paste("./", species, ".RData", sep = ""))

# get the correct y, x, and s
if (species == "cattle_egret") {
  y <- cattle_egret[cv.idx[[cv]]]
} else if (species == "common_nighthawk") {
  y <- common_nighthawk[cv.idx[[cv]]]
} else if (species == "vesper_sparrow") {
  y <- vesper_sparrow[cv.idx[[cv]]]
} else if (species == "western_bluebird") {
  y <- western_bluebird[cv.idx[[cv]]]
} else {
  stop("incorrect species selected")
}

if (knot.percent == 10) {
  knots <- knots.10[[cv]]
} else if (knot.percent == 15) {
  knots <- knots.15[[cv]]
} else if (knot.percent == 20) {
  knots <- knots.20[[cv]]
} else {
  stop("only can use knots.percent = 10, 15, or 20")
}

# extract info about simulation settings
ns     <- length(y)
nt     <- 1
nknots <- nrow(knots)

y.o <- matrix(y, ns, nt)
s.o <- s[cv.idx[[cv]], ]
x.o <- matrix(1, nrow(s.o), 1)
s.p <- s[-cv.idx[[cv]], ]
x.p <- matrix(1, nrow(s.p), 1)
