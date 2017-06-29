#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
SIM=$HOME/repos-git/rare-binary/code/analysis/simstudy
SAMP=$1
N=$2
GENMETH=$3
COPY=$4


echo "Samp: $SAMP, N: $N, Generation Method: $GENMETH, Copy: $COPY"
$RBIN CMD BATCH --vanilla --no-save $SIM/sim-$SAMP-$N-$GENMETH.R $SIM/sim-$SAMP-$N-$GENMETH-$COPY.out 2>&1

exit 0