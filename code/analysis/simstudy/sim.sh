#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
SIM=$HOME/repos-git/rare-binary/code/analysis/simstudy
L=$1
COPY=$2

echo "Setting: $L, Copy: $COPY"
$RBIN CMD BATCH --vanilla --no-save $SIM/sim-$L.R $SIM/sim-$L-$COPY.out 2>&1

exit 0