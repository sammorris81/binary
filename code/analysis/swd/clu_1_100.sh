#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
SWD=$HOME/repos-git/rare-binary/code/analysis/swd
L=$1

echo "Cluster 1, Group $L"
$RBIN CMD BATCH --vanilla --no-save $SWD/fitmodel_clu_1_100_$L.R $SWD/fitmodel_clu_1_100_$L.out 2>&1

exit 0