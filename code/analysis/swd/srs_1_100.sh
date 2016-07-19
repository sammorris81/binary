#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
SWD=$HOME/repos-git/rare-binary/code/analysis/swd
L=$1

echo "SRS 1, Group $L"
$RBIN CMD BATCH --vanilla --no-save $SWD/fitmodel_srs_1_100_$L.R $SWD/fitmodel_srs_1_100_$L.out 2>&1

exit 0