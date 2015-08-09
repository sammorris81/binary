#!/bin/sh
# CWD=`pwd`
OUT=`echo $1 | cut -f1 -d.`.out
R CMD BATCH --vanilla --no-save $1 $OUT & 2>&1