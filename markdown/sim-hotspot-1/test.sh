#!/bin/sh
# CWD=`pwd`
OUT=`echo $1 | cut -f1 -d.`
EXT=".R"
RDATA=".RData"
ORIG="$OUT$EXT"
START=$2
SET=$2
END=$3
BY=$4
GROUP=1
while [ "$SET" -le "$END" ]
do
  if [ "$BY" -eq 1 ]
  then
    echo "option 1"
    INCLUDE="$SET"
    ((SET++))
  elif [ "$BY" -eq 2 ]
  then
    echo "option 2"
    if [ $((SET + 1)) -gt "$END" ]
    then
      INCLUDE="c($SET)"
    else
      INCLUDE="c($SET, $((SET + 1)))"
    fi
    SET=$((SET + BY))
  else
    GROUPSTART=$(((GROUP - 1) * BY + 2))
    GROUPEND=$((GROUP * BY))
    if [ "$GROUPEND" -gt "$END" ]
    then
      GROUPEND="$END"
    fi
    INCLUDE="c($SET"
    for ((i=$GROUPSTART; i<=$GROUPEND; i++))
    do
      INCLUDE="$INCLUDE,$i"
    done
    INCLUDE="$INCLUDE)"

    if [ $((GROUPSTART - 1)) -eq "$END" ]
    then
      INCLUDE="$END"
    fi
    SET=$((SET + BY))
  fi
  MCMC="$OUT$GROUP$EXT"
  OUTPUT="$OUT$GROUP$RDATA"
  echo "OUTPUT: $OUTPUT"
  echo "INCLUDE: $INCLUDE"
  cp "$1" "$MCMC"
  sed "1 a\ filename <- \"$OUTPUT\"" "$ORIG" > "$MCMC.tmp"
  sed "2 a\ sets <- $INCLUDE" "$MCMC.tmp" > "$MCMC"
  rm "$MCMC.tmp"
  bwsubmit_multi 2 r "$MCMC"
  ((GROUP++))
done