#!/bin/bash
# CWD=`pwd`
# proper syntax is filename start end by setting
OUT=`echo $1 | cut -f1 -d.`
EXT=".R"
RDATA=".RData"
ORIG="$OUT$EXT"
START=$2
SET=$2
END=$3
BY=$4
SETTING=$5
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
  MCMC="$OUT-$SETTING-$INCLUDE-bw$EXT"
  OUTPUT="$OUT$GROUP$RDATA"
  echo "OUTPUT: $OUTPUT"
  echo "INCLUDE: $INCLUDE"
  cp "$1" "$MCMC"
  sed "23 a\sets <- $INCLUDE" "$ORIG" > "$MCMC.tmp"
  sed "24 a\setting <- $SETTING" "$MCMC.tmp" > "$MCMC"
  rm "$MCMC.tmp"
  sleep 60
  bwsubmit r "$MCMC"
  ((GROUP++))
done
