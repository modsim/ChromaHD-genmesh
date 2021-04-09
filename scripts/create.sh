#!/bin/bash

###########################3
## Deprecated!!!!
###########################3

# mkdir -p output/${2%.*}/
# nohup ./genmesh $1 $2 > output/${2%.*}/${2%.*}.log 2>&1 & disown

TITLE="${2%.*}"
DIR="output/$TITLE"
CFG="$1"
LOG="$DIR/$TITLE.log"

mkdir -p "$DIR"
nohup ./genmesh "$CFG" "$2" > "$LOG" 2>&1 & disown
sleep 2
watch -n 1 tail -n $(tput lines) "$LOG"
