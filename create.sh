#!/bin/bash

mkdir -p output/${2%.*}/
nohup ./genmesh $1 $2 > output/${2%.*}/${2%.*}.log 2>&1 & disown
