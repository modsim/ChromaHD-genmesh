#!/bin/bash

nohup ./genmesh $1 $2 > logs/${2%.*}.log 2>&1 & disown
