#!/bin/bash

nohup ./mesher $1 $2 > output/${2%.*}.log 2>&1 & disown
