#!/bin/bash

nohup ./mesher $1 > output/${1%.*}.log 2>&1 & disown
