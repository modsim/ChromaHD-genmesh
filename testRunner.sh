#!/bin/bash

## TODO: Find a way to automatically compare results with a reference
set -e

TESTS=(box-1 box-4 cyl-1 cyl-4 per-xy per-xyz)
ROOT=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
TESTDIR=$ROOT/tests

OUT=out-test
EXT=vtk

MODE="RUN"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -c|--clean)
            MODE="CLEAN"; shift ;; 
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ $MODE == "RUN" ]]; then
    for TEST in ${TESTS[@]}; do 
        cd "$TESTDIR/$TEST"
        echo "Running test: $TEST"
        # "$ROOT/genmesh" genmesh.in $OUT.$EXT > $OUT.log 2>&1 & 
        "$ROOT/genmesh" genmesh.in $OUT.$EXT > $OUT.log 
    done
elif [[ $MODE == "CLEAN" ]]; then
    for TEST in ${TESTS[@]}; do 
        echo "Cleaning $TESTDIR/$TEST/$OUT"
        rm -rf $TESTDIR/$TEST/$OUT
        rm $TESTDIR/$TEST/$OUT.log
    done
fi
