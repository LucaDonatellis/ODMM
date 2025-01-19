#!/bin/bash

threads=(1 2 4 8 16 24 32 40 48 56 64 72 80 88 96)

for t in "${threads[@]}"; do
    export OMP_NUM_THREADS=$t
    echo "Running with $t threads:"
    ./test.sh $1
    echo
done
