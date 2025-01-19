#!/bin/bash
nps=(1 2 3 4 9 16)

for np in "${nps[@]}"; do
    echo "Running with $np processors:"
    for i in {4..10}; do
        n=$((2**i))
        echo -ne "$1 $n:\t"
        sum=0
        for j in {1..10}; do
            time_output=$(mpirun -np $np ./$1 $n | grep -oP 'Execution time: \K[0-9.]+')
            echo -n " $time_output"
            sum=$(echo "$sum $time_output" | awk '{print $1 + $2}')
        done
        avg=$(echo "$sum 10" | awk '{print $1 / $2}')
        echo " -> $avg"
    done
    echo
done