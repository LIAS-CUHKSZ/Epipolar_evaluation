#!/bin/bash

max_parallel=12
count=0

for seq in {00..10}; do
    # Check if the program count is less than the max parallel
    if [ $count -lt $max_parallel ]; then
        # Start a new program and increment the count
        ./kitti_analysis $seq eval_optiflow_sac2_99 &
        count=$((count+1))
    else
        # Wait for a program to finish and decrement the count
        wait -n
        count=$((count-1))

        # Start a new program and increment the count
        ./kitti_analysis $seq eval_optiflow_sac2_99 &
        count=$((count+1))
    fi

done

# Wait for all programs to finish
wait


