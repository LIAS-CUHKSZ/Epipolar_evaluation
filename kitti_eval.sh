#!/bin/bash

max_parallel=12
count=0

for seq in {00..10}; do
    # Check if the program count is less than the max parallel
    if [ $count -lt $max_parallel ]; then
        # Start a new program and increment the count
        ./kitti_eval $seq eval &
        count=$((count+1))
    else
        # Wait for a program to finish and decrement the count
        wait -n
        count=$((count-1))

        # Start a new program and increment the count
        ./kitti_eval $seq eval &
        count=$((count+1))
    fi

    # Repeat for other tasks
    if [ $count -lt $max_parallel ]; then
        ./kitti_eval $seq evaltest &
        count=$((count+1))
    else
        wait -n
        count=$((count-1))
        ./kitti_eval $seq evaltest &
        count=$((count+1))
    fi

    if [ $count -lt $max_parallel ]; then
        ./kitti_eval $seq 2_95 &
        count=$((count+1))
    else
        wait -n
        count=$((count-1))
        ./kitti_eval $seq 2_95 &
        count=$((count+1))
    fi
done

# Wait for all programs to finish
wait


