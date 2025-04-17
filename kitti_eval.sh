#!/bin/bash

# Ensure all background tasks are killed if the script is terminated
trap "echo 'Terminating all background tasks...'; kill 0" SIGINT SIGTERM EXIT

# Set max_parallel to the number of logical cores
max_parallel=$(nproc)

datasets=("1_99_sift" "1_99_surf" 
         "2_99_sift" "2.0_99_surf" 
        "3_95_sift" "of2_2_99")

cd build 

count=0

for dataset in "${datasets[@]}"; do
    
    for seq in {00..10}; do
        if [ $count -lt $max_parallel ]; then
            ./kitti_eval $seq $dataset &
            count=$((count+1))
        else
            wait -n
            count=$((count-1))
            ./kitti_eval $seq $dataset &
            count=$((count+1))
        fi
    done
done

wait

# mv all results to the results folder


cd ..


