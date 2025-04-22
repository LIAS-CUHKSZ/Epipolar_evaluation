#!/bin/bash
if [ $# != 5 ]
then
    echo "Usage: ./m-consistent_eval.sh <dataset> <max_parallel> <imgidx1> <imgidx2> <sample_time>"
    exit -1
fi

dataset=$1
max_parallel=$2
imgidx1=$3
imgidx2=$4
sample_time=$5
num_pts=(10 20 40 80 160 320 640 1280 2560 4000 -1)

cd build

# Initialize the counter variable
count=0

# Loop through each subdirectory and run your program with the subdirectory name as an argument
for num in "${num_pts[@]}"
do
    # Check if the program count is less than the max parallel
    if [ $count -lt $max_parallel ]
    then
        # Start a new program and increment the count
        ./MonteCarlo "$dataset" "$num" "$imgidx1" "$imgidx2" "$sample_time" & 
        count=$((count+1))
    else
        # Wait for a program to finish and decrement the count
        wait -n
        count=$((count-1))

        # Start a new program and increment the count
        ./MonteCarlo "$dataset" "$num" "$imgidx1" "$imgidx2" "$sample_time" & 
        count=$((count+1))
    fi
    # if you want to get a accurate time, please comment the above if-else statement and uncomment the following line
    # this wil run the program one by one
    # ./MonteCarlo.exe "$dataset" "$num" "$imgidx1" "$imgidx2" "$sample_time"
done

# Wait for all programs to finish
wait

cd ..
