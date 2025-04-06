#!/bin/bash
if [ $# -ne 2 ]
then
    echo "Usage: $0 <windows_size> <max_parallel>"
    exit 1
fi

windows_size=$1
max_parallel=$2
dataset_path="dataset"

subdirs=$(find $dataset_path -type d)

echo $(echo "$subdirs" | wc -l)
cd build

# Initialize the counter variable
count=0

# Loop through each subdirectory and run your program with the subdirectory name as an argument
for subdir in $subdirs
do
    subdir_name=$(basename "$subdir")

    # Check if the program count is less than the max parallel
    if [ $count -lt $max_parallel ]
    then
        # Start a new program and increment the count
        ./epipolar_eval "$subdir_name" "$windows_size" &
        count=$((count+1))
    else
        # Wait for a program to finish and decrement the count
        wait -n
        count=$((count-1))

        # Start a new program and increment the count
        ./epipolar_eval "$subdir_name" "$windows_size" &
        count=$((count+1))
    fi
done

# Wait for all programs to finish
wait

cd ..