#!/bin/bash

dataset_path="dataset"
subdirs=$(find $dataset_path -type d)
echo $(echo "$subdirs" | wc -l)
cd build

for subdir in $subdirs
do
    subdir_name=$(basename "$subdir")
        # Start a new program and increment the count
    ./epipolar_eval "$subdir_name" -1

done

## kitti_part
dts=("1_99_sift" "1_99_surf" 
         "2_99_sift" "2.0_99_surf" 
        "3_95_sift" "of2_2_99")

for dt in "${datadtssets[@]}"; do
    for seq in {00..10}; do
        ./kitti_eval $seq $dt
    done
done




