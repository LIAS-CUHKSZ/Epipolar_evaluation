#!/bin/bash

for seq in {00..10}; do
    ./kitti_eval $seq eval &
    ./kitti_eval $seq evaltest &
    wait
done


