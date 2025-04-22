# rename
for dir in [0-9][0-9]; do if [ -d "$dir/evaltest" ]; then mv "$dir/evaltest" "$dir/1_99_orb"; echo "Renamed $dir/2_95 to $dir/3_95_sift"; fi; done

# delete
for dir in [0-9][0-9]; do if [ -d "$dir/2_95" ]; then rm -rf "$dir/2_95"; echo "Deleted $dir/2_95"; fi; done