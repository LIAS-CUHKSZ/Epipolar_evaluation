$num_pts = @(10,20,40,80,160,320,640,1280,2560,3000,-1)

$dataset_path = $args[0]
$windows_size = $args[1] 
$absolute_path = Convert-Path dataset/$dataset_path

cd build
# Loop through each subdirectory and run your program with the subdirectory name as an argument
foreach ($item in $num_pts)
{
    .\MonteCarlo.exe $dataset_path  $windows_size $item
}
cd ..
