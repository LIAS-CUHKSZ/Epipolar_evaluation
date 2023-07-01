$windows_size = $args[0]
$max_parallel = $args[1]
$dataset_path = "dataset"

# Get a list of all subdirectories in the dataset folder
# $subdirs = Get-ChildItem -Path $dataset_path -Directory

# only consider high res img
$subdirs = Get-ChildItem -Path $dataset_path -Directory

Write-Host ($subdirs.Count)
Set-Location -Path "build"

# Initialize the counter variable
$count = 0

# Loop through each subdirectory and run your program with the subdirectory name as an argument
foreach ($subdir in $subdirs)
{
    $subdir_name = $subdir.Name

    # Check if the program count is less than the max parallel
    if ($count -lt $max_parallel)
    {
        # Start a new program and increment the count
        Start-Process -FilePath "epipolar_eval.exe" -ArgumentList "$subdir_name $windows_size -NoNewWindow"
        $count++
    }
    else
    {
        # Wait for a program to finish and decrement the count
        Wait-Process -Id (Get-Process -Name "epipolar_eval" | Select-Object -First 1).Id
        $count--

        # Start a new program and increment the count
        Start-Process -FilePath "epipolar_eval.exe" -ArgumentList "$subdir_name $windows_size -NoNewWindow"
        $count++
    }
}

# Wait for all programs to finish
Wait-Process -Name "epipolar_eval"

Set-Location -Path ".."