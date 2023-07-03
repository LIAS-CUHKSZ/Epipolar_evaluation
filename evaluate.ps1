$windows_size = $args[0]
$max_parallel = $args[1]
$dataset_path = "dataset"

# Get a list of all subdirectories in the dataset folder
# $subdirs = Get-ChildItem -Path $dataset_path -Directory

# only consider high res img
# we add a prefix "zlr" to the low res img folders
$subdirs = Get-ChildItem $dataset_path | Where-Object { $_.PSIsContainer -and $_.Name -notlike "zlr*" }
# $subdirs = Get-ChildItem $dataset_path # this will evaluate all scenes

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
        Start-Process powershell.exe -ArgumentList "./epipolar_eval.exe $subdir_name $windows_size"
        $count++
        sleep 1
    }
    while($count -ge $max_parallel)
    {
        $running_jobs = (Get-Process -Name "epipolar_eval" ).Count
        if($running_jobs -lt $max_parallel -and $running_jobs -gt 0)
        {
            $count--
        }
        sleep 1
    }
}

Set-Location -Path ".."