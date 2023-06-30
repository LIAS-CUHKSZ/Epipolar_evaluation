$windows_size = $args[0]

$dataset_path = "dataset"

# Get a list of all subdirectories in the dataset folder
# $subdirs = Get-ChildItem -Path $dataset_path -Directory

# only consider high res img
$subdirs = Get-ChildItem $dataset_path | Where-Object { $_.PSIsContainer -and $_.Name -notlike "zlr*" }

echo ($subdirs | Measure-Object).Count
cd build
# Loop through each subdirectory and run your program with the subdirectory name as an argument
foreach ($subdir in $subdirs)
{
    $subdir_name = $subdir.Name
    Start-Process powershell.exe -ArgumentList "-Command", ".\epipolar_eval.exe $subdir_name $windows_size 1"
    Start-Sleep -Seconds 6
    # .\epipolar_eval.exe $subdir_name $windows_size
}
cd ..

# cd build
# .\ETH3DFormatLoader.exe facade
# .\ETH3DFormatLoader.exe delivery_area
# .\ETH3DFormatLoader.exe kicker
# .\ETH3DFormatLoader.exe courtyard
# .\ETH3DFormatLoader.exe electro
# .\ETH3DFormatLoader.exe meadow
# .\ETH3DFormatLoader.exe office
# .\ETH3DFormatLoader.exe pipes
# .\ETH3DFormatLoader.exe playground
# .\ETH3DFormatLoader.exe relief
# .\ETH3DFormatLoader.exe relief_2
# .\ETH3DFormatLoader.exe terrace
# .\ETH3DFormatLoader.exe terrains
# Start-Process powershell.exe -ArgumentList "-Command", ".\ETH3DFormatLoader.exe terrace $windows_size"
# Start-Process powershell.exe -ArgumentList "-Command", ".\ETH3DFormatLoader.exe terrains $windows_size"
# cd ..