$windows_size = 2

$dataset_path = "dataset"

# Get a list of all subdirectories in the dataset folder
$subdirs = Get-ChildItem -Path $dataset_path -Directory
cd build
# Loop through each subdirectory and run your program with the subdirectory name as an argument
foreach ($subdir in $subdirs)
{
    $subdir_name = $subdir.Name
    Start-Process powershell.exe -ArgumentList "-Command", ".\ETH3DFormatLoader.exe $subdir_name $windows_size"
}

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