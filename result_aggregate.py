import os

def aggregate_errors(input_dir):
    agg_file = "agg.txt"
    with open(agg_file, "a") as agg:
        for folder_name in os.listdir(input_dir):
            folder_path = os.path.join(input_dir, folder_name)
            if os.path.isdir(folder_path) and folder_name.startswith("0406"):
                errors_file = os.path.join(folder_path, "errors.txt")
                if os.path.isfile(errors_file):
                    with open(errors_file, "r") as ef:
                        lines = ef.readlines()
                        if lines:
                            agg.write(lines[-1])  # Append the last line of errors.txt

# Example usage
dir = input("Enter the directory path: ")
aggregate_errors(dir)
