#!/bin/bash

# Check if all three arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <variable_1> <input_dir> <output_dir>"
    exit 1
fi

# Assign variables
A="$1"
input_dir="$2"
output_dir="$3"

# Check if input directory exists
if [ ! -d "$input_dir" ]; then
    echo "Input directory '$input_dir' does not exist."
    exit 1
fi

# Check if output directory exists, if not, create it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# Create an array to store the filenames
files=()

# Loop through files in the folder containing variable A
# while IFS= read -r -d '' file; do
#     if [[ $file == *".bam" ]]; then
#         files+=("$file")
#     fi
# done < <(find "$folder_in" -type f -name "*$A*.bam" -print0)
while IFS= read -r -d '' file; do
    files+=("$file")
done < <(find "$input_dir" -type f -name "*$A*.bam" -print0)

# Check if any files containing variable A were found
if [ ${#files[@]} -eq 0 ]; then
    echo "No files containing '$A' found in $input_dir"
    exit 1
fi

# Create the output filename
output_filename="${output_dir}/bam_${A}.txt"

# Concatenate all filenames separated by comma and save to the output file
printf "%s," "${files[@]}" | sed 's/,$//' > "$output_filename"
printf "\n" >> "$output_filename"

echo "List of files containing '$A' saved to '$output_filename'"