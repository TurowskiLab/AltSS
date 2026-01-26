#!/bin/bash

# Check if all four arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <variable_1> <input_dir> <output_dir> <tool_type>"
    exit 1
fi

# Assign variables
A="$1"
input_dir="$2"
output_dir="$3"
tool_type="$4"

# Validate tool_type
if [ "$tool_type" != "rmats" ] && [ "$tool_type" != "spladder" ]; then
    echo "Error: tool_type must be either 'rmats' or 'spladder', got '$tool_type'"
    exit 1
fi

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

while IFS= read -r -d '' file; do
    files+=("$file")
done < <(find "$input_dir" -type f -name "*$A*.bam" ! -name "*STAR*" -print0)

# Check if any files containing variable A were found
if [ ${#files[@]} -eq 0 ]; then
    echo "No files containing '$A' found in $input_dir"
    exit 1
fi

# Create the output filename
output_filename="${output_dir}/bam_${A}.txt"

# Concatenate based on tool_type
if [ "$tool_type" = "rmats" ]; then
    # For rmats: concatenate separated by comma
    printf "%s," "${files[@]}" | sed 's/,$//' > "$output_filename"
    printf "\n" >> "$output_filename"
elif [ "$tool_type" = "spladder" ]; then
    # For spladder: concatenate in lines (one file per line)
    printf "%s\n" "${files[@]}" > "$output_filename"
else
    echo "Error: Unsupported tool_type '$tool_type'"
    exit 1
fi

echo "List of files containing '$A' saved to '$output_filename'"