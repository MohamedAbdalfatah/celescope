#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory_path> <sample_name>"
    exit 1
fi

# Assign inputs to variables
directory_path="$1"
sample_name="$2"

# Check if the directory exists
if [ ! -d "$directory_path" ]; then
    echo "Directory does not exist: $directory_path"
    exit 1
fi

# List all files in the directory, filter for unique names before "_1.fastq.gz" or "_2.fastq.gz"
file_names=$(find "$directory_path" -type f -name '*_1.fastq.gz' -o -name '*_2.fastq.gz' | sed -E 's|.*/([^/]+)_[12].fastq.gz|\1|' | sort -u)


# Loop through unique file names and append to mapfile.txt
for file in $file_names; do
    echo -e "$file\t$directory_path\t$sample_name" >> mapfile.txt
done
