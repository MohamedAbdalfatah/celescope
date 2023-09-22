#!/bin/bash

# Check if the user provided an input file as an argument
if [ $# -ne 1 ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

input_file="$1"
output_file="${input_file%.sh}.cmd"  # Generate the output filename with .cmd extension

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file '$input_file' does not exist."
  exit 1
fi

# Define the content to be added to the input file
content_to_add=$(cat <<EOF
#!/bin/bash
#SBATCH --job-name=CNAG_81_GEX1
#SBATCH --mail-type=all
#SBATCH --mail-user=mohamed.abdalfttah@cnag.crg.eu
#SBATCH --output=%x.slurm.%J.out
#SBATCH --error=%x.slurm.%J.err
#SBATCH --time=11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --partition=genD
#SBATCH --mem=32G

echo [\$(date "+%Y-%m-%d %T")] started job on \$HOSTNAME

export TENX_IGNORE_DEPRECATED_OS=1
export HDF5_USE_FILE_LOCKING=FALSE

EOF
)

# Append the content to the input file and create the output file
echo "$content_to_add" > "$output_file"
cat "$input_file" >> "$output_file"

# Make the output file executable (if needed)
chmod +x "$output_file"

echo "Script added to '$input_file' and saved as '$output_file'"
