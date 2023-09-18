# Celescope

in this repository we are going to set up and analyze data usng celescope from company 
we are going to follow the instactions here https://github.com/singleron-RD/CeleScope/blob/master/doc/user_guide.md and we are going to analyze the subpeoject SCGTEST_50

## Installation

### Clone the repository 
```{}
git clone https://github.com/singleron-RD/CeleScope.git
cd CeleScope
```
### Create conda enviroment 
```{}
source ~/.bashrc
conda create -n celescope
conda activate celescope
```
### Install required packages
Since we are now in CeleScope directory, we can find a file called conda_pkgs.txt
```{}
cat conda_pkgs.txt
```

This file include conda packages we need to install, we should run:

```{}
conda create -n celescope -y --file conda_pkgs.txt
```
using above command all packages there should be insralled but this is not working as excpected, so we can install all of them manually:

```
conda install -c bioconda star=2.6.1b
conda install -c bioconda subread=2.0.1
# This samtools version seems to be not avilable in conda 
#conda install -c bioconda samtools=1.16.1
conda install -c bioconda samtools
conda install -c bioconda bioconda::igblast
# This bcftools seems to be not avilable in conda 
#conda install -c bioconda bcftools=1.16 
conda install -c bioconda bcftools
conda install -c bioconda gatk4
conda install -c bioconda trust4=1.0.7
conda install -c bioconda snpeff
# We shouldn't install -c conda-forge gcc, it is alredy installed in lunix, I instaled and removed it, becouse the conda version of gcc is dfferent then the gcc++ version 
#conda install -c conda-forge gcc
```
### Install celescope package using pip
```
~/anaconda3/envs/celescope/bin/pip install celescope
```
Now we have everything we need to analyze our data using celescope, since our data we are analyze is Single cell GEX we are going to follow the instructions from this toutorial https://github.com/singleron-RD/CeleScope/blob/master/doc/assay/multi_rna.md

## Create A Refernce Genome 

### Install the Refernce FASTA and GTF
Since our data consists of human genetic information, we need to install the human reference genome in order to map our FASTQ reads accurately. Additionally, we will install the GTF file for the reference genome, which contains the coordinates of genes and chromosomes.
Note: This step may take some minutes
```{}
mkdir hs_ensembl_99
cd hs_ensembl_99

wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz
```

### Filter GTF File
This step not working

```{}
celescope utils mkgtf Homo_sapiens.GRCh38.99.gtf Homo_sapiens.GRCh38.99.filtered.gtf
```
### Create the Refernce Gemnome
```{}
celescope rna mkref \
 --genome_name Homo_sapiens_ensembl_99_filtered \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh38.99.filtered.gtf
```
**WE WILL RUN ALL SCRIPTS IMSIDE THIS DIRECTORY
```{}
mkdir /projects/SCGTEST_50
```

# Get LIMS info

### The script 
```
#!/bin/bash

# Get information for each library (flow cell, lane, sample id, etc.)
# $1  needs to be the name of the project
/home/groups/singlecell/mabdalfttah/projects/scripts/limsq.py -sp $1 | sed 's/;/                                                                                                                                                             \t/g' > "lims_info_"$1".txt"

echo "Created LIMS information file: lims_info.txt"
```
### How to run

```
./1-lims.sh SCGTEST_50
```
Now we have lims_info_SCGTEST_50.txt with all of information of the samples and project

since some fastq files/samples doesn't pass the filters we need to remove them from the lims_info_SCGTEST_50.txt file 

```
awk -F'\t' -v column="LanePassFail" 'BEGIN {OFS=FS} NR==1 {for (i=1; i<=NF; i++) if ($i == column) col=i} $col != "fail"' lims_info_SCGTEST_50.txt > tmp_file && mv tmp_file lims_info_SCGTEST_50.txt
```

We have filterd lims file, we are ready for the next step 

# Get FASTQs Path

### The script 

```
#!/usr/bin/env python

# Writes fastq path by arranging proper flowcell, lane, index and read for a set of libraries

# Load packages
import numpy as np
import pandas as pd
import os
import argparse


# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to transfer feature-barcodes matrices from cluster to lcal")
parser.add_argument("--subproject",
                    dest = "subproject",
                    action = "store",
                    default = None,
                    help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--info_file",
                    dest = "info_file",
                    action = "store",
                    default = None,
                    help = "Tab-delimited file with the information of Illumina sequence of libraries for that subproject")


options = parser.parse_args()
subproject = options.subproject
info_file = options.info_file

# Read file
lims = pd.read_csv(info_file, sep = "\t", header = 0)

# Assemble fastq paths combining flowcell, lane and index
fastq_path = "/scratch/project/production/fastq"
fastq_path_list_r1 = []
fastq_path_list_r2 = []
for idx in lims.index:
    fc = lims.loc[idx, "flowcell"]
    lane = lims.loc[idx, "lane"]
    index = lims.loc[idx, "index"]
    fastq_path_r1 = "{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
    fastq_path_r2 = "{}/{}/{}/fastq/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
    fastq_path_list_r1.append(fastq_path_r1)
    fastq_path_list_r2.append(fastq_path_r2)
library_id_l = list(lims["id"].append(lims["id"]))
p_l = "P" * len(fastq_path_list_r1)
indx_l = list(range(1, len(fastq_path_list_r1) + 1))
pair_id = [p_l[x] + str(indx_l[x]) for x in range(len(indx_l))]
fastq_path_list_r1.extend(fastq_path_list_r2)
pair_id.extend(pair_id)
fastq_path_l = fastq_path_list_r1
read_l = (["R1"] * lims.shape[0]) + (["R2"] * lims.shape[0])
fastq_dict = {"library_id":library_id_l, "fastq_path":fastq_path_l, "read":read_l, "pair_id":pair_id}
fastq_df = pd.DataFrame(fastq_dict)

fastq_df.to_csv("fastq_paths.tab".format(subproject), header = True, index = False, sep="\t")
```
### How to run 
First we need to activate any conda env with python:
```
source ~/.bashrc
conda activate sc_py
```
Run the script 
```
python 2-write_fastq_paths.py --subproject SCGTEST_49 --info_file lims_info_SCGTEST_49.txt
```

# 5- Create a Metadata File
This step should be in R 

```
Path = "../Downloads/"
library(tidyverse)
#===============================================================================
Files = list.files(paste0(Path), pattern = "lims_info_SCGTEST_49.txt")
All_Files = list()
metadata = list()
for (i in seq_along(Files)) {
  All_Files[[i]] = read.table(paste0(Path, Files[i]), sep = "\t", header = T)
  metadata[[i]] = data.frame(subproject = All_Files[[i]]$subproject, gem_id = All_Files[[i]]$SampleName,
                             library_id = All_Files[[i]]$id, library = All_Files[[i]]$library,
                             type = "not_hashed",donor_id = All_Files[[i]]$SampleName, flowcell = All_Files[[i]]$flowcell,
                             lane = All_Files[[i]]$lane, index = All_Files[[i]]$index)
  metadata[[i]]$gem_id = str_replace_all(string = metadata[[i]]$gem_id, pattern = "\\.", replacement = "_")
  
}
write.csv(metadata[[1]],paste0("../Downloads/SCGTEST_49.csv"), row.names = F)
```
Now we have SCGTEST_50.csv metadata file, let's go for the next step

# 6- Create a jobs directories and copy FASTQs to them

In this script we create a directory for each samples and copy the fASTQs files to this directory 

### The script 
```
# This script initializes the filesystem of this project:
# It creates a "jobs" folder which contains as many subdirectories as samples it has
# For each sample directory, it creates the following files/folders:
# 1. fastq: dir with the symlinks pointing to the fastq files
# 2. log: dir which contains standard error and output of cellranger
# 3. (sample_id).cmd: job script to compute the features-barcode matrix using cellranger


# Import required packages
import numpy as np
import pandas as pd
import os
import argparse
import subprocess
import re
import sys
import config_vars as cfg
from utils import *


# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to initialize the filesystem and scripts of this project")
parser.add_argument("--subproject",
                    dest = "subproject",
                    action = "store",
                    default = None,
                    help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--gem_id",
                    dest = "gem_id",
                    action = "store",
                    default = None,
                    help = "Gel Beads in Emulsion id")
parser.add_argument("--verbose",
                    dest = "verbose",
                    action = "store_true",
                    default = False,
                    help = "Print log in standard error")
parser.add_argument("--metadata",
                    dest = "metadata",
                    action = "store",
                    default = None,
                    help = "Metadata csv file for the tonsil atlas project")
parser.add_argument("--fastq_paths",
                    dest = "fastq_paths",
                    action = "store",
                    default = None,
                    help = "File that contains the paths of the fastqs for the subproject libraries")


def create_fastq_symlink_nh(gem_id, fastq_path_df, symlink_path):
    """Creates a symbolic link pointing to a fastq file using cellranger notation

    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      fastq_path_df: pandas dataframe with the fastq paths for that gem_id
      symlink_path: string specifying where to create the symlinks

    Returns:
      None
    """
    pair_ids = np.unique(fastq_path_df["pair_id"])
    for i in range(len(pair_ids)):
        filt = (fastq_path_df["pair_id"] == pair_ids[i])
        pair_df = fastq_path_df.loc[filt, :]
        for j in pair_df.index:
            fastq_path = pair_df.loc[j, "fastq_path"]
            lane = str(i + 1)
            read = pair_df.loc[j, "read"]
            read = read.replace("R", "")
            subprocess.run(["ln", "-s", fastq_path, "{}/{}_S1_L00{}_R{}_001.fastq.gz".format(symlink_path, gem_id, lane, read)])

options = parser.parse_args()
subproject = options.subproject
gem_id = options.gem_id
metadata_path = options.metadata
fastq_paths = options.fastq_paths


# Read data
project_dir = "/home/groups/singlecell/mabdalfttah/projects/{}".format(subproject)
fastq_path_df = pd.read_csv(fastq_paths, sep = "\t", header = 0)
metadata_df = pd.read_csv(metadata_path, sep = ",", header = 0)
if options.verbose:
    sys.stderr.write("Files read successfully!\n")


# For each sample, create directories and jobscript
if not os.path.exists("{}/jobs".format(project_dir)):
    os.mkdir("{}/jobs".format(project_dir))
filt = (metadata_df["gem_id"] == gem_id)
metadata_df = metadata_df.loc[filt]


# Create directories
subproject_dir = "{}/jobs/{}".format(project_dir, gem_id)
fastq_dir = "{}/fastq".format(subproject_dir)
log_dir = "{}/log".format(subproject_dir)
for direct in [subproject_dir, fastq_dir, log_dir]:
    if not os.path.exists(direct):
        os.mkdir(direct)


# Define variables and subset dataframes
library_id = metadata_df.loc[filt, "library_id"]
fastq_sub_df = fastq_path_df.loc[fastq_path_df["library_id"].isin(library_id), :]
type = metadata_df["type"]
type = type.values[0]

# Create symmlinks to fastq files
create_fastq_symlink_nh(gem_id, fastq_sub_df, fastq_dir)
```

### How to run
```
python 3-copy_fastqs  --subproject SCGTEST_50 --fastq_paths fastq_paths.tab --metadata SCGTEST_50.csv --gem_id CNAG_61_1
python 3-copy_fastqs  --subproject SCGTEST_50 --fastq_paths fastq_paths.tab --metadata SCGTEST_50.csv --gem_id CNAG_61_2
```


# Create a Mapfile
Mapfile is a Required tab-delimited text file with at least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.
2nd column: Fastq file directory path.
3rd column: Sample name, which is the prefix of all output files.
4th column: The 4th column has different meaning for each assay. For rna, it means forced cell number and it's an optional column. For other assays, see here.

To generate this file we need to create a bash script:

### The script
```{}
#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <Prefix> <Directory> <Sample>"
    exit 1
fi

prefix="$1"
directory="$2"
sample="$3"

# Generate the output
output="${prefix}       ${directory}    ${sample}"

# Append the output to the mapfile.txt file
echo -e "$output" >> mapfile.txt
```
This file take mainly three inputs: 1- Fastq file prefix 2- Fastq file directory path 3- Sample name
### How to use
```{}
./4-generate_mapfile.sh CNAG_61_1 /home/groups/singlecell/mabdalfttah/projects/SCGTEST_50/jobs/CNAG_61_1/fastq CNAG_61_1
./4-generate_mapfile.sh CNAG_61_2 /home/groups/singlecell/mabdalfttah/projects/SCGTEST_50/jobs/CNAG_61_2/fastq CNAG_61_2
```

This should generate a file looks like:
```{}
CNAG_61_1	/home/groups/singlecell/mabdalfttah/projects/SCGTEST_50/jobs/CNAG_61_1/fastq	CNAG_61_1
CNAG_61_1	/home/groups/singlecell/mabdalfttah/projects/SCGTEST_50/jobs/CNAG_61_1/fastq	CNAG_61_1
```

# Generate scripts for each sample

Since we have a different samples, Celescope have a script to generate a script for each sample seperetly, we need to define where is the mapfile and where is refernce genome

### The Script

```{}
#!/bin/bash
multi_rna\
        --mapfile ./mapfile.txt\
        --genomeDir /home/groups/singlecell/mabdalfttah/CeleScope/hs_ensembl_99\
        --thread 8\
        --mod shell
```

# How To use

```{}
sh 5-run.sh
```
After you sh run.sh, a **shell** directory containing {sample}.sh files will be generated, in our case we have shell/CNAG_61_1.sh shell/CNAG_61_2.sh

# Create a Cluster job
Since this script should run directory but it doesn't assum you are going to run this in a cluster, we need generate a .cmd file to make this file run in the cluster 

### The Script
```{}
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
```

# How To Run
```{}
chmod +x 6-add_slurm_info.sh
./6-add_slurm_info.sh shell/CNAG_61_1.sh
./6-add_slurm_info.sh shell/CNAG_61_2.sh
```
You will find a two script called shell/CNAG_61_1.cmd shell/CNAG_61_2.cmd inside shell directory 

