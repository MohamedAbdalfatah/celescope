# celescope

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
