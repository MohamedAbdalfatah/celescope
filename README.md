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
