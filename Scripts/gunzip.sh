#!/bin/bash
#SBATCH --job-name=elp_01
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/elpetrou/fastq/AK_herring


## ENVIRONMENT SETUP
DATADIR1=/gscratch/scrubbed/elpetrou/fastq/AK_herring
DATADIR12=/gscratch/scrubbed/elpetrou/fastq/WA_herring

## CODE FOR JOB

for FILE in $DATADIR1'/'*.fastq.gz 
do
	gunzip $FILE
done


for FILE in $DATADIR2'/'*.fastq.gz 
do
	gunzip $FILE
done
