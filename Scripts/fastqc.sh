#!/bin/bash
#SBATCH --job-name=elp_06
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=6:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ##########
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/fastq/AK_herring
OUTDIR=/mmfs1/gscratch/scrubbed/elpetrou/fastqc
SAMPLELIST=fastq_list.txt
MYSINGULARITY=/gscratch/merlab/singularity_sif/singularity-fastqc_0.11.9.sif

## Load modules
module load singularity

#### CODE FOR JOB #####

## go into the data directory
cd $DATADIR 

## save list of fastq files in this directory to a text file (for looping later)
ls *.fastq > fastq_list.txt 

## Use the singularity exec command to run fastqc program. The for loop will do this for all files in this directory

for SAMPLEFILE in `cat $SAMPLELIST`
do
	singularity exec \
	$MYSINGULARITY \
	fastqc -f fastq --extract \
	$SAMPLEFILE
done

## Move all the fastqc results files (ending in _fastqc, .zip, .html) to the output directory

for FILE in $DATADIR'/'*.html 
do
	mv $FILE $OUTDIR
done



for FILE in $DATADIR'/'*.zip 
do
	mv $FILE $OUTDIR
done



for FILE in $DATADIR'/'*_fastqc
do
	mv $FILE $OUTDIR
done



