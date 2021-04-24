#!/bin/bash
#SBATCH --job-name=elp_bowtie2_AK
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=6-8:00:00
## Memory per node
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/mmfs1/gscratch/scrubbed/elpetrou/fastq_trimmed/AK

## Job-specific Variables
###############################################
# Specify the path to your scipt and its name
MYSCRIPT=/mmfs1/home/elpetrou/scripts/parallel_bowtie2.sh

# Specify suffix of files to analyze
SUFFIX1=_R1_001.trim.fastq

################################################
# Make your script executable
chmod +x $MYSCRIPT

# Run the script on each file in the current directory
for i in *$SUFFIX1
do
	echo $i
	$MYSCRIPT $i
done

