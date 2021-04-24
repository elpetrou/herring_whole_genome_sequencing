#!/bin/bash
##### ENVIRONMENT SETUP ##########
GENOMEDIR=/gscratch/merlab/genomes/atlantic_herring #location of genome
GENOME_PREFIX=GCF_900700415.1_Ch_v2.0.2 #prefix of .bt2 files made by bowtie2
SUFFIX1=_R1_001.trim.fastq # Suffix to trimmed fastq files. The forward reads with paired-end data.
SUFFIX2=_R2_001.trim.fastq # Suffix to trimmed fastq files. The reverse reads with paired-end data.

###################################

# Save the base name of each input file
MYBASE=$(basename --suffix=$SUFFIX1 "$1")

# Sanity check
echo "$1"
echo $MYBASE
echo ${MYBASE}$SUFFIX1
echo ${MYBASE}$SUFFIX2
echo ${MYBASE}.sam

# Run bowtie
bowtie2 -x $GENOMEDIR'/'$GENOME_PREFIX \
--phred33 -q \
-1 ${MYBASE}$SUFFIX1 \
-2 ${MYBASE}$SUFFIX2 \
-S ${MYBASE}.sam \
--very-sensitive \
--minins 0 --maxins 1500 --fr \
--threads ${SLURM_JOB_CPUS_PER_NODE} \
--rg-id ${MYBASE} --rg SM:${MYBASE} --rg LB:${MYBASE} --rg PU:Lane1 --rg PL:ILLUMINA


