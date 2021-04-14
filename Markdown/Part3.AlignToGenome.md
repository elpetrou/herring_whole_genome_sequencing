# Align the sequences to Atlantic Herring genome using Bowtie2

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#adding-to-path) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.

## Download latest version of Atlantic herring genome from NCBI

I navigated to the NCBI RefSeq website and found the annotated Atlantic herring genome: https://www.ncbi.nlm.nih.gov/genome/15477?genome_assembly_id=495882
- Assembly: Ch_v2.0.2
- RefSeq: GCF_900700415.1

I downloaded all files (fasta = GCF_900700415.1_Ch_v2.0.2_genomic.fna, gff = GCF_900700415.1_Ch_v2.0.2_genomic.gff) and uploaded them to Klone here: 

``` /gscratch/merlab/genomes/atlantic_herring ```

 
## Install bowtie2 on Klone

I installed bowtie2 on the MerLab node by logging into a Klone terminal and typing these commands:

```
cd /gscratch/merlab/software/miniconda3/bin
conda install -c bioconda bowtie2
```

## Index the reference genome with bowtie2
The first and basic step of running Bowtie2 is to build Bowtie2 index from a reference genome sequence. The basic usage of the command bowtie2-build is:

```
bowtie2-build -f input_reference.fasta index_prefix

```
input_reference.fasta is an input file of sequence reads in fasta format, and index_prefix is the prefix of the generated index files. The option -f that is used when the reference input file is a fasta file.

Here is how I indexed the genome:

```
srun -p compute-hugemem -A merlab --nodes=1 \
--ntasks-per-node=1 --time=02:00:00 \
--mem=40G --pty /bin/bash

# Specify files
GENOME=/gscratch/merlab/genomes/atlantic_herring/GCF_900700415.1_Ch_v2.0.2_genomic.fna
INDEX_PREFIX=GCF_900700415.1_Ch_v2.0.2

# Run bowtie2-build
bowtie2-build -f $GENOME $INDEX_PREFIX

```
The files ending in ".bt2" are the indexed genome.


## Align sequences to reference genome with bowtie2

The command bowtie2 takes a Bowtie2 index and set of sequencing read files and outputs set of alignments in SAM format. The script below is still a work in progress - I haven't yet figured out how to parallelize this process in the most efficient way. 

```
#!/bin/bash
#SBATCH --job-name=elp_bowtie2WA
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ##########
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/fastq_trimmed/WA #directory containing the trimmed fastq files
GENOMEDIR=/gscratch/merlab/genomes/atlantic_herring #directory containing the genome
GENOME_PREFIX=GCF_900700415.1_Ch_v2.0.2 #prefix of .bt2 files made by bowtie2
SUFFIX1=_R1_001.trim.fastq # Suffix to trimmed fastq files. The forward reads with paired-end data.
SUFFIX2=_R2_001.trim.fastq # Suffix to trimmed fastq files. The reverse reads with paired-end data.
OUTDIR=/mmfs1/gscratch/scrubbed/elpetrou/sam #where to store output (sam) files


############################################################################
cd $DATADIR

## I am trying to use xargs to parallelize this task. Some notes:
## find - search for files in a directory hierarchy
##  xargs - build and execute command lines from standard input
## basename --suffix=SUFFIX: remove a trailing SUFFIX
## xargs -I: Replaces occurrences of replace-str in the initial-arguments with names read from standard input


find *$SUFFIX1 | xargs basename --suffix=$SUFFIX1 | xargs -I{} bowtie2 \
-x $GENOMEDIR'/'$GENOME_PREFIX \
--phred33 -q \
-1 {}$SUFFIX1 \
-2 {}$SUFFIX2 \
-S {}.sam \
--very-sensitive \
--minins 0 --maxins 1500 --fr \
--threads 30 \
--rg-id {} --rg SM:{} --rg LB:{} --rg PU:Lane1 --rg PL:ILLUMINA



#############################################################################
## Move the results files to the output directory

mv *sam $OUTDIR
    
```



