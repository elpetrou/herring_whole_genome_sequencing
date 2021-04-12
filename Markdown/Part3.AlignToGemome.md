# Align the sequences to Atlantic Herring genome using Bowtie 2

## Download latest version of Atlantic herring genome from NCBI
 Notes here
 
## Install bowtie2 on Klone

I installed bowtie2 on the MerLab node by logging into a Klone terminal and typing these commands:

```
cd /gscratch/merlab/software/miniconda3/bin
conda install -c bioconda bowtie2
```

## Some notes about Bowtie2

Bowtie2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.

### Index the reference genome with bowtie2
The first and basic step of running Bowtie2 is to build Bowtie2 index from a reference genome sequence. The basic usage of the command bowtie2-build is:

```
bowtie2-build -f input_reference.fasta index_prefix

```
input_reference.fasta is an input file of sequence reads in fasta format, and index_prefix is the prefix of the generated index files. The option -f that is used when the reference input file is a fasta file.

### Align sequences to reference genome with bowtie2

The command bowtie2 takes a Bowtie2 index and set of sequencing read files and outputs set of alignments in SAM format. An example of how to run Bowtie2 local alignment with paired-end fastq files and 8 CPUs is shown below

```
#!/bin/bash
#SBATCH --job-name=elp_bowtie2
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=2-12:00:00
## Memory per node
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu

##### ENVIRONMENT SETUP ##########
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/fastq_trimmed

bowtie2 -x index_prefix -q \
-1 input_reads_pair_1.fastq \
-2 input_reads_pair_2.fastq \
-S bowtie2_alignments.sam \
--sensitive \
-p $SLURM_NTASKS_PER_NODE



```



