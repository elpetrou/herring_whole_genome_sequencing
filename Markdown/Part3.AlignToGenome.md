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

The command bowtie2 takes a Bowtie2 index and set of sequencing read files and outputs set of alignments in SAM format. I wrote a bash script containing the bowtie2 commands for aligning paired end data. Here it is:

```
!/bin/bash
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

```

I did some tests and found that it took ~30 min to align a single sample to the genome, using 20 or 32 threads. Additionally, memory is not an issue for bowtie2 because it has a small memory footprint (3-4 Gb). Given that our node on Klone has 40 cores, I realized that I should stack "2 jobs" with 20 cores each and run them in parallel. To do this, I divided the fastq files into two folders and used a separate sbatch script to submit the job for each folder, like this : ` sbatch slurm_parallel_bowtie2.sh`. It's not the most elegant solution in the world but it reduced overall analysis time to 15 min per sample. 

Here is the sbatch script to do this:

```
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

```

On average, it looks like ~86% of sequences are aligning to the Atlantic herring genome and 55% are aligning uniquely. I hope this translates into some nice data downstream!

## Use Samtools to convert sam to bam, format filter for quality & alignment, and sort bame files

I first had to install samtools on Klone, and that was a pain because I was having trouble with the version of openssl that samtools uses. To get around this issue, I decided to create a special conda environment for samtools. This how I did it:

```
# Install samtools using a conda environment (to get around the open ssl bug):
#The syntax for the first command says “conda” runs conda, “create -n” creates a new environment
conda create -n samtools_env

# To activate or enter the environments you just created simply type:
conda activate samtools_env

# Once in the environment, install the versions of htslib samtools and openssl that you want:
conda install -c bioconda htslib samtools openssl=1.0

```
Once samtools was installed and functionining properly, I had to figure out how to call it from within an sbatch script. It turns out you use the command `source activate` rather than `conda activate` (that will cause your script to fail. So arcane, this knowledge!

Finally, I was able to run samtools, using this script (samtools_sbatch.sh):


```
#!/bin/bash
#SBATCH --job-name=elp_samtools
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
## Specify the directory containing data
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/bam #directory with sam files
SUFFIX1=.sam #file suffix
MYCONDA=samtools_env #name of the conda environment that has samtools

## activate the conda environment that runs samtools
source activate $MYCONDA


###################################################################################################################
## Move into the working directory and run script
cd $DATADIR

## Run samtools commands. This takes about 5 min per sample (so like 2 days total for whole data set?)
for MYSAMPLEFILE in *$SUFFIX1
do
    echo $MYSAMPLEFILE
    MYBASE=`basename --suffix=$SUFFIX1 $MYSAMPLEFILE`
    samtools view -bS -F 4 $MYBASE'.sam' > $MYBASE'.bam'
    samtools view -h -q 20 $MYBASE'.bam' | samtools view -buS - | samtools sort -o $MYBASE'_minq20_sorted.bam'
    samtools index $MYBASE'_minq20_sorted.bam'
done

## Flag explanations for samtools view:
## -b       output BAM
## -h       include header in SAM output
## -q INT   only include reads with mapping quality >= INT [0]
##-F INT   only include reads with none of the bits set in INT set in FLAG [0] (aka when this is set to 4, you remove unmapped reads)


```

