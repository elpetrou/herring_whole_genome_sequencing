#!/bin/bash
#SBATCH --job-name=elp_09
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=18:00:00
## Memory per node
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/fastq/AK_herring
SAMPLELIST=trimmomatic_list.txt
SUFFIX1=_R1_001.fastq # Suffix to raw fastq files. The forward reads with paired-end data.
SUFFIX2=_R2_001.fastq # Suffix to raw fastq files. The reverse reads with paired-end data. 
ADAPTERFILE=/mmfs1/home/elpetrou/scripts/NexteraPE_EP.fa # File with adapter sequences. This is based on the file that was provided in the Physalia course, and double-checked using the Illumina Adapter Sequences pdf
OUTDIR=/mmfs1/gscratch/scrubbed/elpetrou/fastq_trimmed #where to store output files


##############################################################################
## Save the sample names of the forward reads into a text file (for looping thru samples later)
mkdir $OUTDIR

cd $DATADIR
ls *$SUFFIX1 > $SAMPLELIST

##############################################################################
## Run Trimmomatic. For an explanation of trimmomatic commands, see: 
## http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
## https://datacarpentry.org/wrangling-genomics/03-trimming/index.html


for infile in `cat $SAMPLELIST`
do
        base=`basename --suffix=$SUFFIX1 $infile`
        trimmomatic PE \
        -threads 10 -phred33 \
        ${infile} \
        ${base}$SUFFIX2 \
        ${base}_R1_001.trim.fastq \
        ${base}_R1_001.unpaired.trim.fastq \
        ${base}_R2_001.trim.fastq \
        ${base}_R2_001.unpaired.trim.fastq \
        ILLUMINACLIP:$ADAPTERFILE:2:30:10:1:true \
        SLIDINGWINDOW:4:20 \
        MINLEN:40
done

#############################################################################
## Move the results files to the output directory

mv *_R1_001.trim.fastq $OUTDIR
mv *_R2_001.trim.fastq $OUTDIR
mv *_R1_001.unpaired.trim.fastq $OUTDIR
mv *_R2_001.unpaired.trim.fastq $OUTDIR

