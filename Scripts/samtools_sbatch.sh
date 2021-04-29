#!/bin/bash
#SBATCH --job-name=elp_samtools
#SBATCH --account=merlab
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=4-12:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elpetrou@uw.edu


##### ENVIRONMENT SETUP ##########
# Specify the directory containing data
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/bam #directory with sam files
SUFFIX1=.sam #file suffix
MYCONDA=samtools_env #name of the conda environment that has samtools

# activate the conda environment that runs samtools
source activate $MYCONDA


###################################################################################################################
#Move into the working directory and run script
cd $DATADIR

# Run samtools commands. This takes about 5 min per sample (so like 2 days total for whole data set?)
for MYSAMPLEFILE in *$SUFFIX1
do
    echo $MYSAMPLEFILE
    MYBASE=`basename --suffix=$SUFFIX1 $MYSAMPLEFILE`
    echo $MYBASE
    samtools view -b -h -q 20 -F 4 $MYSAMPLEFILE | \
    samtools sort - | \
    samtools index - $MYBASE'_sorted.bam'
done

# Flag explanations for samtools view:
# -b       output BAM
# -h       include header in SAM output
# -q INT   only include reads with mapping quality >= INT [0]
#-F INT   only include reads with none of the bits set in INT set in FLAG [0] (aka when this is set to 4, you remove unmapped reads)


