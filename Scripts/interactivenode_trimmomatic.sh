##############################################################################################################
# Trial on 21210406

# Use interactive node to check that singularity is working
srun -p compute-hugemem -A merlab --nodes=1 \
--ntasks-per-node=1 --time=02:00:00 \
--mem=100G --pty /bin/bash


# Specify the paths to input files and directories
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/fastq/WA_herring
SAMPLELIST=trimmomatic_list.txt
SUFFIX1=_R1_001.fastq # Suffix to raw fastq files. The forward reads with paired-end data.
SUFFIX2=_R2_001.fastq # Suffix to raw fastq files. The reverse reads with paired-end data. 
ADAPTERFILE=NexteraPE_NT.fa # File with adapter sequences. This is the file that Nina provided in the physalia course. Should check
OUTDIR=/mmfs1/gscratch/scrubbed/elpetrou/fastq_trimmed #where to store output files


##############################################################################

# Load the singularity module
module load singularity
MYSINGULARITY=/mmfs1/gscratch/merlab/singularity_sif/BioinfoContainers_trimmomatic.sif 

# Check software version
singularity exec \
$MYSINGULARITY \
trimmomatic --version

##############################################################################
# Save the sample names of the forward reads into a text file (for looping thru samples later)
mkdir $OUTDIR

cd $DATADIR
ls *$SUFFIX1 > $SAMPLELIST  

##############################################################################
# Run Trimmomatic. For an explanation of trimmomatic commands, see: 
# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
# https://datacarpentry.org/wrangling-genomics/03-trimming/index.html

# Testing this piece of code
for infile in `cat $SAMPLELIST` 
do
	base=$(basename --suffix=$SUFFIX1 $SAMPLE) # this cute piece of code gets the basename of the fastq files
	echo $base \
	echo ${base}$SUFFIX2 \
	singularity exec \
	$MYSINGULARITY \
	trimmomatic PE \
	-threads 1 -phred33 \
	${infile} \
	${base}$SUFFIX2 \
	${base}_R1_001.trim.fastq \
	${base}_R1_001.unpaired.trim.fastq \
	${base}_R2_001.trim.fastq \
	${base}_R2_001.unpaired.trim.fastq \
	ILLUMINACLIP:$ADAPTERFILE:2:30:10:1:true \
	MINLEN:40 
done


#############################################################################
# Move the results files to the output directory

# Test this to make sure it works (I think it will)
mv *_R1_001.trim.fastq $OUTDIR
mv *_R2_001.trim.fastq $OUTDIR
mv *_R1_001.unpaired.trim.fastq $OUTDIR
mv *_R2_001.unpaired.trim.fastq $OUTDIR

