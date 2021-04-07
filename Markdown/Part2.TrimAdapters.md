# Trim adapters using Trimmomatic

I will attempt to trim adapters using [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

Here are some highlights from the user manual:

Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop Illumina (FASTQ) data as well as to remove adapters.Trimmomatic works with FASTQ files (using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used). 

## Running Trimmomatic

### Processing Order
The different processing steps occur in the order in which the steps are specified on the command line. It is recommended in most cases that adapter clipping, if required, is done as early as possible.

### Paired End Mode Input and Output files
For paired-end data, two input files, and 4 output files are specified, 2 for the 'paired' output where both reads survived the processing, and 2 for corresponding 'unpaired' output where a read survived, but the partner read did not.

#### Syntax: 
See [user manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for explanation of terms
``` 
PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] #[-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>]
```

### Trimming commands:

  - ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
  - MINLEN: Drop the read if it is below a specified length


#### Syntax:
See user manual for explanation of terms
```
ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
```

#### File with Nextera adapters, as specified in Physalia course:
Saved as "NexteraPE_NT.fa"
```
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>Prefix_PCR/1
AATGATACGGCGACCACCGAGATCTACAC
>Prefix_PCR/2
CAAGCAGAAGACGGCATACGAGAT
>PCR_i5
AATGATACGGCGACCACCGAGATCTACAC
>PCR_i5_rc
GTGTAGATCTCGGTGGTCGCCGTATCATT
>PCR_i7
CAAGCAGAAGACGGCATACGAGAT
>PCR_i7_rc
ATCTCGTATGCCGTCTTCTGCTTG
```

##### My code in progress (have not run yet)

``` bash
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
	base=$(basename --suffix=$SUFFIX1 $infile) # this cute piece of code gets the basename of the fastq files
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





```
