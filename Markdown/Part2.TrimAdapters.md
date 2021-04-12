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

### File with Nextera adapters, as specified in Physalia course and Illumina Adapter Sequences document:

I made a [fasta file containing Illumina Nextera adapter sequences](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/6156aca1bec94cb8261570e0636fa7d9a3c236f5/Scripts/NexteraPE_EP.fa). These are the sequences that I will trim away from my raw sequencing data using Trimmomatic. 

### Script to run Trimmomatic over a directory of fastq files on Klone
 
Next, I wrote a [sbatch script to run trimmomatic on Klone](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/4f5ca5f8b5585f7ed5a01cf4e0459be749d0ae75/Scripts/trimmomatic_sbatch.sh)

This script removes Illumina adapters, trims sequences if their phred score drops below 20 (SLIDINGWINDOW:4:20), and removes sequences that are shorter than 40 bp (MINLEN:40). 

I tested and debugged this script on 20210407, and it ran -  hooray! Here is what was printed to the terminal:

``` Input Read Pairs: 6825779 Both Surviving: 5932557 (86.91%) Forward Only Surviving: 468865 (6.87%) Reverse Only Surviving: 215887 (3.16%) Dropped: 208470 (3.05%) ```




