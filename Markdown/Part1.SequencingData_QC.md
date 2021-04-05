
# WGS data

## Data download and backup

Sequencing data were downloaded from the NW Genomics Center using Globus software. The raw data (.tar files) have been saved to three locations:
  1. An external hard drive in my office
  2. The HauserLab fireproof external hard drive
  3. The /gscratch/scrubbed/elpetrou directory on Klone

*NB: We should also back up the data on LOLO archive

For each of these data backups, I used the MD5sum file to verify that the data were not corrupted. 

## Plot distribution of raw sequences per sample

I plotted the distribution of raw sequencing reads for each herring sample [using this R script](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/6131687bd581f88f8a08c3e23cd06f516001c82b/Scripts/plot_distro_raw_seqs.R).
These data were provided to me by the NW Genomics Sequencing Center in .csv format. 

![raw seq distro](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/11a515129c73adc8c18a78f0db3a0f224e851bee/Markdown/raw_seq_distro.jpeg) 

Here is a summary of the raw sequencing data :
  - We sequenced 556 herring (number of WA samples = 281; number of AK samples = 275)
  - Average number of reads per sample = 12.88 million reads
  - Range = 0 to 80.06 million reads
  - Standard deviation = 4.41 million reads
  - Seven samples sequenced poorly (reads per sample less than 2 sd from the mean)
  - Nine samples sequenced very deeply (reads per sample more than 2 sd from mean)

## Unzip data

I unzipped the data on Klone [using this script](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/102c8f2fcdacae63b32a074be61be1d13fdb52a1/Scripts/gunzip.sh)

## Run FastQC

To check the quality of the raw sequence data I ran the software FastQC [using this script](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/102c8f2fcdacae63b32a074be61be1d13fdb52a1/Scripts/fastqc.sh) on Klone

## Visualize FastQC output using MultiQC

I visualized the voluminous FastQC output using MultiQC software on Klone [using this script](https://github.com/EleniLPetrou/herring_whole_genome_sequencing/blob/102c8f2fcdacae63b32a074be61be1d13fdb52a1/Scripts/multiqc.sh)

