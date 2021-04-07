# This file documents how I checked the md5sum hashes for the raw sequencing data after downloading it from the UW Genomics Core using Globus. 
# I also document how I transferred these large files to klone.


# Navigate to folder containing raw sequences and MD5 file
cd /media/ubuntu/herring_LCWGS/raw_sequences
ls

# View the MD5 file for the WA herring
cat Eleni_Lane1_WA_herring.md5 

#The hash is: 97c30345cff1381395dfa08adc1ff2b3  Eleni_Lane1_WA_herring.tar

# To cerify that the file is not corrupted, type: md5sum <filename>. The hash should match the hash above.
md5sum Eleni_Lane1_WA_herring.tar

# Here is the output of this command in the terminal:
#97c30345cff1381395dfa08adc1ff2b3  Eleni_Lane1_WA_herring.tar
# Woohoo!! It matches! No corruption here!



# View the MD5 file for the AK herring
cat Eleni_Lane2_AK_herring.md5

#The hash is: c38adac2275099474b828a2fa46b77ce  Eleni_Lane2_AK_herring.tar


# To cerify that the file is not corrupted, type: md5sum <filename>. The hash should match the hash above.
md5sum Eleni_Lane2_AK_herring.tar

# Here is the output of this command in the terminal:
c38adac2275099474b828a2fa46b77ce  Eleni_Lane2_AK_herring.tar
# Woohoo!! It matches! No corruption here!

# Awesome!

################################################
# Transfer the all the sequencing files (.tar files) and md5sum files to /gscratch/scrubbed/elpetrou directory on klone
# using scp. The -r argument of scp works just like the -r arg in cp, it will transfer your entire folder and all the files and subdirectories inside.

# In local machine terminal, type:
cd /media/ubuntu/herring_LCWGS/raw_sequences

DIR1=/media/ubuntu/herring_LCWGS/raw_sequences #where the file lives now
DIR2=elpetrou@klone.hyak.uw.edu:/gscratch/scrubbed/elpetrou #where I want the file to go

scp -r $DIR1 \
$DIR2

# Then I logged into a klone terminal and requested a compute node to run md5sum commands
# This was done using the srun command. The -p argument specifies the partition name, while the -A argument is our group's name (merlab) on Klone.

# Take a peek at what partitions/compute nodes are available (and what their names are) on klone:
sinfo

# request compute node
srun -p compute-hugemem -A merlab --nodes=1 \
--ntasks-per-node=1 --time=01:00:00 \
--mem=20G --pty /bin/bash


md5sum Eleni_Lane1_WA_herring.tar
md5sum Eleni_Lane2_AK_herring.tar

#Output of this command: 97c30345cff1381395dfa08adc1ff2b3  Eleni_Lane1_WA_herring.tar #woohoo! It matches!
#                        c38adac2275099474b828a2fa46b77ce  Eleni_Lane2_AK_herring.tar #woohoo! It matches!
			

