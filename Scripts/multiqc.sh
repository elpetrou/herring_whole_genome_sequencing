# As an srun script to check that fastQC singularity is working
srun -p compute-hugemem -A merlab --nodes=1 \
--ntasks-per-node=1 --time=02:00:00 \
--mem=100G --pty /bin/bash


# Specify the path and the name of the singularity you want to use
DATADIR=/mmfs1/gscratch/scrubbed/elpetrou/fastqc

cd $DATADIR


# Load the singularity module
module load singularity
MYSINGULARITY=/mmfs1/gscratch/merlab/singularity_sif/multiqc-srf_1.9.sif

# Use the singularity exec command to use the singularity and run commands that are specific to the software

singularity exec \
$MYSINGULARITY \
multiqc . 

