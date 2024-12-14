#!/bin/bash

#SBATCH --ntasks=12                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

# load functions coded in separate files
# added rRNA sequences to the blacklist file
# this will ignore ribosomal reads when deeptools calculates scaling factor
source /home/robin/sra2bw/sra2bw_functions.sh
SCRATCH="/scratch"

# ---------------- CONFIGURE VARIABLES ----------------- #

FLAG=1540
MIN_MAPQ=255
THREADS=12
THREAD_MEM=4G
# BLACKLIST="/data/reference_genomes/adapters/empty.bed"

# ------------------ CALL FUNCTIONS -------------------- #

setupVariables $1
trimReads 5
run_fastQC
run_fastqScreen
alignSTAR /data/reference_genomes/GRCm39
groomSam
run_stringTie /data/reference_genomes/GRCm39/Mus_musculus.GRCm39.112.gtf
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
rename_cleanup

