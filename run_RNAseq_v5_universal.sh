#!/bin/bash

#SBATCH --ntasks=8                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

# load functions coded in separate files
# added rRNA sequences to the blacklist file
# this will ignore ribosomal reads when deeptools calculates scaling factor
source /home/robin/sra2bw/sra2bw_functions.sh
SCRATCH="/scratch2"

# ---------------- CONFIGURE VARIABLES ----------------- #

FLAG=1540
MIN_MAPQ=255
THREADS=8
THREAD_MEM=4G
BLACKLIST="/data/reference_genomes/adapters/empty.bed"

# ------------------ CALL FUNCTIONS -------------------- #

setupVariables $1
trimReads 5
run_fastQC
run_fastqScreen
alignSTAR $2
groomSam
rRNAmetrics_mm10
trueStats
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
rename_cleanup

