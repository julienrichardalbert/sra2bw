#!/bin/bash
#SBATCH --ntasks=2                      # number of MPI processes
#SBATCH --mem-per-cpu=2G               # memory; default unit is megabytes

# load functions coded in separate files
source ./sra2bw_functions.sh
SCRATCH="/scratch"

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=10
THREADS=2
THREAD_MEM=2G
BLACKLIST="/data/reference_genomes/hg19/hg19_blacklist_Anshul_Duke_combined.bed"

# ------------------ CALL FUNCTIONS -------------------- #
# takes as input the full path to the fastq file
setupVariables $1
trimReads 5
run_fastQC
alignBowtie2 /data/reference_genomes/hg19/hg19
groomSam
trueStats
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
rename_cleanup
