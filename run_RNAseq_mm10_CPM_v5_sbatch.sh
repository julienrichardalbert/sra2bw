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
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET_and_rRNA.bed"

# ------------------ CALL FUNCTIONS -------------------- #

setupVariables $1
trimReads 5
run_fastQC
run_fastqScreen
alignSTAR /data/reference_genomes/mm10
groomSam
rRNAmetrics_mm10
trueStats
run_preseq
run_stringTie /data/reference_genomes/mm10/mm10_refseq_genes.gtf
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
rename_cleanup

