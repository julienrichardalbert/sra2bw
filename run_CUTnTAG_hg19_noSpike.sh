#!/bin/bash

# usage
# run_CUTnTAG_mm10_noSpike.sh <full path to file>.NO.EXTENSION
# example:
# ./run_CUTnTAG_mm10_noSpike.sh /media/teamgreenberg/Stagiaire/2021-Line-CUT-TAG/201230_X591_FCHF2YLCCX2_L4_CF1-FLAG
# extension in these cases MUST BE "_R1.fastq.gz"

# load functions coded in separate files
source ./sra2bw_functions.sh

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=10
THREADS=8
THREAD_MEM=4G
BLACKLIST="/data/reference_genomes/hg19/hg19_blacklist_Anshul_Duke_combined.bed"



# ------------------ CALL FUNCTIONS -------------------- #

# takes as input the full path to the fastq file
setupVariables $1
trimReads 5
run_fastQC
run_fastqScreen
alignBowtie2_cutNtag_shortreads /data/reference_genomes/hg19/hg19
groomSam
trueStats
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
rename_cleanup


