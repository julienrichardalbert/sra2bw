#!/bin/bash

# ----------------------- USAGE ------------------------ #
# run_ChIPseq_mm10_CPM.sh <full path to file>.NO.EXTENSION
# example:
# ./run_ChIPseq_mm10_CPM.sh /data/fastq/study/dataset_replicate1
# extension in these cases MUST BE "_R1.fastq.gz"


# load functions coded in separate files
source ./sra2bw_functions.sh
SCRATCH="/scratch"

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=10
THREADS=4
THREAD_MEM=5G
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"

# ------------------ CALL FUNCTIONS -------------------- #
# takes as input the full path to the fastq file
setupVariables $1
trimReads 1
# run_fastQC
# run_fastqScreen
alignBowtie2 /data/reference_genomes/mm10/mm10
groomSam
# trueStats
# run_preseq
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
rename_cleanup
