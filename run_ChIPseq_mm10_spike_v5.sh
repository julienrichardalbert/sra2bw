#!/bin/bash
#SBATCH --ntasks=8                      # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes
#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL

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
MIN_MAPQ_SPIKE=10
THREADS=12
THREAD_MEM=5G
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"

# ------------------ CALL FUNCTIONS -------------------- #
# takes as input the full path to the fastq file
setupVariables $1
trimReads 5
run_fastQC
alignBowtie2 /data/reference_genomes/dm6/dm6
calculateSpikeFactor #will make global variable $SPIKE_SCALEFACTOR and remove sam file
alignBowtie2 /data/reference_genomes/mm10/mm10
groomSam
trueStats
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
bamToBigWigDeeptoolsSpike 1 0
callPeaks broad
calculateEnrichment
rename_cleanup
