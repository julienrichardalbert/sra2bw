#!/bin/bash
#SBATCH --ntasks=8                      # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes
#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL

if [ $# -ne 1 ]; then
        script_name=$(basename $0)
        echo "Usage: $script_name /full/path/to/<fastq_file>.NO.EXTENSION"
        echo "Note: FASTQ file extension MUST end in _R1.fastq.gz and _R2.fastq.gz or .fastq.gz"
        echo "Example: $script_name /data/priscillia_data/fastq/PL01"
        echo "Would align: /data/priscillia_data/fastq/PL01_R1.fastq.gz and /data/priscillia_data/fastq/PL01_R2.fastq.gz"
        exit 1
fi


# ------------------ LOAD FUNCTIONS -------------------- #
source /home/robin/sra2bw/sra2bw_functions.sh

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=10
THREADS=8
THREAD_MEM=4G
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"
CHROM_SIZES="/data/reference_genomes/mm10/mm10.sizes"

# ------------------ CALL FUNCTIONS -------------------- #

# takes as input the full path to the fastq file
setupVariables $1
trimReads 5
run_fastQC
alignBowtie2_cutNtag /data/reference_genomes/mm10/mm10
groomSam
trueStats
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
callPeaks narrow
calculateEnrichment
rename_cleanup


