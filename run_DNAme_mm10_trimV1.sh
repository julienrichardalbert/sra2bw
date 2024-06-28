#!/bin/bash
#SBATCH --ntasks=4                      # number of MPI processes
#SBATCH --mem-per-cpu=3G               # memory; default unit is megabytes
#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL

# usage
# run_DNAme_mm10.sh <full path to file>.NO.EXTENSION
# example:
# extension MUST BE "_R1.fastq.gz"

# NOTE: USES 80% OF AVAILABLE RAM FOR METHYLATION EXTRACTOR STEP

# load functions coded in separate files
source /home/robin/sra2bw/sra2bw_functions.sh

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=1
THREADS=4
THREAD_MEM=3G
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"
# ------------------ CALL FUNCTIONS -------------------- #

# takes as input the full path to the fastq file
setupVariables $1
trimReads 1
run_fastQC
# run_fastqScreen_convertedDNA

alignBismark /data/reference_genomes/mm10/
collectBismarkStats
mergeTwoStrandMethylation
convertMethylationToBigWig 1 /data/reference_genomes/mm10/mm10.sizes
convertMethylationToBigWig 5 /data/reference_genomes/mm10/mm10.sizes

cleanupBismark
