#!/bin/bash


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
THREADS=12
THREAD_MEM=3G
BUFFER_SIZE_BISEXTRACT="40%"
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"
CHROM_SIZES="/data/reference_genomes/mm10/mm10.sizes"
# ------------------ CALL FUNCTIONS -------------------- #

# takes as input the full path to the fastq file
setupVariables $1
trimReads test
run_fastQC
# run_fastqScreen_convertedDNA
alignBismark /data/reference_genomes/mm10/
collectBismarkStats
mergeTwoStrandMethylation
convertMethylationToBigWig 1
convertMethylationToBigWig 5
convertMethylationToBigWig 10

callPeaks narrow
calculateEnrichment
cleanupBismark
