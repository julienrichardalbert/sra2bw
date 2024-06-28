#!/bin/bash
#SBATCH --ntasks=6 # number of MPI processes
#SBATCH --mem-per-cpu=4G  # memory; default unit is megabytes


# usage
# /home/robin/run_epimethyltag_mm10.sh <full path to file>.NO.EXTENSION
# extension MUST BE "_R1.fastq.gz"
# e.g.
# sbatch /home/robin/run_epimethyltag_mm10.sh /data/priscillia_data/fastq/PL01


if [ $# -ne 1 ]; then
        script_name=$(basename $0)
        echo "Usage: sbatch $script_name /full/path/to/<fastq_file>.NO.EXTENSION"
	echo "Note: FASTQ file extension MUST end in _R1.fastq.gz and _R2.fastq.gz or .fastq.gz"
	echo "Example: sbatch $script_name /data/priscillia_data/fastq/PL01"
	echo "Would align: /data/priscillia_data/fastq/PL01_R1.fastq.gz and /data/priscillia_data/fastq/PL01_R2.fastq.gz"
        exit 1
fi


# ------------------ LOAD FUNCTIONS -------------------- #
source /home/robin/sra2bw/sra2bw_functions.sh

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=1
THREADS=12
let BISMARK_THREAD=$THREADS/4
THREAD_MEM=3G
BUFFER_SIZE_BISEXTRACT="20%"
BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"
CHROM_SIZES="/data/reference_genomes/mm10/mm10.sizes"

# ------------ SET OUTPUT FOLDER NAMES ----------------- #
SCRATCH="/scratch"
OUTPUT_MAIN="/data"
OUTPUT_BAM_FOLDER="$OUTPUT_MAIN"/bams
OUTPUT_BIGWIG_FOLDER="$OUTPUT_MAIN"/bigWigs
OUTPUT_STATS_FOLDER="$OUTPUT_MAIN"/stats
OUTPUT_CPG_REPORT_FOLDER="$OUTPUT_MAIN"/cpg_reports
OUTPUT_PEAKS_FOLDER="$OUTPUT_MAIN"/peaks


# ------------------ CALL FUNCTIONS -------------------- #

# takes as input the full path to the fastq file
setupVariables $1
trimReads 6
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
