#!/bin/bash
#SBATCH --ntasks=12 # number of MPI processes
#SBATCH --mem-per-cpu=4G  # memory; default unit is megabytes

if [ $# -ne 1 ]; then
        script_name=$(basename $0)
        echo "Usage: sbatch $script_name /full/path/to/<fastq_file>.NO.EXTENSION"
	echo "Note: FASTQ file extension MUST end in _R1.fastq.gz and _R2.fastq.gz or .fastq.gz"
        exit 1
fi


# ------------------ LOAD FUNCTIONS -------------------- #
source /home/robin/sra2bw/sra2bw_functions_scaffolds.sh

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=1
THREADS=12
let BISMARK_THREAD=$THREADS/4
THREAD_MEM=4G
BUFFER_SIZE_BISEXTRACT="10%"
BLACKLIST="/data/reference_genomes/adapters/empty.bed"
CHROM_SIZES="/data/reference_genomes/rheMac10/rheMac10.sizes"

# ------------------ CALL FUNCTIONS -------------------- #

# takes as input the full path to the fastq file
setupVariables $1
trimReads 1
run_fastQC
alignBismark /data/reference_genomes/rheMac10/
collectBismarkStats
mergeTwoStrandMethylation
convertMethylationToBigWig 1
convertMethylationToBigWig 5
cleanupBismark
