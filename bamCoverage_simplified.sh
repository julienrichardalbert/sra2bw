#!/bin/bash

#SBATCH --ntasks=8                      # number of MPI processes
#SBATCH --mem-per-cpu=8G               # memory; default unit is megabytes
#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL

# JRA October 2020

source /home/robin/sra2bw/sra2bw_functions.sh

# Check the number of command line arguments
if [[ $# -ne 4 ]]; then
	script_name=$(basename $0)
	echo "Usage: sbatch $script_name input.bam min_mapq bin_size smooth_length"
	exit 1
fi

BLACKLIST="/data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed"

function bamToBigWigDeeptoolsCPMsmoothKpDup {

	# where $1 = first parameter
	# $2 = second parameter...
	# ./bamCoverage_simplified.sh AM_2i-tKO1-CTCF.bam 10 25 100

	local FILE=$1
	local MIN_MAPQ=$2
	local BIN_SIZE=$3
	local SMOOTH_LEN=$4

	printProgress "[bamToWigCPMsmooth] Started with bin size: $BIN_SIZE and smoothing over $SMOOTH_LEN bp"

#	cd /scratch/
#	ln -s $FILE .

	$BAMCOVERAGE \
	--binSize $BIN_SIZE \
	--smoothLength $SMOOTH_LEN \
	--minMappingQuality $MIN_MAPQ \
	-p 8 \
	--normalizeUsing CPM \
	--outFileFormat bigwig \
	--blackListFileName $BLACKLIST \
	--ignoreForNormalization chrX chrM chrY \
	-b $FILE \
	--outFileName "${FILE//.bam/}"_kpDup_q"$MIN_MAPQ"_b"$BIN_SIZE"_s"$SMOOTH_LEN"_CPM.bw

	printProgress "[bamToWigCPMsmooth] finished successfully"
}

bamToBigWigDeeptoolsCPMsmoothKpDup $1 $2 $3 $4
