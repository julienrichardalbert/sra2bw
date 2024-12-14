#!/bin/bash

# check to make sure the correct number of parameters are given
if [ $# -ne 1 ]; then
    script_name=$(basename $0)
    echo "Usage: $script_name genomeBuild"
    echo "E.g.: mm10, hg19 or sacCer3"
    exit 1
fi


source /home/robin/sra2bw/sra2bw_functions.sh # load standard HTS functions
DEPENDENCIES=($TWOBIT2FA $BISMARK_INDEX $SAMTOOLS $BOWTIE2)
SCRATCH="/scratch"
#OUTPUT_BUILD_FOLDER
checkDependencies # a function in sra2bw_functions.sh
download_genome $1

