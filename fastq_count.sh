#!/bin/bash
# JRA June 2022

source /home/robin/sra2bw/sra2bw_functions.sh

# Check the number of command line arguments
if [[ $# -ne 0 ]]; then
  script_name=$(basename $0)
  echo "Usage: $script_name"
  exit 1
fi

function count_fastq_reads {
  INPUT_FASTQ=$1
  echo -ne $INPUT_FASTQ"\t" >> ~/fastq_counts.txt #n removes newline after echo. e allows the backslash to be interpreted and the tab to be written
  zcat $INPUT_FASTQ | echo $((`wc -l`/4)) | numfmt --grouping >> ~/fastq_counts.txt #numfmt groups numbers by 3, which looks nice to JRA's eyes
}

# by default, run on all fastq files in working directory
for x in *fastq.gz; do
  count_fastq_reads $x
done
