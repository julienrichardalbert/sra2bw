#!/bin/bash
source /home/robin/sra2bw/epipax.config

TRIMVER=crop36
TRIM_PARAMS="CROP:36"
THREADS=12

for p in *R1.fastq.gz; do
  FILE=${p//_R1.fastq.gz}
  $TRIMMOMATIC PE -threads $THREADS \
  "$FILE"_R1.fastq.gz "$FILE"_R2.fastq.gz \
  "$FILE"_"$TRIMVER"_R1.fastq.gz trash_R1.fastq.gz "$FILE"_"$TRIMVER"_R2.fastq.gz trash_R2.fastq.gz \
  $TRIM_PARAMS
done


