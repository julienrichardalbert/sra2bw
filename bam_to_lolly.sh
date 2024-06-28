#!/bin/bash

# Takes bismark-aligned WGBS and makes a bigLolly file for UCSC Browser
# Bismark MapQ for single-matched reads = 44? Trying to filter to avoid blacklist regions
# Maybe discard PCR duplicates at some point?
# Also limiting total pileup to 1000
# Made by Aaron
# Last updated 2023-05-01
# Can pass minDist [20] and size [2] variables into awk script

AWK_SCRIPT="/home/robin/sra2bw/bamlol.awk"
CHR_SIZES="/data/reference_genomes/mm10/mm10.sizes"
BIGLOLLY_AS="/data/reference_genomes/mm10/bigLolly-size.as"
INPUT=$1
OUTPUT=${INPUT//.bam/.bb}

mkdir tempdir
samtools view -q 30 $INPUT | awk -f $AWK_SCRIPT | sort -k1,1 -k2,2n -T ./tempdir > temp.bed
bedToBigBed -as=$BIGLOLLY_AS -type=bed9+1 temp.bed $CHR_SIZES $OUTPUT
rm -r tempdir temp.bed
