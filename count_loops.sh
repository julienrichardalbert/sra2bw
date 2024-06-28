#!/bin/bash
# JRA 2024

source ~/sra2bw/sra2bw_functions.sh
source activate meme

# Check the number of command line arguments
if [ $# -ne 2 ]; then
        script_name=$(basename $0)
        echo "Usage: $script_name loops.bedpe peaks.bed"
        exit 1
fi

INPUT=$1
PEAKS=$2


INPUT_COUNT=$(wc -l $INPUT | cut -f1 -d ' ')
PEAKS_COUNT=$(wc -l $PEAKS | cut -f1 -d ' ')

cut -f1-3 $INPUT > bin1.bed
cut -f4-6 $INPUT > bin2.bed
headRest 1 bin1.bed > bin1.1.bed
headRest 1 bin2.bed > bin2.1.bed
cat bin1.1.bed bin2.1.bed | sort -k1,1 -k2,2n > allbins.bed
BINS=$(wc -l allbins.bed | cut -f1 -d ' ')
PEAK_BINS=$(bedtools intersect -a allbins.bed -b $PEAKS -u | wc -l | cut -f1 -d ' ')
echo "LOOP BINS $INPUT overlapping peaks $PEAKS"
echo "BINS: $BINS"
echo "PEAK_BINS: $PEAK_BINS"
echo -e "FILE\tLOOPS_BINS\tOVERLAPPING_bins\tPEAKS" >> counts.txt
echo -e "$INPUT\t$BINS\t$PEAK_BINS\t$PEAKS_COUNT" >> counts.txt

#e.g.
#BINS: 80
#PEAK_BINS: 75
#would equate to 5 loops have one bin that doesn't intersect a peak; the rest of the loops overlap a peak at both ends.
