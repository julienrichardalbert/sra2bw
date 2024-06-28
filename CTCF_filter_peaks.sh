#!/bin/bash
# JRA 2024

source ~/sra2bw/sra2bw_functions.sh
source activate meme

# Check the number of command line arguments
if [ $# -ne 3 ]; then
        script_name=$(basename $0)
        echo "Usage: $script_name peaks.bed reference.fasta motif.meme"
        exit 1
fi

INPUT=$(basename $1)
REFERENCE=$2
MOTIF=$3

printProgress "Processing file $INPUT and motif $MOTIF"
# remove header, keep first 3 columns, sort by coordinate
grep -v "#" $INPUT | cut -f1-3 | sort -k1,1 -k2,2n > $INPUT"_tmp"

# pull out the underlying sequences
bedtools getfasta -fi $REFERENCE -bed $INPUT"_tmp" -fo "${INPUT%.*}".fasta
rm $INPUT"_tmp"

printProgress "Running FIMO..."
fimo --thresh 1e-2 --max-stored-scores 9999999 -o "${INPUT%.*}"_"${MOTIF%.*}"_fimo $MOTIF "${INPUT%.*}".fasta

awk '{OFS=FS="\t"}{ if ($8 < 0.0005623413252) print $0}' "${INPUT%.*}"_"${MOTIF%.*}"_fimo/best_site.narrowPeak | cut -f1-3,6,8 | sort -k1,1 -k2,2n > "${INPUT%.*}"_1e325_"${MOTIF%.*}"_motifs.bed
PEAK_COUNT=$(wc -l $INPUT | cut -f1 -d ' ')
MOTIF_COUNT=$(wc -l "${INPUT%.*}"_1e325_"${MOTIF%.*}"_motifs.bed | cut -f1 -d ' ')
PERCENTAGE=$(echo $(( 100 * $MOTIF_COUNT/$PEAK_COUNT )))

echo "Found $MOTIF_COUNT motifs in $PEAK_COUNT peaks using a pval cutoff of 1e-3.25"
echo "Percent of peaks with motif: $PERCENTAGE"
printProgress "Outputting a list of peaks with at least one significant motif"
bedtools intersect -wa -a $INPUT -b "${INPUT%.*}"_1e325_"${MOTIF%.*}"_motifs.bed > "${INPUT%.*}"_1e325_"${MOTIF%.*}"_peaks.bed
rm -r "${INPUT%.*}"_"${MOTIF%.*}"_fimo
rm "${INPUT%.*}".fasta
printProgress "Done!"
