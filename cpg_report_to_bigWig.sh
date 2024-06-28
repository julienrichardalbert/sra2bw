#!/bin/bash

# Check the number of command line arguments
if [ $# -ne 3 ]; then
	script_name=$(basename $0)
	echo "Usage: $script_name cpg_report.txt min_depth reference.sizes"
	exit 1
fi

# source /home/teamgreenberg/sra2bw/greenberg05.config

INPUT_TWOSTRANDMERGED_CPG_REPORT=$1
METH_TO_BIGWIG_MIN_DEPTH=$2
CHROM_SIZES=$3

echo "[CpGReportToBigWig] Started on file $INPUT_TWOSTRANDMERGED_CPG_REPORT with min coverage x$METH_TO_BIGWIG_MIN_DEPTH" | tee -a $INPUT_TWOSTRANDEDMERGED_CPG_REPORT"CpG_counts.txt"
awk -v MIN_DEPTH=$METH_TO_BIGWIG_MIN_DEPTH '
	BEGIN {
		print "#chr\tstart\tend\tvalue"
		FS = "\t"
		OFS = "\t"
		if (MIN_DEPTH < 1)
			MIN_DEPTH = 1
	}
$4 + $5 >= MIN_DEPTH {
		print $1, $2-1, $2, $3
}' $INPUT_TWOSTRANDMERGED_CPG_REPORT | grep -v "J02459.1" > $INPUT_TWOSTRANDMERGED_CPG_REPORT"_x"$METH_TO_BIGWIG_MIN_DEPTH

echo "converting to BigWig"
bedGraphToBigWig $INPUT_TWOSTRANDMERGED_CPG_REPORT"_x"$METH_TO_BIGWIG_MIN_DEPTH $CHROM_SIZES ${INPUT_TWOSTRANDMERGED_CPG_REPORT/.txt/_x$METH_TO_BIGWIG_MIN_DEPTH.bw}

# count number of covered CpGs
echo "[CpGReportToBigWig] Counting the number of covered CpGs by $METH_TO_BIGWIG_MIN_DEPTH" | tee -a $INPUT_TWOSTRANDMERGED_CPG_REPORT"CpG_counts.txt"
CPG_COUNT=$( grep -v "J02459.1" $INPUT_TWOSTRANDMERGED_CPG_REPORT"_x"$METH_TO_BIGWIG_MIN_DEPTH |  wc -l  |  awk '{print $1}' )
echo "$CPG_COUNT CpGs covered by $METH_TO_BIGWIG_MIN_DEPTH reads in $INPUT_TWOSTRANDMERGED_CPG_REPORT" | tee -a $INPUT_TWOSTRANDMERGED_CPG_REPORT"CpG_counts.txt"
checkFileExists ${INPUT_TWOSTRANDMERGED_CPG_REPORT/.txt/_x$METH_TO_BIGWIG_MIN_DEPTH.bw}
# rm $INPUT_TWOSTRANDMERGED_CPG_REPORT"_x"$METH_TO_BIGWIG_MIN_DEPTH
echo "[CpGReportToBigWig] $INPUT_TWOSTRANDMERGED_CPG_REPORT x$METH_TO_BIGWIG_MIN_DEPTH CpG_report -> bigWig done!" | tee -a $INPUT_TWOSTRANDMERGED_CPG_REPORT"CpG_counts.txt"
echo "" | tee -a $INPUT_TWOSTRANDMERGED_CPG_REPORT"CpG_counts.txt"



