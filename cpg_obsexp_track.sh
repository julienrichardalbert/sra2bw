#!/bin/bash
# JRA 2024

# Check the number of command line arguments
if [ $# -ne 1 ]; then
        script_name=$(basename $0)
        echo "Usage: $script_name reference_genome.fa"
	echo "Requires: faidx, bedtools, bedGraphToBigWig, python with pandas installed"
        exit 1
fi



GENOME_PATH=$1
GENOME=$(basename $GENOME_PATH)
echo "Calculating CpG density track for genome: $GENOME_PATH"
faidx  $GENOME_PATH -i chromsizes > ${GENOME//.fa/.sizes}
echo "Making windows"
bedtools makewindows -w 100 -g ${GENOME//.fa/.sizes} > ${GENOME//.fa/.100w}
bedtools makewindows -w 1000 -s 100 -g ${GENOME//.fa/.sizes} > ${GENOME//.fa/.1000w.100s}
echo "Counting CpGs"
bedtools nuc -fi $GENOME_PATH -bed  ${GENOME//.fa/.1000w.100s}  -pattern CG -C > ${GENOME//.fa/.1000w.100s.data}
echo "Smoothing CpG obs/exp signals"
awk '{OFS=FS="\t"}{ if ($7*$8==0) print $1, $2, $3, "0";  else print $1, $2, $ 3, ($13*12*100)/($7*$8)'} ${GENOME//.fa/.1000w.100s.data}  > ${GENOME//.fa/.1000w.100s.data.obs.exp}
# would be nice to combine regions of contiguous 0s here
bedtools map -a ${GENOME//.fa/.100w} -b ${GENOME//.fa/.1000w.100s.data.obs.exp} -c 4 -o mean | sort -k1,1 -k2,2n | awk '{OFS=FS="\t"}{ if ($4!=0) print $0}' > ${GENOME//.fa/.1000w.100s.data.obs.exp.smoothed}
echo "Converting to bigwig format"
bedGraphToBigWig ${GENOME//.fa/.1000w.100s.data.obs.exp.smoothed} ${GENOME//.fa/.sizes} ${GENOME//.fa/CpG.obsexp.1000w.100s.bw}
echo "Done!"
rm ${GENOME//.fa/.sizes} ${GENOME//.fa/.100w} ${GENOME//.fa/.1000w.100s} ${GENOME//.fa/.1000w.100s.data} ${GENOME//.fa/.1000w.100s.data.obs.exp} ${GENOME//.fa/.1000w.100s.data.obs.exp.smoothed}
