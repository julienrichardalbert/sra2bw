#!/bin/bash
# JRA 2022
# Check the number of command line arguments
source /home/robin/sra2bw/sra2bw_functions.sh

if [ $# -ne 2 ]; then
        script_name=$(basename $0)
        echo "Usage: $script_name /full/path/to/input_bam chromosome.sizes"
        exit 1
fi

# REQUIREMENTS
# conda install -c bioconda seacr
# conda install -c bioconda bedtools
# conda install -c bioconda samtools

THREADS=4

function bam2bed_SEACR {

  local PATH_TO_FILE=$1
  local FILE=$(basename $PATH_TO_FILE)
  local SIZES=$2

  printProgress "[bam2bed] Started on file: $FILE"

#  cd /scratch/
#  ln -s $PATH_TO_FILE .
  checkFileExists $FILE
  $BEDTOOLS intersect -v -a $FILE -b /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed > $FILE"_tmp0.bam"
  $SAMTOOLS sort -@ $THREADS -n $FILE"_tmp0.bam" -o $FILE"_tmp.bam"
  $BEDTOOLS bamtobed -bedpe -i $FILE"_tmp.bam" > $FILE"_tmp.bed"
  awk '$1==$4 && $6-$2 < 1000 && $2>0 {print $0}' $FILE"_tmp.bed" > $FILE"_tmp.clean.bed"
  cut -f 1,2,6 $FILE"_tmp.clean.bed" | sort -k1,1 -k2,2n -k3,3n > $FILE"_tmp.clean.fragments.bed"
  $BEDTOOLS genomecov -bg -i $FILE"_tmp.clean.fragments.bed" -g $SIZES > ${FILE//.bam/.fragments.bedgraph}
  rm $FILE"_tmp0.bam" $FILE"_tmp.bam" $FILE"_tmp.bed" $FILE"_tmp.clean.bed" $FILE"_tmp.clean.fragments.bed"
  checkFileExists ${FILE//.bam/.fragments.bedgraph}
  printProgress "[bam2bed] Finished converting $FILE -> ${FILE//.bam/.fragments.bedgraph}"
#  unlink $FILE

}

bam2bed_SEACR $1 $2
