#!/bin/bash
source /home/robin/sra2bw/sra2bw_functions.sh

function doThis {
#  FILE=/data/fastq/D0_WT_HiChIP_rep1-2_R1.fastq.gz
  FILE=$1
  PREFIX=$(basename $FILE)
  FILE_R1="$FILE"_R1.fastq.gz
  FILE_R2="$FILE"_R2.fastq.gz
  printProgress "Processing files $FILE and $FILE2"
  printProgress "aligning reads in SE mode"
#  bwa mem -t 55 -A1 -B4 -E50 -L0 /data/reference_genomes/mm10/mm10.fa $FILE_R1 | samtools view -Shb -@ 55 - > /scratch/"$PREFIX"_R1.bam
#  bwa mem -t 55 -A1 -B4 -E50 -L0 /data/reference_genomes/mm10/mm10.fa $FILE_R2 | samtools view -Shb -@ 55 - > /scratch/"$PREFIX"_R2.bam
#  checkBamExists /scratch/"$PREFIX"_R1.bam
#  checkBamExists /scratch/"$PREFIX"_R2.bam

  printProgress "creating hic.h5 matrix file"
  hicBuildMatrix  \
    --samFiles /scratch/"$PREFIX"_R1.bam /scratch/"$PREFIX"_R2.bam \
    --binSize 5000 \
    --restrictionCutFile /data/jra/ARIMA/cut_sites_GATC.bed /data/jra/ARIMA/cut_sites_GANTC.bed \
    --danglingSequence GATC A.TC \
    --restrictionSequence GATC GA.TC \
    -o /scratch/"$PREFIX".h5 \
    --QCfolder /scratch/"$PREFIX"_hichipQC

  checkFileExists /scratch/"$PREFIX".h5

}

#doThis /data/fastq/D0_TKO_HiChIP_rep1-2
#doThis /data/fastq/D0_WT_HiChIP_rep1-2
#doThis /data/fastq/D4_TKO_HiChIP_rep1-2
doThis /data/fastq/D4_WT_HiChIP_rep1-2

#doThis /data/fastq/D0_TKO_HiChIP_shallow_intermediate_rep1
#doThis /data/fastq/D0_TKO_HiChIP_shallow_intermediate_rep2
#doThis /data/fastq/D0_WT_HiChIP_shallow_intermediate_rep1
#doThis /data/fastq/D0_WT_HiChIP_shallow_intermediate_rep2
#doThis /data/fastq/D4_TKO_HiChIP_shallow_intermediate_rep1
#doThis /data/fastq/D4_TKO_HiChIP_shallow_intermediate_rep2
#doThis /data/fastq/D4_WT_HiChIP_shallow_intermediate_rep1
#doThis /data/fastq/D4_WT_HiChIP_shallow_intermediate_rep2

#doThis /data/fastq/V6.5_mESC_serum_cohesin_HiChIP_bioRep1_techRep1_Mumbach2016_GSM2138328
#doThis /data/fastq/V6.5_mESC_serum_cohesin_HiChIP_bioRep1_techRep2_Mumbach2016_GSM2138329
#doThis /data/fastq/V6.5_mESC_serum_cohesin_HiChIP_bioRep2_techRep1_Mumbach2016_GSM2138330
#doThis /data/fastq/V6.5_mESC_serum_cohesin_HiChIP_bioRep2_techRep2_Mumbach2016_GSM2138331

