#!/bin/bash


# samplename=mESC_serum_Kdm2ab_JmjCfl_UNT_ChIP_K36me2-input_Turberfield2019_rep1_GSM3615604_trimV5_mm10

samplename=alignments_peaks

#fragment lengths
samtools view -F 16 -F 4 -f 1 -S ${samplename}.bam | cut -f 9 > ${samplename}.fraglengths.txt
#calculate fragment length histogram
Rscript /scratch/fragment_sizes_tmp/fragment_length_histogram_from_raw_lengths.r ${samplename}.fraglengths.txt ${samplename}
cat ${samplename}_fragment_size_summary_stats.txt


samplename=alignments_background
samtools view -F 16 -F 4 -f 1 -S ${samplename}.bam | cut -f 9 > ${samplename}.fraglengths.txt
#calculate fragment length histogram
Rscript /scratch/fragment_sizes_tmp/fragment_length_histogram_from_raw_lengths.r ${samplename}.fraglengths.txt ${samplename}
cat ${samplename}_fragment_size_summary_stats.txt


