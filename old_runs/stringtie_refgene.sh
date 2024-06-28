#!/bin/bash
source /home/robin/sra2bw/epipax.config

stringtie /data/bams/E14_TKO_D0_RNA_D129T05_rep1_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_TKO_D0_RNA_D129T05_rep1_trimV5_mm10 -o /data/bams/E14_TKO_D0_RNA_D129T05_rep1_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_TKO_D0_RNA_D129T13_rep2_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_TKO_D0_RNA_D129T13_rep2_trimV5_mm10 -o /data/bams/E14_TKO_D0_RNA_D129T13_rep2_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_TKO_D2_RNA_D129T06_rep1_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_TKO_D2_RNA_D129T06_rep1_trimV5_mm10 -o /data/bams/E14_TKO_D2_RNA_D129T06_rep1_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_TKO_D2_RNA_D129T14_rep2_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_TKO_D2_RNA_D129T14_rep2_trimV5_mm10 -o /data/bams/E14_TKO_D2_RNA_D129T14_rep2_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_TKO_D4_RNA_D129T07_rep1_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_TKO_D4_RNA_D129T07_rep1_trimV5_mm10 -o /data/bams/E14_TKO_D4_RNA_D129T07_rep1_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_TKO_D4_RNA_D129T15_rep2_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_TKO_D4_RNA_D129T15_rep2_trimV5_mm10 -o /data/bams/E14_TKO_D4_RNA_D129T15_rep2_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_TKO_D7_RNA_D129T08_rep1_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_TKO_D7_RNA_D129T08_rep1_trimV5_mm10 -o /data/bams/E14_TKO_D7_RNA_D129T08_rep1_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_TKO_D7_RNA_D129T16_rep2_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_TKO_D7_RNA_D129T16_rep2_trimV5_mm10 -o /data/bams/E14_TKO_D7_RNA_D129T16_rep2_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_WT_D0_RNA_D129T01_rep1_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_WT_D0_RNA_D129T01_rep1_trimV5_mm10 -o /data/bams/E14_WT_D0_RNA_D129T01_rep1_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_WT_D0_RNA_D129T09_rep2_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_WT_D0_RNA_D129T09_rep2_trimV5_mm10 -o /data/bams/E14_WT_D0_RNA_D129T09_rep2_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_WT_D2_RNA_D129T02_rep1_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_WT_D2_RNA_D129T02_rep1_trimV5_mm10 -o /data/bams/E14_WT_D2_RNA_D129T02_rep1_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_WT_D2_RNA_D129T10_rep2_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_WT_D2_RNA_D129T10_rep2_trimV5_mm10 -o /data/bams/E14_WT_D2_RNA_D129T10_rep2_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_WT_D4_RNA_D129T03_rep1_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_WT_D4_RNA_D129T03_rep1_trimV5_mm10 -o /data/bams/E14_WT_D4_RNA_D129T03_rep1_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_WT_D4_RNA_D129T11_rep2_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_WT_D4_RNA_D129T11_rep2_trimV5_mm10 -o /data/bams/E14_WT_D4_RNA_D129T11_rep2_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_WT_D7_RNA_D129T04_rep1_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_WT_D7_RNA_D129T04_rep1_trimV5_mm10 -o /data/bams/E14_WT_D7_RNA_D129T04_rep1_trimV5_mm10.gtf -p 20
stringtie /data/bams/E14_WT_D7_RNA_D129T12_rep2_trimV5_mm10.bam -G /data/reference_genomes/mm10/mm10.refGene.gtf -l /data/bams/E14_WT_D7_RNA_D129T12_rep2_trimV5_mm10 -o /data/bams/E14_WT_D7_RNA_D129T12_rep2_trimV5_mm10.gtf -p 20

