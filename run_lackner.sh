#!/bin/bash

#SBATCH --ntasks=24                    # number of MPI processes
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes

./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_10h_RNA_Lackner2021_GSM4323423
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_14h_RNA_Lackner2021_GSM4323425
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_18h_RNA_Lackner2021_GSM4323427
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_20h_RNA_Lackner2021_GSM4323428
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_22h_RNA_Lackner2021_GSM4323429
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_24h_RNA_Lackner2021_GSM4323430
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_26h_RNA_Lackner2021_GSM4323431
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_2h_RNA_Lackner2021_GSM4323419
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_2i_LIF_RNA_Lackner2021_rep2_GSM4323439
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_2i_RNA_Lackner2021_rep1_GSM4323436
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_32h_RNA_Lackner2021_rep1_GSM4323434
./run_RNAseq_mm10_CPM_v5_sbatch.sh /data/fastq/Lackner2021/RC9_mESCs_4h_RNA_Lackner2021_GSM4323420
