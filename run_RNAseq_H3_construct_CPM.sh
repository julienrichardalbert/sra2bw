#!/bin/bash
# load functions coded in separate files

source /home/robin/sra2bw/sra2bw_functions.sh
SCRATCH="/scratch/"

# ---------------- CONFIGURE VARIABLES ----------------- #
FLAG=1540
MIN_MAPQ=255
THREADS=3
THREAD_MEM=8G

# ------------------ CALL FUNCTIONS -------------------- #
setupVariables $1
trimReads 5
alignSTAR /data/reference_genomes/H3_ROSA
groomSam
rm *sam

# Remove duplicate reads and low-quality mapping reads
echo "-q $MIN_MAPQ -F $FLAG "$FILE".bam" > "$FILE"_counts.txt
samtools view -bh -q $MIN_MAPQ -F $FLAG "$FILE".bam > "$FILE"_q"$MIN_MAPQ"_F"$FLAG".bam
samtools coverage "$FILE"_q"$MIN_MAPQ"_F"$FLAG".bam >> "$FILE"_counts.txt

# Keep all aligned reads
echo "-q 1 -F none "$FILE".bam" >> "$FILE"_counts.txt
samtools view -bh -q 1 "$FILE".bam > "$FILE"_tmp.bam
samtools coverage "$FILE"_tmp.bam >> "$FILE"_counts.txt
rm "$FILE"_tmp.bam


mv "$FILE"_counts.txt /scratch/ROSA/
mv "$FILE".bam /scratch/ROSA/

rm -r $SCRATCH/$FILE
