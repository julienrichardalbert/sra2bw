source ~/sra2bw/epipax.config

TRIMVER=v1
TRIM_PARAMS="ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36"
THREADS=40

for p in *R1.fastq.gz; do
  FILE=${p//_R1.fastq.gz}
  ($TRIMMOMATIC PE -threads $THREADS \
  -summary "$FILE"_trimSummary.txt \
  "$FILE"_R1.fastq.gz "$FILE"_R2.fastq.gz \
  "$FILE"_"$TRIMVER"_R1.fastq.gz "$FILE"_unpaired_trim_R1.fastq.gz "$FILE"_"$TRIMVER"_R2.fastq.gz "$FILE"_unpaired_trim_R2.fastq.gz \
  $TRIM_PARAMS ) 2>> "$FILE"_trimLog.txt
done


