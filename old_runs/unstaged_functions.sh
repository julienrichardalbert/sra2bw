#! /bin/bash
# JRA 2022

# example on an output bam
# $FILE.bam


function bam2bed_SEACR {

  # this function convert a bowtie2-aligned BAM file into a SEACR2-compatible bedGraph file
  # annoyingly, we must include a chromosome sizes file for the bedtools genomcov step.
  CHROM_SIZES="$PATH_TO_GENOME".sizes
  checkFileExists $CHROM_SIZES
  printProgress "[bam2bed] Started on file: $FILE"
  checkFileExists "$FILE".bam
  printProgress "[bam2bed] Removing reads aligned to blacklisted  regions $BLACKLIST"
  $BEDTOOLS intersect -v -a "$FILE".bam -b $BLACKLIST > $FILE"_tmp0.bam"
  printProgress "[bam2bed] Sorting bam file by read name"
  $SAMTOOLS sort -@ $THREADS -n $FILE"_tmp0.bam" -o $FILE"_tmp.bam"
  printProgress "[bam2bed] Converting bam to paired-end bed file"
  $BEDTOOLS bamtobed -bedpe -i $FILE"_tmp.bam" > $FILE"_tmp.bed"
  printProgress "[bam2bed] Keeping intrachromosomal fragments with insert size <1000"
  awk '$1==$4 && $6-$2 < 1000 && $2>0 {print $0}' $FILE"_tmp.bed" > $FILE"_tmp.clean.bed"
  printProgress "[bam2bed] Sorting bed file"
  cut -f 1,2,6 $FILE"_tmp.clean.bed" | sort -k1,1 -k2,2n -k3,3n > $FILE"_tmp.clean.fragments.bed"
  printProgress "[bam2bed] Converting paired-end bed file to genome coverage bedgraph"
  $BEDTOOLS genomecov -bg -i $FILE"_tmp.clean.fragments.bed" -g $CHROM_SIZES > "$FILE".fragments.bedgraph
  printProgress "[bam2bed] Removing temporary bam2bed function files"
  rm $FILE"_tmp0.bam" $FILE"_tmp.bam" $FILE"_tmp.bed" $FILE"_tmp.clean.bed" $FILE"_tmp.clean.fragments.bed"
  checkFileExists "$FILE".fragments.bedgraph
  printProgress "[bam2bed] Finished converting $FILE -> "$FILE".fragments.bedgraph"
}






function callPeaks {
  PEAK_WIDTH=$1 # broad or narrow

  if [[ $FILE == *"ChIP"* ]]; then
    if [[ $PEAK_WIDTH == "narrow" ]]; then
      printProgress "Calling peaks using MACS2 narrow"
      $MACS2_ACTIVATE
      $MACS2 callpeak -t "$FILE".bam
      $MACS2_DEACTIVATE
      mv NA_peaks.narrowPeak "$FILE"_macs2_narrow_peaks.bed

    elif [[ $PEAK_WIDTH == "broad" ]]; then
      printProgress "Calling peaks using MACS2 broad"
      $MACS2_ACTIVATE
      # run it with mouse genome size, because mice are the best
      $MACS2 callpeak -t "$FILE".bam --broad -g 1.87e9
      $MACS2_DEACTIVATE
      mv NA_peaks.broadPeak "$FILE"_macs2_broad_peaks.bed
    else
      printProgress "please use 'narrow' or 'broad' for peak width option, e.g."
      printProgress "callPeaks narrow"
      return 0
    fi

  elif [[ $FILE == *"CUT"* ]] || [[ $FILE == *"CnR"* ]] || [[ $FILE == *"CnT"* ]] ; then
    printProgress "Calling peaks using SEACR"
    bam2bed_SEACR
    $SEACR_ACTIVATE
    $SEACR "$FILE".fragments.bedgraph 0.1 non stringent "$FILE"_seacr_0.1_non_peaks
    $SEACR "$FILE".fragments.bedgraph 0.01 non stringent "$FILE"_seacr_0.01_non_peaks
    $SEACR "$FILE".fragments.bedgraph 0.001 non stringent "$FILE"_seacr_0.001_non_peaks
    $SEACR "$FILE".fragments.bedgraph 0.0001 non stringent "$FILE"_seacr_0.0001_non_peaks
#    $SEACR "$FILE".fragments.bedgraph 0.00001 non stringent "$FILE"_seacr_0.00001_non_peaks
    for o in *peaks.stringent.bed; do mv $o ${o//peaks.stringent/stringent_peaks}; done
#    mv tmp.stringent.bed "$FILE"_seacr_0.01_non_stringent_peaks.bed
    $SEACR_DEACTIVATE
  else
      printProgress "File name did not contain ChIP, CUT, CnR or CnT, exiting"
      return 0
  fi
}


function calculateEnrichment {

  if [[ $PEAK_WIDTH == "narrow" ]]; then
    # this creates a matrix file of 200 columns
    UPSTREAM=10000
    DOWNSTREAM=10000
    MATRIX_BIN_SIZE=100
  else # broad
    UPSTREAM=50000
    DOWNSTREAM=50000
    MATRIX_BIN_SIZE=500
  fi

  # oh my god these fucking programs that don't recognize >3 column bed files
  for p in "$FILE"*peaks.bed; do cut -f1-3 $p > "$p"3; done
  PEAKS=$(ls "$FILE"*peaks.bed3)

  printProgress "Calculating the fraction of reads in peaks (FRiP)"
  ALN_COUNT=$($SAMTOOLS view -c "$FILE".bam)

  printProgress "Number of alignments in $FILE.bam: $ALN_COUNT"
  for p in $PEAKS; do
    FRIP=$($BEDTOOLS intersect -u -a "$FILE".bam -b $p | $SAMTOOLS view -c)
    printProgress "Number of alignments in peak file $p: $FRIP"
  done

  printProgress "Converting bam to bigwig for enrichment calculation"
  bamToBigWigDeeptoolsCPMsmoothKeepDup 10 0
  checkFileExists "$FILE"_kpDup_q"$MIN_MAPQ"_b10_s0_CPM.bw

  printProgress "Computing matrix using deeptools over $PEAKS using file "$FILE"_kpDup_q"$MIN_MAPQ"_b10_s0_CPM.bw"
  $COMPUTEMATRIX reference-point \
    --referencePoint center \
    -b $UPSTREAM \
    -a $DOWNSTREAM \
    --binSize $MATRIX_BIN_SIZE \
    --missingDataAsZero \
    -R $PEAKS \
    -S "$FILE"_kpDup_q"$MIN_MAPQ"_b10_s0_CPM.bw \
    --numberOfProcessors $THREADS \
    -o "$FILE"_matrix.gz \
    --outFileNameMatrix "$FILE"_matrix.txt
    rm "$FILE"_kpDup_q"$MIN_MAPQ"_b10_s0_CPM.bw
  checkFileExists "$FILE"_matrix.txt

  printProgress "Generating heatmap"
  $PLOTHEATMAP -m  "$FILE"_matrix.gz \
    --colorMap RdYlBu_r \
    --missingDataColor "#000000" \
    --yMin 0 \
    -out "$FILE"_heatmap.pdf
  checkFileExists "$FILE"_heatmap.pdf

  PEAK_COUNT=$(wc -l $PEAKS)
  AREA=$(awk '{ SUM+=($3-$2)} END {print SUM }' $PEAKS)

  printProgress "Calculating enrichment of IP over background"
  ENRICHMENT=$(cat "$FILE"_matrix.txt | awk -F $'\t' 'BEGIN {background = 0; peak = 0} {background = background + $7; peak = peak + $100} END {enrichment = peak / background; printf("%.2f", enrichment)}');
  printProgress ""
  printProgress "---- Some stats ----"
  printProgress "Count: $PEAK_COUNT"
  printProgress "Area: $AREA"
  printProgress "Enrichment: $ENRICHMENT"
  printProgress "Number of alignments in $FILE.bam: $ALN_COUNT"
  printProgress "FRiP: $FRIP"

}




function epicypherSpike {

  printProgress "Counting the number of synthetic spike-in nucleosomes in dataset"
  printProgress "File: "$FILE"_R1.fastq.gz"
	gunzip -c "$FILE"_R1.fastq.gz > sample1_R1.fastq
  for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA ;
  do
  	grep -c $barcode sample1_R1.fastq  | tee -a $LOG_FILE
  done
  rm sample1_R1.fastq

  printProgress "File: "$FILE"_R2.fastq.gz"
  gunzip -c "$FILE"_R2.fastq.gz > sample1_R2.fastq
	for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA ;
  do
  	grep -c $barcode sample1_R2.fastq | tee -a $LOG_FILE
	done
  rm sample1_R2.fastq

  # Barcode identities
  # Unmodified (A & B)
  # TTCGCGCGTAACGACGTACCGT
  # CGCGATACGACCGCGTTACGCG

  # H3K4me1 (A & B)
  # CGACGTTAACGCGTTTCGTACG
  # CGCGACTATCGCGCGTAACGCG

  # H3K4me2 (A & B)
  # CCGTACGTCGTGTCGAACGACG
  # CGATACGCGTTGGTACGCGTAA

  # H3K4me3 (A & B)
  # TAGTTCGCGACACCGTTCGTCG
  # TCGACGCGTAAACGGTACGTCG

  # H3K9me1 (A & B)
  # TTATCGCGTCGCGACGGACGTA
  # CGATCGTACGATAGCGTACCGA

  # H3K9me2 (A & B)
  # CGCATATCGCGTCGTACGACCG
  # ACGTTCGACCGCGGTCGTACGA

  # H3K9me3 (A & B)
  # ACGATTCGACGATCGTCGACGA
  # CGATAGTCGCGTCGCACGATCG

  # H3K27me1 (A & B)
  # CGCCGATTACGTGTCGCGCGTA
  # ATCGTACCGCGCGTATCGGTCG

  # H3K27me2 (A & B)
  # CGTTCGAACGTTCGTCGACGAT
  # TCGCGATTACGATGTCGCGCGA

  # H3K27me3 (A & B)
  # ACGCGAATCGTCGACGCGTATA
  # CGCGATATCACTCGACGCGATA

  # H3K36me1 (A & B)
  # CGCGAAATTCGTATACGCGTCG
  # CGCGATCGGTATCGGTACGCGC

  # H3K36me2 (A & B)
  # GTGATATCGCGTTAACGTCGCG
  # TATCGCGCGAAACGACCGTTCG

  # H3K36me3 (A & B)
  # CCGCGCGTAATGCGCGACGTTA
  # CCGCGATACGACTCGTTCGTCG

  # H4K20me1 (A & B)
  # GTCGCGAACTATCGTCGATTCG
  # CCGCGCGTATAGTCCGAGCGTA

  # H4K20me2 (A & B)
  # CGATACGCCGATCGATCGTCGG
  # CCGCGCGATAAGACGCGTAACG

  # H4K20me3 (A & B)
  # CGATTCGACGGTCGCGACCGTA
  # TTTCGACGCGTCGATTCGGCGA

}

