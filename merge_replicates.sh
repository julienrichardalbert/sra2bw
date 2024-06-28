#!/bin/bash
# JRA 2021

# Check the number of command line arguments
if [ $# -ne 3 ]; then
	script_name=$(basename $0)
	echo "Usage: $script_name extension common_string STAY/GO"
  echo "Example: $script_name.sh .bw rep GO"
  echo "Results: Merging BC_E7.5_epi_RNA_Hanna2019_*.bw -> BC_E7.5_epi_RNA_Hanna2019_rep1-6.bw"
  echo "Example: $script_name .bam Hanna2019_ STAY"
  echo "Results: Merging BC_E7.5_epi_RNA_Hanna2019_*.bw -> BC_E7.5_epi_RNA_Hanna2019_rep1-6.bw"
	exit 1
fi

EXTENSION=$1
LAST_COMMON_STRING=$2
SHOULD_I_STAY_OR_SHOULD_I_GO=$3


# Make a list of file families that are made up of replicates
if [[ $SHOULD_I_STAY_OR_SHOULD_I_GO == "STAY" ]]; then
  PREFIX=$(for x in *$LAST_COMMON_STRING*$EXTENSION; do # get all files that end with EXTENSION
    y=${x//$LAST_COMMON_STRING*/$LAST_COMMON_STRING} # for each file, find the last common string
    echo $y # print the last common string
  done | uniq) # get list of unique strings

#same as above but remove the last common string for final output filename
elif [[ $SHOULD_I_STAY_OR_SHOULD_I_GO == "GO" ]]; then
PREFIX=$(for x in *$EXTENSION; do y=${x//$LAST_COMMON_STRING*}; echo $y; done | uniq)

else
    echo "Argument 3 needs to be STAY or GO. STAY if you'd rather keep the last common string. GO if not. Retry"
    exit 1
fi

#  Now, we have a list of all the unique file names (families)
#  Go through that list of families one by one and concatenate the members into their families

for n in $PREFIX; do
  REP_COUNT_TMP=$(ls -l $n*$EXTENSION | wc -l )
  REP_COUNT=$(basename $REP_COUNT_TMP)
#  echo "prefix:$n replicate count:$REP_COUNT"
  OUTPUT_FILE="$n"rep1-"$REP_COUNT""$EXTENSION"


  echo "Merging $n*$EXTENSION -> $OUTPUT_FILE"
  if [[ $EXTENSION == ".bam" ]]; then
    samtools merge $OUTPUT_FILE $n*$EXTENSION -@ 8
    samtools index -@ 8 $OUTPUT_FILE
  else
    cat $n*$EXTENSION > $OUTPUT_FILE
  fi
done


# now run the type-specific script, like merge CpG reports or BAMs
# do in separate script... enough debauchery here.
