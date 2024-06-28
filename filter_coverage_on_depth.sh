#!/bin/bash
# JRA 2022

# Check the number of command line arguments
if [ $# -ne 2 ]; then
        script_name=$(basename $0)
        echo "Usage: $script_name inputfile.txt depth"
        echo "Example: $script_name WT2i_K27ac_MCnT_rep1_2208_trimV6_mm10_SEPE_raw_rmDup.CpG_report_mergeTwoStrands.txt 5"
        exit 1
fi

MINIMUM_DEPTH=$2
INPUT_FILE=$1

awk -v MIN_DEPTH=$MINIMUM_DEPTH '
        BEGIN {
                print "#chr\tstart\tend\tpercentMeth\tmeth\tunmeth"
                FS = "\t"
                OFS = "\t"
                if (MIN_DEPTH < 1)
                        MIN_DEPTH = 1
        }
$4 + $5 >= MIN_DEPTH {
                print $1, $2-1, $2, $3, $4, $5
}' $INPUT_FILE \
| grep -v "J0" \
>  $INPUT_FILE"_x"$MINIMUM_DEPTH
