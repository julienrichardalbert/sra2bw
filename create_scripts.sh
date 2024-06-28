#!/bin/bash
# JRA 2021

# Check the number of command line arguments
if [ $# -ne 2 ]; then
    script_name=$(basename $0)
    echo "Usage: $script_name DNAme_file.bed5 ROI_coordinate_sorted.bed"
    exit 1
fi

SCRIPT_FILE="script_"$(date '+%y-%m-%d')".sh"

echo "bedmapCoverage.sh $1 $2 count" >> $SCRIPT_FILE
echo "bedmapCoverage.sh $1 $2 mean" >> $SCRIPT_FILE
echo "bedmapCoverage.sh $1 $2 echo-map-range" >> $SCRIPT_FILE

chmod 777 $SCRIPT_FILE

echo "created $SCRIPT_FILE successfully!"
