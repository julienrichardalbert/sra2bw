#!/bin/bash
SEACR="/home/robin/bin/miniconda3/bin/SEACR_1.3.sh"

for FILE in *bedgraph; do
	$SEACR "$FILE" 0.1 non stringent "$FILE"_seacr_0.1_non_peaks
	$SEACR "$FILE" 0.01 non stringent "$FILE"_seacr_0.01_non_peaks
	$SEACR "$FILE" 0.001 non stringent "$FILE"_seacr_0.001_non_peaks
	$SEACR "$FILE" 0.0001 non stringent "$FILE"_seacr_0.0001_non_peaks
	$SEACR "$FILE" 0.00001 non stringent "$FILE"_seacr_0.00001_non_peaks

        $SEACR "$FILE" 0.1 non relaxed "$FILE"_seacr_0.1_non_peaks
        $SEACR "$FILE" 0.01 non relaxed "$FILE"_seacr_0.01_non_peaks
        $SEACR "$FILE" 0.001 non relaxed "$FILE"_seacr_0.001_non_peaks
        $SEACR "$FILE" 0.0001 non relaxed "$FILE"_seacr_0.0001_non_peaks
        $SEACR "$FILE" 0.00001 non relaxed "$FILE"_seacr_0.00001_non_peaks

done
