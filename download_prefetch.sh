#! /bin/bash
# by ABB (UBC). Couldn't have done it without JRA. November 2018
# takes filename as only aggument (tab delimited SRA + description)
# Updated 3 December 2021 to fix the PAIRED variable (see commented line for old version)
# Updated 26 January 2022 to add NCBI API key after getting 403: forbidden errors
# Updated 21 June 2022 to increase the max size of downloads from 20 to 80Gb. Also updated sra-toolkit and swiched to fasterq-dump
# Updated 13 May 2022 to compress the fastq in a different directory; so that a spinning disk doesn't have to read and write
# Updated 27 Feb 2024 undo previous upudate

source /home/robin/sra2bw/sra2bw_functions.sh
export NCBI_API_KEY=5001bf72922d25a715ac5b4203dc71863d08 # JRA NCBI API key (how many acronyms is too many acronyms?)
SRA_DEFAULT_PATH="/home/robin/ncbi/public/sra"


function download () {
    OUTPUT="/data/fastq"
    SCRATCH="/scratch"
#    SCRATCH2="/scratch2"
    FILE=$(basename $1)
    LOG_FILE="$FILE"_"$(date '+%y-%m-%d')"_log.txt
    SCRATCH_DOWNLOADS=$SCRATCH/"${FILE//.txt/}"_"$(date '+%y-%m-%d')_downloads"
    mkdir $SCRATCH_DOWNLOADS
    cp $1 $SCRATCH_DOWNLOADS/
    cd $SCRATCH_DOWNLOADS
    checkFileExists $FILE

    while read n
    do
	#declares m as each line of the SRA:name file as an array
        declare -a m=($n)
	#populate TMP with each entry of the first column of the SRA:name file as a space delimited list
        TMP=$TMP" "${m[0]}
    done < $FILE

# TMP becomes the list of SRA codes you're going to download
    for n in $TMP
    do

# populates SRA:name code from table
        SRACODE=$(grep -e "$n" $FILE | cut -f1)
        NAME=$(grep -e "$n" $FILE | cut -f2)

        echo "downloading "$SRACODE" at [`date`] ..."
#       $ESEARCH -db sra -query $SRACODE | $EFETCH -format runinfo | cut -d ',' -f 10 | grep https | xargs wget
        SRR_CODES=$($ESEARCH -db sra -query $SRACODE | $EFETCH -format runinfo | cut -d ',' -f 1 | headRest 1 stdin)
        for j in $SRR_CODES ; do
            $PREFETCH $j --max-size 80G --output-directory $SCRATCH_DOWNLOADS
#            $VALIDATE $SCRATCH/"${FILE//.txt/}"_"$(date '+%y-%m-%d')_downloads"/$j
            mv $j/$j".sra" $SCRATCH_DOWNLOADS/$j".sra"
            rm -r $j
#           mv $SRA_DEFAULT_PATH/$j".sra" ./
        done

        echo "finished "$SRACODE" -> "$NAME" download! at [`date`] "
#       PAIRED=`$ESEARCH -db sra -query $SRACODE | $EFETCH -format runinfo | cut -d ',' -f 16 | tail -n 2 | head -n 1`
        PAIRED=`$ESEARCH -db sra -query $SRACODE | $EFETCH -format runinfo | cut -d ',' -f 16 | tail -n 1`

       echo "PAIRED VARIABLE: $PAIRED"
        if [[ $PAIRED = SINGLE ]] ; then
            PAIRED_END=false
            echo "Data are single-end."
        else
            PAIRED_END=true
            echo "Data are paired-end."
        fi

        if $PAIRED_END ; then
            FILE_FASTQ1="$NAME"_R1.fastq
            FILE_FASTQ2="$NAME"_R2.fastq

            echo "Dumping .fastq files at [`date`] ..."
            for SRA_FILE in [DS]RR*
            do
#                $FASTQ_DUMP --split-files ./"$SRA_FILE"
                $FASTERQ_DUMP --split-files -e $THREADS -m $THREAD_MEM ./"$SRA_FILE"
            done
            echo "concatenating .fastq files at [`date`] ..."
            COUNT=`ls -1 [DS]RR*_*1.fastq | wc -l`
            echo "count: $COUNT"
            if [[ $COUNT > 1 ]] ; then # If only single file, move, don't copy
                cat [DS]RR*_*1.fastq > "$NAME"_R1.fastq
                cat [DS]RR*_*2.fastq > "$NAME"_R2.fastq
                checkFileExists "$NAME"_R1.fastq && checkFileExists "$NAME"_R2.fastq
                rm -r [DS]RR*
            else
                mv [DS]RR*_*1.fastq "$NAME"_R1.fastq
                mv [DS]RR*_*2.fastq "$NAME"_R2.fastq
                checkFileExists "$NAME"_R1.fastq && checkFileExists "$NAME"_R2.fastq
                rm -r [DS]RR*
            fi
            echo "compressing .fastq files at [`date`] ..."
#            gzip "$NAME"_R1.fastq
#            gzip "$NAME"_R2.fastq
            pigz -c -p $THREADS "$NAME"_R1.fastq > "$NAME"_R1.fastq.gz
            rm "$NAME"_R1.fastq
            pigz -c -p $THREADS "$NAME"_R2.fastq > "$NAME"_R2.fastq.gz
            rm "$NAME"_R2.fastq
            mv "$NAME"_R*.fastq.gz $OUTPUT/

        else #Single-End
            FILE_FASTQ="$NAME".fastq

            echo "dumping .fastq files at [`date`] ..."
            for SRA_FILE in [DS]RR*
            do
#                $FASTQ_DUMP ./"$SRA_FILE"
                $FASTERQ_DUMP -e $THREADS -m $THREAD_MEM ./"$SRA_FILE"

            done
            echo "concatenating .fastq files at [`date`] ..."
            COUNT=`ls -1 [DS]RR*.fastq | wc -l`
            if [[ $COUNT > 1 ]] ; then # If only single file, move, don't copy
                cat [DS]RR*.fastq > "$NAME".fastq
                checkFileExists "$NAME".fastq
                rm -r [DS]RR*
            else
                mv [DS]RR*.fastq "$NAME".fastq
                checkFileExists "$NAME".fastq
                rm -r [DS]RR*
            fi
            echo "compressing .fastq files at [`date`] ..."
            gzip "$NAME".fastq
            mv "$NAME".fastq.gz $OUTPUT/
        fi
    done
    echo "All done! at [`date`]"
}

download $1
