#! /bin/bash
# by ABB (UBC). Couldn't have done it without JRA. November 2018
# takes filename as only aggument (tab delimited SRA + description)
# Updated 3 December 2021 to fix the PAIRED variable (see commented line for old version)
# Updated 26 January 2022 to add NCBI API key after getting 403: forbidden errors
# Updated 7 June 2022 to add the --no-check-certificate parameter to wget

source /home/robin/sra2bw/sra2bw_functions.sh
export NCBI_API_KEY=5001bf72922d25a715ac5b4203dc71863d08 # JRA NCBI API key (how many acronyms is too many acronyms?)
export api_key=5001bf72922d25a715ac5b4203dc71863d08

function download () {
    OUTPUT="/data/fastq"
    SCRATCH="/scratch"
#    SCRATCH2="/scratch2"
    FILE=$(basename $1)
    mkdir $SCRATCH/"${FILE//.txt/}"_"$(date '+%y-%m-%d')_downloads"
    cp $1 $SCRATCH/"${FILE//.txt/}"_"$(date '+%y-%m-%d')_downloads"
    cd $SCRATCH/"${FILE//.txt/}"_"$(date '+%y-%m-%d')_downloads"
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

        if [[ $NAME == *[rR]ep* ]] ; then
            X=${NAME%_*}
            FOLDER_NAME=${X##*_}
        else
            FOLDER_NAME=${NAME##*_}
        fi
        echo "downloading "$SRACODE" at [`date`] ..."

        # NEW add -q and --no-check-certificate parameters to wget. Hope to resolve the 403 forbidden error
        $ESEARCH -db sra -query $SRACODE | $EFETCH -format runinfo | cut -d ',' -f 10 | grep https | xargs wget --no-check-certificate

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
                $FASTQ_DUMP --split-files ./"$SRA_FILE"
            done
            echo "concatenating .fastq files at [`date`] ..."
            COUNT=`ls -1 [DS]RR*_*1.fastq | wc -l`
            echo "count: $COUNT"
            if [[ $COUNT > 1 ]] ; then # If only single file, move, don't copy
                cat [DS]RR*_*1.fastq > "$NAME"_R1.fastq
                cat [DS]RR*_*2.fastq > "$NAME"_R2.fastq
                checkFileExists "$NAME"_R1.fastq && checkFileExists "$NAME"_R2.fastq
                rm [DS]RR*
            else
                mv [DS]RR*_*1.fastq "$NAME"_R1.fastq
                mv [DS]RR*_*2.fastq "$NAME"_R2.fastq
                checkFileExists "$NAME"_R1.fastq && checkFileExists "$NAME"_R2.fastq
                rm [DS]RR*
            fi
            echo "compressing .fastq files at [`date`] ..."

            # try this
            gzip "$NAME"_R1.fastq
#            rm "$NAME"_R1.fastq
            gzip "$NAME"_R2.fastq
#            rm "$NAME"_R2.fastq
            mv $SCRATCH/"$NAME"_R*.fastq.gz $OUTPUT/
            # end try this

        else #Single-End
            FILE_FASTQ="$NAME".fastq

            echo "dumping .fastq files at [`date`] ..."
            for SRA_FILE in [DS]RR*
            do
                $FASTQ_DUMP ./"$SRA_FILE"
            done
            echo "concatenating .fastq files at [`date`] ..."
            COUNT=`ls -1 [DS]RR*.fastq | wc -l`
            if [[ $COUNT > 1 ]] ; then # If only single file, move, don't copy
                cat [DS]RR*.fastq > "$NAME".fastq
                checkFileExists "$NAME".fastq
                rm [DS]RR*
            else
                mv [DS]RR*.fastq "$NAME".fastq
                checkFileExists "$NAME".fastq
                rm [DS]RR*
            fi
            echo "compressing .fastq files at [`date`] ..."
            gzip "$NAME".fastq
            mv "$NAME".fastq.gz $OUTPUT/
        fi
    done
    echo "All done! at [`date`]"
}

download $1
