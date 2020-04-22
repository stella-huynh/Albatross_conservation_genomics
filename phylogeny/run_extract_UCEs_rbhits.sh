#!/bin/bash



POSITIONAL=()

while [[ $# -gt 0 ]]
do
        key="$1"

case $key in
        -r|--REF) #List of UCEs sequences in outgroup
        LIST="$2"
        shift # past argument
        shift # past value
        ;;
        -i|--infile) #Output from Megablast
        INFILE="$2"
        shift # past argument
        shift # past value
        ;;
        -o|--outfile) #New output filename
        OUTFILE="$2"
        shift # past argument
        shift # past value
        ;;
esac

done

set -- "${POSITIONAL[@]}" # restore positional parameters



#filename="/group/sbs_ssin/stella/albatross/05_reseq/initial_pipeline/06_ANGSDv2/05_RBlast/UCEs/List_UCEs_unaligned_${OUTG}.txt"

#INFILE="${wDIR}/Megablast_${SP1}_${SP2}.txt"
#OUTFILE="${wDIR}/Blastallp_${SP1}_${SP2}.besthits.txt"


echo -e "############################################\n"
echo -e "Running script \"run_extract_UCEs_rbhits.sh\" with following parameters: \nLIST = $LIST \nINFILE = $INFILE \nOUTFILE = $OUTFILE \n"
echo -e "\n############################################\n\n\n"



rm -f $OUTFILE


while IFS= read -r line
do


        PROT=$( echo $line | cut -d" " -f1 )
        LEN=$( echo $line | rev | cut -d" " -f1 | rev )


        echo "----  Searching $PROT  ----"
        # extract best hit for OUTG

        # if Megablast_${OUTG}_${SP}
        LOC1=$( cat ${INFILE} | grep -v "^#" | grep -w "^$PROT" | sort -k12rn | head -n 1 )
        # check corresponding hit in SP2 is unique or find its best hit in same LOC1
        L=$(echo $LOC1 | cut -d" " -f2); S=$(echo $LOC1 | cut -d" " -f9); E=$(echo $LOC1 | cut -d" " -f10)
        LOC2=$( awk -v loc=$L -v start=$S -v end=$E -F'\t' '$2 == loc && $9 >= start && $10 <= end {print $0}' ${INFILE} | sort -k12rn | head -n 1 )

        # if Megablast_${SP}_${OUTG}
        if [ -z $L ]; then
         LOC1=$( cat ${INFILE} | grep -v "^#"| grep -w "$PROT" | sort -k12rn | head -n 1 )
         # check corresponding hit in SP2 is unique or find its best hit in same LOC1
         L=$(echo $LOC1 | cut -d" " -f1); S=$(echo $LOC1 | cut -d" " -f7); E=$(echo $LOC1 | cut -d" " -f8)
         LOC2=$( awk -v loc=$L -v start=$S -v end=$E -F'\t' '$1 == loc && $7 >= start && $8 <= end {print $0}' ${INFILE} | sort -k12rn | head -n 1 )
        fi


        # write UCE pair to output file
        if [[ $LOC1==$LOC2 && ! -z $LOC1 && ! -z $LOC2 ]]; then

          echo "----  Writing $PROT to the output  ----"
          echo $LOC1 | sed 's/ /\t/g' >> ${OUTFILE}

        fi


done < $LIST

