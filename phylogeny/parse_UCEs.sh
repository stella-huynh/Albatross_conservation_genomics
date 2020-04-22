#!/bin/bash


SP=$1
iDIR=$2
oDIR=$3


rm -f ${oDIR}/UCEs_${SP}.fasta


for file in `ls ${iDIR}/*fasta`
do


        PREFIX=$( echo $file | sed 's/_s.fasta//g' | rev | cut -d"/" -f1 | rev )

        echo -e "\n... Extracting ${SP} for UCE : ${PREFIX} ...\n"

        seqkit grep -n -p "${SP}" $file | seqkit replace --pattern "${SP}" --replacement "${SP}_${PREFIX}" >> ${oDIR}/UCEs_${SP}.fasta

        seqkit seq -g -w 0 UCEs_${SP}.fasta > UCEs_unaligned_${SP}.fasta

done

