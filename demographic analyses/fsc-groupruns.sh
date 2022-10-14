#!/bin/bash

wDIR=$1
PREFIX=$2

oDIR="${wDIR}/${PREFIX}/sum_bestlhoods"
mkdir -p -m775 ${oDIR}


MODEL=$( echo $wDIR | rev | cut -d "/" -f1 | rev )


#n=0

#for i in ${wDIR}/${PREFIX}/run*
for i in {1..100}
do

   #N=$( echo $i | rev | cut -d"/" -f1 | rev )
   #cp ${i}/${PREFIX}/*.bestlhoods ${oDIR}/${PREFIX}_${N}.bestlhoods
   cp ${wDIR}/${PREFIX}/run${i}/${PREFIX}/*.bestlhoods ${oDIR}/${PREFIX}_run${i}.bestlhoods


#   if [ $n == 0 ]; then
   if [[ $i == 1 ]]; then
        awk -v model="$MODEL" -v group="$PREFIX" \
        '{OFS="\t"} NR==1 {print "MODEL","GROUP",$0} NR>1 {print model,group,$0}' \
                ${oDIR}/${PREFIX}_run${i}.bestlhoods > ${oDIR}/${PREFIX}_all.bestlhoods
   else
        awk -v model="$MODEL" -v group="$PREFIX" \
        '{OFS="\t"} NR>1 {print model,group,$0}' \
                ${oDIR}/${PREFIX}_run${i}.bestlhoods >> ${oDIR}/${PREFIX}_all.bestlhoods
   fi

#   (( n++ ))

done
