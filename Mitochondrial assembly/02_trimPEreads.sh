#!/bin/bash


for SP in WAL WAAL STAL LAAL
do

    if [ $SP == "BFAL" ];       then BC="TATCAG"
    elif [ $SP == "LAAL" ];     then BC="GTGAAA"
    elif [ $SP == "STAL" ];     then BC="ATCGTG"
    elif [ $SP == "WAAL" ];     then BC="TAATCG"
    elif [ $SP == "WAL" ];      then BC="GTAGAG"
    fi

        if [ ${SP} != "BFAL" ]; then
        echo -e "\nRun Trimmomatic for 20M data & species : [ $SP ].\n"

        trimmomatic PE \
                -phred33 \
                -threads ${NCORES} \
                -trimlog ${wDIR2}/${SP}_220bp_${BC}.20M.log \
                -summary ${wDIR2}/${SP}_220bp_${BC}.20M.statsSummaryFile.txt \
                -validatePairs ${wDIR2}/${SP}_220bp_${BC}.R1.20M.fastq.gz \
                ${wDIR2}/${SP}_220bp_${BC}.R2.20M.fastq.gz \
                ${wDIR2}/${SP}_220bp_${BC}.R1.20M.paired.fastq.gz \
                ${wDIR2}/${SP}_220bp_${BC}.R1.20M.unpaired.fastq.gz \
                ${wDIR2}/${SP}_220bp_${BC}.R2.20M.paired.fastq.gz \
                ${wDIR2}/${SP}_220bp_${BC}.R2.20M.unpaired.fastq.gz \
                ILLUMINACLIP:TruSeq2-PE_SuppAdapters.fa:2:30:10 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        &> ${wDIR2}/trim_${SP}_20M.out
        fi

done


