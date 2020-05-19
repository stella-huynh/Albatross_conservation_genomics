#!/bin/bash


for SP in BFAL LAAL STAL WAAL WAL
do

    echo -e "\nStep 1: Initial mapping assembly using MIRA 4.0\n"
    mira manifest_GRCg6a_${SP}.conf >& log_mapping_GRCg6a-${SP}.out &

done



for SP in BFAL LAAL STAL WAAL WAL
do

    echo -e "\nStep 2: Baiting and iterative mapping using MITObim 1.7\n"
    
    if [[ $SP == "BFAL" ]]; then BC="TATCAG"
    elif [[ $SP == "LAAL" ]]; then BC="GTGAAA"
    elif [[ $SP == "STAL" ]]; then BC="ATCGTG"
    elif [[ $SP == "WAAL" ]]; then BC="TAATCG"
    elif [[ $SP == "WAL" ]]; then BC="GTAGAG"
    fi
    
    interleave-fastqgz-MITOBIM.py \
            ${SP}_220bp_${BC}.R1.20M.paired.fastq.gz
            ${SP}_220bp_${BC}.R2.20M.paired.fastq.gz
            > ${SP}_220bp_${BC}.20M.paired.fastq.gz
    
    MITObim.pl -start 1 -end 10 --clean \
               -sample ${SP} \
               -ref GRCg6a_genome \
               -readpool ${wDIR}/data/${SP}_220bp_${BC}.20M.paired.fastq.gz \
               -maf ${wDIR}/MIRA-GRCg6a-${SP}_assembly/MIRA-GRCg6a-${SP}_d_results/MIRA-GRCg6a-${SP}_out.maf \
               --redirect_tmp /data/huynhs \
               --mirapath /group/sbs_ssin/envs/mira_v4.0.2/bin/ \
               &> ${oDIR}/MITObim_run_GRCg6a-${SP}.log

    #redirect temporary files to unmounted disk
    #use MIRA 4.0 for compatibility with MITObim 1.7

done
