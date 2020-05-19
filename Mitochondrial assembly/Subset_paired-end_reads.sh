#!/bin/bash 


for SP in BFAL LAAL STAL WAAL WAL
do

    ((j=j%nparallel)); ((j++==0)) && wait

    #########################################################

    FILE1_frag="${wDIR}/data/${SP}_220bp_*R1.fastq.gz"
    FILE2_frag="${wDIR}/data/${SP}_220bp_*R2.fastq.gz"
    BC_frag=$(echo ${FILE1_frag} | rev | cut -d"/" -f1 | rev | cut -d"_" -f3 | cut -d"." -f1 )

    echo -e "\nA1. Subsampling fragment library at 20M R1 reads for species [ ${SP} ].\n"
    zcat ${FILE1_frag} | seqkit sample -s 200 -n 20000000 -o ${wDIR}/data/${SP}_220bp_${BC_frag}.R1.20M.fastq.gz
    echo -e "\nA2. Subsampling fragment library at 20M R2 reads for species [ ${SP} ].\n"
    zcat ${FILE2_frag} | seqkit sample -s 200 -n 20000000 -o ${wDIR}/data/${SP}_220bp_${BC_frag}.R2.20M.fastq.gz

    #########################################################

    FILE1_jump="${wDIR}/data/${SP}_3kb_*R1.fastq.gz"
    FILE2_jump="${wDIR}/data/${SP}_3kb_*R2.fastq.gz"
    BC_jump=$(echo ${FILE1_jump} | rev | cut -d"/" -f1 | rev | cut -d"_" -f3 | cut -d"." -f1 )

    echo -e "\nA3. Subsampling jumping library at 20M raw R1 reads for species [ ${SP} ].\n"
    zcat ${FILE1_jump} | seqkit sample -s 100 -n 20000000 -o ${wDIR}/data/${SP}_3kb_${BC_jump}.R1.20M.fastq.gz
    echo -e "\nA4. Subsampling jumping library at 20M raw R2 reads for species [ ${SP} ].\n"
    zcat ${FILE2_jump} | seqkit sample -s 100 -n 20000000 -o ${wDIR}/data/${SP}_3kb_${BC_jump}.R2.20M.fastq.gz

    #########################################################

done
wait



