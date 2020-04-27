#!/bin/bash


# 1a # Blastn Megablast - OUTG to SP

for OUTG in EURHE FULGL GAVST GALGA
do

  # Extract UCEs of 2-3 close outgroups (data from Jarvis et al. 2014, Science)
  ./parse_UCEs.sh ${OUTG} uce-filtered-alignments-without-gator/ ${PWD} > nohup_EURHE.out 2>&1 &
  # Unalign sequences
  seqkit seq -g -w 0 UCEs_${OUTG}.fasta > UCEs_unaligned_${OUTG}.fasta
  # Build database for RBH
  formatdb -p F -i UCEs_unaligned_${SP}.fasta -o F

  # Run reciprocal blast (megablast)
  for SP in BFAL LAAL STAL WAAL WAL
  do
  
    qsub PBS_Rblastn.cmd

    ./run_extract_UCEs_rbhits.sh \
        -r List_UCEs_unaligned_${OUTG}.txt \
        -i Megablast_UCEs_${OUTG}_${SP}.txt \
        -o Megablast_UCEs_${OUTG}_${SP}.besthits.txt 
        > nohup_megaB_${OUTG}_${SP}.txt 2>&1 &

  done
  
  
# 1b # Blastn Megablast - SP to OUTG
  
for SP in BFAL LAAL STAL WAAL WAL
do

  for OUTG in EURHE FULGL GAVST GALGA
  do
  
    qsub PBS_Rblastn.cmd
    
    ./run_extract_UCEs_rbhits.sh \
        -r List_UCEs_unaligned_${OUTG}.txt \
        -i Megablast_UCEs_${SP}_${OUTG}.txt \
        -o Megablast_UCEs_${SP}_${OUTG}.besthits.txt 
        > nohup_megaB_${SP}_${OUTG}.txt 2>&1 &

  done
  
  
  
####################################################################
# Extract UCEs commonly found among species

???????????????
???????????????
???????????????


####################################################################

#step 1 : Extract UCEs found in all outgroup for each species

for SP in BFAL LAAL STAL WAAL WAL; 
do 
    awk '{FS=OFS="\t"} NR==FNR {a[$1 FS $3]=$0; next} {idx=$1 FS $3 ; if(idx in a) print a[idx], $0}' Megablast_UCEs_${SP}_EURHE.ubesthits.txt Megablast_UCEs_${SP}_FULGL.ubesthits.txt \
    | awk '{FS=OFS="\t"} NR==FNR {a[$1 FS $3]=$0; next} {idx=$1 FS $3 ; if(idx in a) print a[idx], $0}' - Megablast_UCEs_${SP}_GAVST.ubesthits.txt \
    awk '{FS=OFS="\t"} NR==FNR {a[$1 FS $3]=$0; next} {idx=$1 FS $3 ; if(idx in a) print a[idx], $0}' - Megablast_UCEs_${SP}_GALGA.ubesthits.txt \
    > Megablast_UCEs_${SP}_all.ubesthits.txt    
done


#step 2 : Extract for each species the whole region that could aligned to any outgroup

for SP in BFAL LAAL STAL WAAL WAL
do
    awk '{FS=OFS="\t"; min=$4; max=$5; for(i=9; i<=NF; i+=5) if($i < min) min=$i; for(j=10; j<=NF; j+=5) if($j > max) max=$j} {print $0, min, max, max-min}' \
    Megablast_UCEs_${SP}_all.ubesthits.txt \
    | cut -f1,3,21-23 \
    > Megablast_UCEs_${SP}_all_fullregion.ubesthits.txt
done


#step 3 : Only keep regions that aligned at least 50% of their full length onto EURHE/FULGL

for SP in BFAL LAAL STAL WAAL WAL
do 
    awk '{FS=OFS="\t"} NR==FNR {a[$1]=$0; tot[$1]=$5; next} {idx=$1; if (idx in a) printf "%s\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n", a[idx], $5-$4, ($5-$4)/tot[idx], $10-$9, ($10-$9)/tot[idx], $15-$14, ($15-$14)/tot[idx], $20-$19, ($20-$19)/tot[idx]}' \
    Megablast_UCEs_${SP}_all_fullregion.ubesthits.txt \
    Megablast_UCEs_${SP}_all.ubesthits.txt \
    | awk '{FS=OFS="\t"} {if($6 >=0.50 && $8 >=0.50) print $0}' \
    > Megablast_UCEs_${SP}_all_overlapEF.ubesthits.txt
done


#step 4 : Extract (final) UCEs that are commonly retained in all species

cut -f1 Megablast_UCEs_BFAL_all_overlapEF.ubesthits.txt \
  | cut -d"_" -f2,3 | \
  awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {split($1,a,"_"); idx=a[2]"_"a[3]; if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_LAAL_all_overlapEF.ubesthits.txt \
  | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {split($1,a,"_"); idx=a[2]"_"a[3]; if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_STAL_all_overlapEF.ubesthits.txt \
  | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {split($1,a,"_"); idx=a[2]"_"a[3]; if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_WAAL_all_overlapEF.ubesthits.txt \
  | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {split($1,a,"_"); idx=a[2]"_"a[3]; if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_WAL_all_overlapEF.ubesthits.txt \
  > listF_uniq_UCEs.txt


#step 5 : Retrieve the whole region of these final UCEs for each species

for SP in BFAL LAAL STAL WAAL WAL
do 
    awk -v SP=$SP '{FS=OFS="\t"} NR==FNR {scaff[$1]=$0; next} {idx=SP"_"$1; if(idx in scaff) print $0, scaff[idx]}' \
    Megablast_UCEs_${SP}_all_overlapEF.ubesthits.txt listF_uniq_UCEs.txt \
    > Megablast_UCEs_EF_${SP}.txtF
done 




  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
