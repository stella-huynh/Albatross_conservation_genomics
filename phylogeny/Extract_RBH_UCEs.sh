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


for SP in BFAL LAAL STAL WAAL WAL
do
    for OUTG in EURHE FULGL GAVST GALGA
    do 
        awk -v SP=$SP '{FS=OFS="\t"} NR==FNR { outg[$2]=$2; sp[$2]=$1; start[$2]=$7; end[$2]=$8; next } ($1 in outg) { i=index(outg[$1],"_"); s=substr(outg[$1],i+1); print SP"_"s, outg[$1], sp[$1], start[$1], end[$1] }' \
        Megablast_UCEs_${SP}_${OUTG}.besthits.txt \
        Megablast_UCEs_${OUTG}_${SP}.besthits.txt \
        > Megablast_UCEs_${SP}_${OUTG}.ubesthits.txt 
        
        #checking
        diff <(cut -f1 Megablast_UCEs_${SP}_${OUTG}.ubesthits.txt | sed "s/^${SP}_//g") \
             <(cut -f2 Megablast_UCEs_${SP}_${OUTG}.ubesthits.txt | sed "s/${OUTG}_//g")
    done
done







####################################################################

#step 1 : Extract UCEs found in all outgroup for each species

for SP in BFAL LAAL STAL WAAL WAL; 
do 
    awk '{FS=OFS="\t"} NR==FNR {a[$1 FS $3]=$0; next} {idx=$1 FS $3 ; if(idx in a) print a[idx], $0}' Megablast_UCEs_${SP}_EURHE.ubesthits.txt Megablast_UCEs_${SP}_FULGL.ubesthits.txt \
    | awk '{FS=OFS="\t"} NR==FNR {a[$1 FS $3]=$0; next} {idx=$1 FS $3 ; if(idx in a) print a[idx], $0}' - Megablast_UCEs_${SP}_GAVST.ubesthits.txt \
    | awk '{FS=OFS="\t"} NR==FNR {a[$1 FS $3]=$0; next} {idx=$1 FS $3 ; if(idx in a) print a[idx], $0}' - Megablast_UCEs_${SP}_GALGA.ubesthits.txt \
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
  | cut -d"_" -f2- \
  | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {i=index($1,"_"); idx=substr($1,i+1); if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_LAAL_all_overlapEF.ubesthits.txt \
  | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {i=index($1,"_"); idx=substr($1,i+1); if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_STAL_all_overlapEF.ubesthits.txt \
  | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {i=index($1,"_"); idx=substr($1,i+1); if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_WAAL_all_overlapEF.ubesthits.txt \
  | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {i=index($1,"_"); idx=substr($1,i+1); if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_WAL_all_overlapEF.ubesthits.txt \
  > listF_uniq_UCEs.txt

#command below NOT useful anymore, solution found by replacing split(a,$1,"_") by index+substr in awk command.
#sed -i '/^chrun_random$/d' listF_uniq_UCEs.txt


#step 5 : Retrieve the whole region of these final UCEs for each species

for SP in BFAL LAAL STAL WAAL WAL
do 
    awk -v SP=$SP '{FS=OFS="\t"} NR==FNR {scaff[$1]=$0; next} {idx=SP"_"$1; if(idx in scaff) print $0, scaff[idx]}' \
    Megablast_UCEs_${SP}_all_overlapEF.ubesthits.txt listF_uniq_UCEs.txt \
    > Megablast_UCEs_EF_${SP}.txtF
done 


#step 6: extract DNA sequences for each UCEs (one file per gene, all 4 species for each gene)

rDIR="/group/sbs_ssin/stella/albatross/00_genomes/01_WGS/01_dna"
mkdir -p -m775 gene_set

./run_parse_UCEs.sh -i Megablast_UCEs_EF_BFAL.txtF -f ${rDIR}/BFAL_genome_v2.fasta -o gene_set/ -n yes > nohup_parseUCE_BFAL.out &
sleep 20

for SP in LAAL STAL WAAL WAL;
do
    ./run_parse_UCEs.sh \
        -i Megablast_UCEs_EF_${SP}.txtF \
        -f /group/sbs_ssin/stella/albatross/00_genomes/01_WGS/01_dna/${SP}_genome.fasta \
        -o gene_set/ \
        -n no \
        > nohup_parseUCE_${SP}.out &
    sleep 20
done
wait 

  #check all gene files have all species and one sequence per species
  rm -f check_nUCEs.txt
  for file in `ls gene_set/*`; do grep -c "^>" $file >> check_nUCEs.txt ; done
  sort -V check_nUCEs.txt | uniq 
  
  for SP in BFAL LAAL STAL WAAL WAL; do
    rm -f check_UCEs_${SP}.txt check_UCEs_${SP}.err
    for file in `ls gene_set/*`; do
      grep -c "^>${SP}" $file >> check_UCEs_${SP}.txt 
    done
    occ=$(cat check_UCEs_${SP}.txt | paste -sd+ | bc) ; echo "${SP} : ${occ} orthologous sequences extracted."
    if [[ $occ != $(wc -l listF_uniq_UCEs.txt | cut -d" " -f1 ) ]]; then grep -v -n 1 check_UCEs_${SP}.txt > check_UCEs_${SP}.err ; fi #retrieve UCEs (lines) that has failed
  done
  
  
  
#step 7: Run PASTA

nohup ./run_pasta_UCEs.sh -p 5 > nohup_pasta_UCEs.out &


