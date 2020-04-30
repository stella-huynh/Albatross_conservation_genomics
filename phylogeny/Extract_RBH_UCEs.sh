#!/bin/bash


# 1a # Blastn Megablast - OUTG to SP

for OUTG in EURHE FULGL GAVST GALGA; do

  # Extract UCEs of 2-3 close outgroups (data from Jarvis et al. 2014, Science)
  ./parse_UCEs.sh ${OUTG} uce-filtered-alignments-without-gator/ ${PWD} > nohup_EURHE.out 2>&1 &
  # Unalign sequences
  seqkit seq -g -w 0 UCEs_${OUTG}.fasta > UCEs_unaligned_${OUTG}.fasta
  # Build database for RBH
  formatdb -p F -i UCEs_unaligned_${SP}.fasta -o F

  # Run reciprocal blast (megablast)
  for SP in BFAL LAAL STAL WAAL WAL; do
  
    #qsub PBS_Rblastn.cmd
    ./run_extract_UCEs_rbhits.sh \
        -r List_UCEs_unaligned_${OUTG}.txt \
        -i Megablast_UCEs_${OUTG}_${SP}.txt \
        -o Megablast_UCEs_${OUTG}_${SP}.besthits.txt 
        > nohup_megaB_${OUTG}_${SP}.txt 2>&1 &

  done
done

  
# 1b # Blastn Megablast - SP to OUTG
  
for SP in BFAL LAAL STAL WAAL WAL; do
  for OUTG in EURHE FULGL GAVST GALGA; do
  
    #qsub PBS_Rblastn.cmd
    ./run_extract_UCEs_rbhits.sh \
        -r List_UCEs_unaligned_${OUTG}.txt \
        -i Megablast_UCEs_${SP}_${OUTG}.txt \
        -o Megablast_UCEs_${SP}_${OUTG}.besthits.txt 
        > nohup_megaB_${SP}_${OUTG}.txt 2>&1 &
  done
done 
  
  
####################################################################
# Extract UCEs commonly found among species

for SP in BFAL LAAL STAL WAAL WAL; do
    for OUTG in EURHE FULGL GAVST GALGA; do 
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

## step 1 : Extract UCEs found in all outgroup for each species

for SP in BFAL LAAL STAL WAAL WAL; do 
    awk '{FS=OFS="\t"} NR==FNR {a[$1 FS $3]=$0; next} {idx=$1 FS $3 ; if(idx in a) print a[idx], $0}' Megablast_UCEs_${SP}_EURHE.ubesthits.txt Megablast_UCEs_${SP}_FULGL.ubesthits.txt \
    | awk '{FS=OFS="\t"} NR==FNR {a[$1 FS $3]=$0; next} {idx=$1 FS $3 ; if(idx in a) print a[idx], $0}' - Megablast_UCEs_${SP}_GAVST.ubesthits.txt \
    | awk '{FS=OFS="\t"} NR==FNR {a[$1 FS $3]=$0; next} {idx=$1 FS $3 ; if(idx in a) print a[idx], $0}' - Megablast_UCEs_${SP}_GALGA.ubesthits.txt \
    > Megablast_UCEs_${SP}_all.ubesthits.txt
done


## step 2 : Extract for each species the whole region that could aligned to any outgroup

  #based on all 4 outgroups
  for SP in BFAL LAAL STAL WAAL WAL; do
      awk '{FS=OFS="\t"; min=$4; max=$5; for(i=9; i<=NF; i+=5) if($i < min) min=$i; for(j=10; j<=NF; j+=5) if($j > max) max=$j} {print $0, min, max, max-min}' \
      Megablast_UCEs_${SP}_all.ubesthits.txt \
      | cut -f1,3,21-23 \
      > Megablast_UCEs_${SP}_all_fullregion.ubesthits.txt
  done

  #based on the 2 closest outgroup (EURHE/FULGL)
  for SP in BFAL LAAL STAL WAAL WAL; do
      awk '{FS=OFS="\t"; min=$4; max=$5; for(i=9; i<=10; i+=5) if($i < min) min=$i; for(j=10; j<=NF; j+=5) if($j > max) max=$j} {print $0, min, max, max-min}' \
      Megablast_UCEs_${SP}_all.ubesthits.txt \
      | cut -f1,3,21-23 \
      > Megablast_UCEs_${SP}_all_fullregionEF.ubesthits.txt
  done

## step 3a : Extract intersect region between all outgroups

  #based on all 4 outgroups
  for SP in BFAL LAAL STAL WAAL WAL; do
      bedtools intersect -a <(awk '{OFS="\t"; print $3,$4,$5,$1,$2}' Megablast_UCEs_${SP}_all.ubesthits.txt) \
                         -b <(awk '{OFS="\t"; print $8,$9,$10,$6,$7}' Megablast_UCEs_${SP}_all.ubesthits.txt) \
      | bedtools intersect -a - -b <(awk '{OFS="\t"; print $13,$14,$15,$11,$12}' Megablast_UCEs_${SP}_all.ubesthits.txt) \
      | bedtools intersect -a - -b <(awk '{OFS="\t"; print $18,$19,$20,$16,$17}' Megablast_UCEs_${SP}_all.ubesthits.txt) \
      | awk '{OFS="\t"; print $4"__"$1, $2, $3, $1, $4, $5}' | sort -V \
      | bedtools merge -d 10000 -i - \
      | awk '{OFS="\t"; split($1,a,"__"); print a[1], a[2], $2, $3, $3-$2 }' \
      | awk '{OFS="\t"} NR==FNR {scaff[$1]=$0; next} {if($1 in scaff) print scaff[$1]}' - Megablast_UCEs_${SP}_all_fullregion.ubesthits.txt \
      > Megablast_UCEs_${SP}_all_intersectregion.ubesthits.txt
  done

  #based on the 2 closest outgroup (EURHE/FULGL)
  for SP in BFAL LAAL STAL WAAL WAL; do
      bedtools intersect -a <(awk '{OFS="\t"; print $3,$4,$5,$1,$2}' Megablast_UCEs_${SP}_all.ubesthits.txt) \
                         -b <(awk '{OFS="\t"; print $8,$9,$10,$6,$7}' Megablast_UCEs_${SP}_all.ubesthits.txt) \
      | awk '{OFS="\t"; print $4"__"$1, $2, $3, $1, $4, $5}' | sort -V \
      | bedtools merge -d 10000 -i - \
      | awk '{OFS="\t"; split($1,a,"__"); print a[1], a[2], $2, $3, $3-$2 }' \
      | awk '{OFS="\t"} NR==FNR {scaff[$1]=$0; next} {if($1 in scaff) print scaff[$1]}' - Megablast_UCEs_${SP}_all_fullregionEF.ubesthits.txt \
      > Megablast_UCEs_${SP}_all_intersectregionEF.ubesthits.txt
  done


## step 3b : Only keep regions that aligned at least 50% of their full length onto EURHE($6)/FULGL($8)

  #based on all 4 outgroups
  for SP in BFAL LAAL STAL WAAL WAL; do 
      awk '{FS=OFS="\t"} NR==FNR {a[$1]=$0; tot[$1]=$5; next} {idx=$1; if (idx in a) printf "%s\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n", a[idx], $5-$4, ($5-$4)/tot[idx], $10-$9, ($10-$9)/tot[idx], $15-$14, ($15-$14)/tot[idx], $20-$19, ($20-$19)/tot[idx]}' \
      Megablast_UCEs_${SP}_all_fullregion.ubesthits.txt \
      Megablast_UCEs_${SP}_all.ubesthits.txt \
      | awk '{FS=OFS="\t"} {if($7 >= 0.5 && $9 >= 0.5 && $11 >= 0.5 && $13 >= 0.5) print $0}' \
      | awk '{FS=OFS="\t"} NR==FNR {len[$1]=$5; perc[$1]=$13; next} {if($1 in len) print $0, len[$1], perc[$1]}' \
      - Megablast_UCEs_${SP}_all_intersectregion.ubesthits.txt \
      > Megablast_UCEs_${SP}_all_overlap.ubesthits.txt
  done

  #based on the 2 closest outgroup (EURHE/FULGL)
  for SP in BFAL LAAL STAL WAAL WAL; do
      awk '{FS=OFS="\t"} NR==FNR {scaff[$1]=$5; next} {if($1 in scaff) {gap=$5/scaff[$1]; printf "%s\t%d\t%.2f\n", $0, scaff[$1], gap} }' \
      Megablast_UCEs_${SP}_all_fullregionEF.ubesthits.txt \
      Megablast_UCEs_${SP}_all_intersectregionEF.ubesthits.txt \
      | awk '{FS=OFS="\t"} $7 >= 0.5 {print $0}' \
      > Megablast_UCEs_${SP}_all_overlapEF.ubesthits.txt
  done

## step 4 : Extract (final) UCEs that are commonly retained in all species

  suffix="all_overlapEF.ubesthits.txt" ; outfile="listEF_uniq_UCEs.txt"
  #suffix="all_overlap.ubesthits.txt" ; outfile="listF_uniq_UCEs.txt"

  cut -f1 Megablast_UCEs_BFAL_$suffix \
    | cut -d"_" -f2- \
    | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {i=index($1,"_"); idx=substr($1,i+1); if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_LAAL_$suffix \
    | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {i=index($1,"_"); idx=substr($1,i+1); if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_STAL_$suffix \
    | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {i=index($1,"_"); idx=substr($1,i+1); if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_WAAL_$suffix \
    | awk '{FS=OFS="\t"}  NR==FNR {scaff[$1]=$0; next} {i=index($1,"_"); idx=substr($1,i+1); if(idx in scaff) print scaff[idx]}' - Megablast_UCEs_WAL_$suffix \
    > $outfile

  #command below NOT useful anymore, solution found by replacing split(a,$1,"_") by index+substr in awk command.
  #sed -i '/^chrun_random$/d' listF_uniq_UCEs.txt


#step 5 : Retrieve the whole region of these final UCEs for each species

  for SP in BFAL LAAL STAL WAAL WAL; do 
      awk -v SP=$SP '{FS=OFS="\t"} NR==FNR {scaff[$1]=$0; next} {idx=SP"_"$1; if(idx in scaff) print $0, scaff[idx]}' \
      Megablast_UCEs_${SP}_all_overlapEF.ubesthits.txt listEF_uniq_UCEs.txt \
      > Megablast_UCEs_EF_${SP}.txtF
  done 


#step 6: extract DNA sequences for each UCEs (one file per gene, all 4 species for each gene)

  rDIR="/group/sbs_ssin/stella/albatross/00_genomes/01_WGS/01_dna"
  mkdir -p -m775 gene_set gene_set_EF

  ./run_parse_UCEs.sh -i Megablast_UCEs_EF_BFAL.txtF -f ${rDIR}/BFAL_genome_v2.fasta -o gene_set_EF/ -n yes > nohup_parseUCE_EF_BFAL.out &
  ./run_parse_UCEs.sh -i Megablast_UCEs_F_BFAL.txtF -f ${rDIR}/BFAL_genome_v2.fasta -o gene_set/ -n yes > nohup_parseUCE_F_BFAL.out &
  sleep 20

  for SP in LAAL STAL WAAL WAL ; do
      ./run_parse_UCEs.sh -i Megablast_UCEs_EF_${SP}.txtF \
                          -f ${rDIR}/${SP}_genome.fasta -o gene_set_EF -n no \
                          > nohup_parseUCE_EF_${SP}.out &
      ./run_parse_UCEs.sh -i Megablast_UCEs_F_${SP}.txtF \
                          -f ${rDIR}/${SP}_genome.fasta -o gene_set -n no \
                          > nohup_parseUCE_F_${SP}.out &
      sleep 20
  done
  wait 

  #check all gene files have all species and one sequence per species
  rm -f check_nUCEs.txt check_nUCEs_EF.txt
  for file in `ls gene_set/*`; do grep -c "^>" $file >> check_nUCEs.txt ; done
  for file in `ls gene_set_EF/*`; do grep -c "^>" $file >> check_nUCEs_EF.txt ; done
  sort -V check_nUCEs.txt | uniq  ;  sort -V check_nUCEs_EF.txt | uniq
  
  for SP in BFAL LAAL STAL WAAL WAL; do
    rm -f check_UCEs_${SP}.txt check_UCEs_${SP}.err
    for file in `ls gene_set/*.fasta`; do
      grep -c "^>${SP}" $file >> check_UCEs_${SP}.txt 
    done
    occ=$(cat check_UCEs_${SP}.txt | paste -sd+ | bc) ; echo "${SP} : ${occ} orthologous F sequences extracted."
    if [[ $occ != $(wc -l listF_uniq_UCEs.txt | cut -d" " -f1 ) ]]; then grep -v -n 1 check_UCEs_${SP}.txt > check_UCEs_${SP}.err ; fi #retrieve UCEs (lines) that has failed
    
    rm -f check_UCEs_EF_${SP}.txt check_UCEs_EF_${SP}.err
    for file in `ls gene_set_EF/*.fasta`; do
      grep -c "^>${SP}" $file >> check_UCEs_EF_${SP}.txt 
    done
    occ=$(cat check_UCEs_EF_${SP}.txt | paste -sd+ | bc) ; echo "${SP} : ${occ} orthologous EF sequences extracted."
    if [[ $occ != $(wc -l listEF_uniq_UCEs.txt | cut -d" " -f1 ) ]]; then grep -v -n 1 check_UCEs_EF_${SP}.txt > check_UCEs_EF_${SP}.err ; fi #retrieve UCEs (lines) that has failed
  done
  
  
  
#step 7: Run PASTA + trimAl

  nohup ./run_pasta_UCEs.sh -p 5 -i gene_set/ > nohup_pasta_F_UCEs.out && wait
  find gene_set/ -not -empty -type f -name "*.marker001.*" -ls | wc -l #check all alignment have been generated correctly  
  nohup ./run_pasta_UCEs.sh -p 5 -i gene_set_EF/ > nohup_pasta_EF_UCEs.out && wait
  find gene_set_EF/ -not -empty -type f -name "*.marker001.*" -ls | wc -l #check all alignment have been generated correctly  

  list=$(ls gene_set/*.fasta) ; dir="gene_set"
  list=$(ls gene_set_EF/*.fasta) ; dir="gene_set"

  for file in $list ; do
   PREFIX=$(echo $file | rev | cut -d"/" -f1 | rev | sed 's/.fasta//')
   oDIR="${dir}/$PREFIX"
   trimal  -in ${oDIR}/${PREFIX}.marker001.${PREFIX}.aln \
            -out ${oDIR}/${PREFIX}.trim.aln \
            -htmlout ${oDIR}/${PREFIX}.trim.html \
            -automated1 &
   sleep 2
  done
  
  find gene_set/ -not -empty -type f -name "*trim.aln" -ls | wc -l      #check all trimmed alignment have been generated correctly
  find gene_set_EF/ -not -empty -type f -name "*trim.aln" -ls | wc -l   #check all trimmed alignment have been generated correctly

  #copy all ALN files to a new directory and simpler names
  mkdir -p -m775 alnF alnF_trim
  for file in $(find gene_set/ -type f -name "*trim.aln"); do
    seqkit replace -p "_ch.+$" -r "" $file > alnF_trim/$(echo $file | rev | cut -d"/" -f1 | rev); done  #find gene_set/ -type f -name "*trim.aln" -exec cp {} alnF_trim \;
  for file in $(find gene_set/ -type f -name "*.marker001*.aln") ; do
    seqkit replace -p "_ch.+$" -r "" $file > alnF/$(echo $file | rev | cut -d"." -f1-2 | rev); done     #echo $file | sed -e "p;s/.*\/.*\/.*UCE/alnF\/UCE/g" | xargs -n2 cp ; done
  mv alnF_trim gene_set/ ; mv alnF gene_set/

  #copy all ALN files to a new directory and simpler names
  mkdir -p -m775 alnEF alnEF_trim
  for file in $(find gene_set_EF/ -type f -name "*trim.aln"); do
    seqkit replace -p "_ch.+$" -r "" $file > alnEF_trim/$(echo $file | rev | cut -d"/" -f1 | rev); done #find gene_set_EF/ -type f -name "*trim.aln" -exec cp {} alnEF_trim \;
  for file in $(find gene_set_EF/ -type f -name "*.marker001*.aln") ; do
    seqkit replace -p "_ch.+$" -r "" $file > alnEF/$(echo $file | rev | cut -d"." -f1-2 | rev); done    #echo $file | sed -e "p;s/.*\/.*\/.*UCE/alnEF\/UCE/g" | xargs -n2 cp ; done
  mv alnEF_trim gene_set_EF/ ; mv alnEF gene_set_EF/


#step 8 : run gene trees (RAxML-ng) using ParGenes

module load mpich/gcc/3.1.4
module load gcc/8.2.0


pargenes.py -a gene_set_EF/aln_EF/ -o gene_set_EF/aln_EF/parGenes/ \
            -d nt -c ${NTHREADS} -m \
            --per-msa-modeltest-parameters --modeltest-criteria "AICc" --modeltest-perjob-cores 4 \

for file in `ls gene_set/aln_EF/*.aln`; do
  PREFIX=$(echo $file | rev | cut -d"/" -f1 | rev | sed 's/.aln//')
  oDIR="gene_set/aln_EF"

  #pargenes.py -a ${oDIR}/${PREFIX}.trim.aln -o ${oDIR}/${PREFIX}.trim.tree \
  #            -d nt -c ${NTHREADS} -R "--model GTR" -b 100 --use-astral \
  #            &> ${oDIR}/${PREFIX}.pargenes.log
  #export relevant data
  export.py -i ${oDIR}/${PREFIX}.trim.tree -o ${oDIR}/${PREFIX}.trim.tree.export \
            -best-ml-tree --best-ml-model --bootstrap-trees --support-values-tree \
            &> ${oDIR}/${PREFIX}.pargexport.log
  
done


#step 9 : run species tree (ASTRAL)

find gene_set/ -name "*.bestTree" -type -f -exec -c sh 'cat "{}" >> gene_set/gene_trees.trim.bestTrees' ';' 

#run species tree
java -jar astral.5.7.3.jar -i gene_set/gene_trees.trim.bestTrees \
                           -o gene_set/gene_trees.trim.species.tree \
                           -t 1 \
                           2> gene_set/gene_trees.trim.species.log &
wait

#estimate node supports
java -jar astral.5.7.3.jar -q gene_set/gene_trees.trim.species.tree \
                           -i gene_set/gene_trees.trim.bestTrees \
                           -o gene_set/gene_trees.trim.species2.tree
                           -t 2 \
                           2> gene_set/gene_trees.trim.species2.log &
wait



#step 10: Do bootstrap analyses
#ASTRAL recommand local posterior probability and quadripartition (the four clusters around a branch)

#prepare bootstrap file
???????
gene_set/gene_trees.trim.bootList

#bootstrap with ASTRAL
#* As ASTRAL performs bootstrapping, it continually outputs the bootstrapped ASTRAL tree for each replicate. So, if the number of replicates is set to 100, it first outputs 100 trees. Then, it outputs a greedy consensus of all the 100 bootstrapped trees (with support drawn on branches). Finally, it performs the main analysis (i.e., on trees provided using -i option) and draws branch support on this main tree using the bootstrap replicates. Therefore, in this example, the output file will include 102 trees.
#* What to use: The most important tree is the tree outputted at the end; this is the ASTRAL tree on main input trees, with support values drawn based on bootstrap replicates. Our analyses have shown this tree to be better than the consensus tree in most cases.

java -jar astral.5.7.3.jar -i gene_set/gene_trees.trim.bestTrees \
                           -b gene_set/gene_trees.trim.bootList \
                           -o gene_set/gene_trees.trim.boot.species.tree \
                           -r 100 &
#Visualize in DensiTree


