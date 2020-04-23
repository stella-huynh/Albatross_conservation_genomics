#!/bin/bash


#==============================================================================================================
### Download UCEs from Avian Phylogenomics Project
#==============================================================================================================
#http://gigadb.org/dataset/101041

tar -xf FASTA_files_of_loci_datasets.tar.gz -C ./
tar -xf uce-filtered-alignments-without-gator.tar.gz -C ./   <--------- this is the 3679 UCEs to be used  

### Reference UCEs species
### 1. Gallus gallus  (chicken)              --> GALGA
### 2. Gavia stellata (red-throated loon)    --> GAVST
### 3. Fulmarus glacialis (northern fulmar)  --> FULGL
### 4. Eurypygia helias (little egret)       --> EURHE


# Extract UCEs for selected outgroups + create NCBI database
for OUTG in GALGA GAVST FULGL EURHE
do
    ./parse_UCEs.sh ${OUTG} uce-filtered-alignments-without-gator ${PWD}  
    makeblastdb -in UCEs_unaligned_${OUTG}.fasta -dbtype nucl -out UCEs_unaligned_${OUTG}.fasta
done


# Generate NCBI database for ingroup species 
for SP in BFAL LAAL STAL WAAL WAL
do
    formatdb -p F -i ${SP}_genome.fasta -o F
done


# Run megablast reciprocally (see bash script)
./run_megablast.sh


# Retrieve the best hit for each UCE and each pairwise megablast
for SP in BFAL LAAL STAL WAAL WAL
do
    for OUTG in GALGA GAVST FULGL EURHE
    do
        ./run_extract_UCEs_rbhits.sh \
           -r List_UCEs_unaligned_${OUTG}.txt \
           -i Megablast_UCEs_${SP}_${OUTG}.txt \
           -o Megablast_UCEs_${SP}_${OUTG}.besthits.txt \
           > nohup_megaB_${SP}_${OUTG}.txt 2>&1 &
           
         ./run_extract_UCEs_rbhits.sh \
           -r List_UCEs_unaligned_${OUTG}.txt \
           -i Megablast_UCEs_${OUTG}_${SP}.txt \
           -o Megablast_UCEs_${OUTG}_${SP}.besthits.txt \
           > nohup_megaB_${OUTG}_${SP}.txt 2>&1 &
    done
done


# Keep only UCEs that passed the filters for each reciprocal blast pair
for SP in BFAL LAAL STAL WAAL WAL
do 
    for OUTG in EURHE FULGL GAVST GALGA
    do 
        awk -v SP=$SP '{FS=OFS="\t"} NR==FNR { outg[$2]=$2; sp[$2]=$1; start[$2]=$7; end[$2]=$8; next } 
                      ($1 in outg) { split(outg[$1],s,"_"); print SP"_"s[2]"_"s[3], outg[$1], sp[$1], start[$1], end[$1] }' \
        Megablast_UCEs_${SP}_${OUTG}.besthits.txt \
        Megablast_UCEs_${OUTG}_${SP}.besthits.txt \
        > Megablast_UCEs_${SP}_${OUTG}.ubesthits.txt      
    done
done


# Check that each UCE (among all OUTG) have expectedly the same hit in the SP (QUERY)
 #check list of UCEs that has found reciprocal best hits (in all 5 queries - TO BE DONE)
 cat *ubesthits.txt | cut -f2 | cut -d"_" -f2,3 | sort -V | uniq > list_uniq_UCEs.txt
 #check query ranges are consistently similar for same UCE across outgroup species
 awk '{FS=OFS="\t"} NR==FNR {outg[$2]=$2; sp[$2]=$3; start[$2]=$4; end[$2]=$5; next} ($2 in outg && $ ) { print }'


# align the UCEs using Phyluce (trimming mode)
phyluce_align_seqcap_align \
--fasta all-taxa-incomplete.fasta \
--output mafft-nexus-edge-trimmed \
--taxa 4 \
--aligner mafft \
--cores 12 \
--incomplete-matrix \
--log-path log






