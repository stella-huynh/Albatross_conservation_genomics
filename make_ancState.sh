#!/bin/bash


## For ingroup and outgroup species ##

# step 1: get coverage per nucleotide per site across each species #

angsd -ref $REF \
      -bam $BAMLIST \
      -rf list_autosomes_scafs_forANGSD.txt \
      -out 01_depth/recal/$PREFIX.ACGTpersite \
      -nThreads 20 \
      -nInd $NIND \
      -baq 1 \
      -doCounts 1 -doDepth 1 -dumpCounts 3

# step 2: merge the .pos and .counts file together #

paste <(zcat 01_depth/recal/$PREFIX.ACGTpersite.pos.gz) \
      <(zcat 01_depth/recal/$PREFIX.ACGTpersite.counts.gz) \
      > 01_depth/recal/$PREFIX.ACGTpersite.txt

# step 3: convert coverage into counts of individual per nucleotide per site #

awk -v NIND=$NIND \
    'NR==1 {print $0} 
     NR>1 {printf "%s\t%d\t%d\t%.0f,%.0f,%.0f,%.0f\n", $1, $2, $3, $4/$3*NIND, $5/$3*NIND, $6/$3*NIND, $7/$3*NIND}' \
     01_depth/recal/${SP}_ACGTpersite.txt > 01_depth/recal/${SP}_ACGTpersite.txtF

    # merge the two ingroup species into one if needed #
    awk '{FS=OFS="\t"} FNR==NR {scaff[$1 SUBSEP $2]=$0; next}
                {idx=$1 SUBSEP $2; if(idx in scaff) print scaff[$1 SUBSEP $2], $0}' \
           01_depth/recal/tmp/BFAL_ACGTpersite_${scaff}.tmp 01_depth/recal/tmp/LAAL_ACGTpersite_${scaff}.tmp \
           | awk '{FS=OFS="\t"} {split($4,a,","); split($8,b,",")}
                  {print $1,$2, $3+$7, a[1]+b[1]","a[2]+b[2]","a[3]+b[3]","a[4]+b[4]}' \
           > 01_depth/recal/tmp/BFAL_LAAL_ACGTpersite_${scaff}.tmp

# step 4: split file by scaffolds #

cat 01_depth/recal/BFAL_ACGTpersite.txtF \
    | awk -v PREFIX=$PREFIX '{FS=OFS="\t"} { split($1,a,":"); prevfile=ofile; 
                             ofile="01_depth/recal/tmp/"PREFIX".ACGTpersite_"a[1]".tmp"; 
                             if(NR > 1 && ofile != prevfile) close(prevfile); print $0 >> ofile }' \ 
      01_depth/recal/BFAL_ACGTpersite.txtF &

# step 5: for each scaffold file, aggregate ingroup with outgroup files together #

INGROUP="BFAL_LAAL"
OUTG1="STAL"
OUTG2="WAAL"
OUTG3="WAL"

awk '{FS=OFS="\t"} FNR==NR {scaff[$1 SUBSEP $2]=$0; next} {idx=$1 SUBSEP $2; if(idx in scaff) print scaff[$1 SUBSEP $2], $4}' \
  01_depth/recal/tmp/${INGROUP}.ACGTpersite_${scaff}.tmp \
  01_depth/recal/tmp/${OUTG1}.ACGTpersite_${scaff}.tmp \
  | awk '{FS=OFS="\t"} FNR==NR {scaff[$1 SUBSEP $2]=$0; next} {idx=$1 SUBSEP $2; if(idx in scaff) print scaff[$1 SUBSEP $2], $4}' \
  - 01_depth/recal/tmp/${OUTG2}.ACGTpersite_${scaff}.tmp \
  | awk '{FS=OFS="\t"} FNR==NR {scaff[$1 SUBSEP $2]=$0; next} {idx=$1 SUBSEP $2; if(idx in scaff) print scaff[$1 SUBSEP $2], $4}' \
  - 01_depth/recal/tmp/${OUTG3}.ACGTpersite_${scaff}.tmp \
  > 01_depth/recal/tmp/${INGROUP}_outG.ACGTpersite_${scaff}.tmp &

# step 6: concatenate all scaffold files together #

list=$(for s in `cat list_allscafs.bed | cut -f1`; do (echo "01_depth/recal/tmp/${INGROUP}_outG_ACGTpersite_$s.tmp "); done)
cat $list > 01_depth/recal/${INGROUP}_outG.ACGTpersite.txt
cut -f4,5,6,7 01_depth/recal/${INGROUP}_outG.ACGTpersite.txt > 01_depth/recal/${INGROUP}_outG.ACGTpersite.txtF

# step 7: run est-sfs on splitted files (every 1M lines) #

split -d -a2 -l100000000 01_depth/recal/${INGROUP}_outG.ACGTpersite.txtF \
                         01_depth/recal/${INGROUP}_outG.ACGTpersite.txtF-

cd 01_depth/recal
for file in `ls ${INGROUP}_outG.ACGTpersite.txtF-*`
do
        ((j=j%nparallel)); ((j++==0)) && wait
        est-sfs config-rate6.txt \
                ${file} \
                seedfile.txt \
                ${file}.sfs \
                ${file}.p_anc &
done
wait


# step 8: merge est-sfs outputs with scaffolds' coordinates #

cat ${INGROUP}_outG.ACGTpersite.txtF-*.p_anc > ${INGROUP}_outG.ACGTpersite.p_anc
paste <(cut -f1,2 ${INGROUP}_outG.ACGTpersite.txt) 
      ${INGROUP}_outG.ACGTpersite.p_anc \
      > ${INGROUP}_outG.ACGTpersite.p_anc.txt


# step 9: recode probabilities into allele calls #

awk '{FS=OFS="\t"} 
     $6 >= 0.95 || $7 >= 0.95 || $8 >= 0.95 || $9 >= 0.95 {nc="A"} 
     $10 >= 0.95 || $11 >= 0.95 || $12 >= 0.95 || $13 >= 0.95 {nc="C"} 
     $14 >= 0.95 || $15 >= 0.95 || $16 >= 0.95 || $17 >= 0.95 {nc="G"} 
     $18 >= 0.95 || $19 >= 0.95 || $20 >= 0.95 || $21 >= 0.95 {nc="T"} 
     $6 > 0.05 && $6 < 0.95 || $7 > 0.05 && $7 < 0.95 || $8 > 0.05 && $8 < 0.95|| $9 > 0.05 && $9 < 0.95 {n1="A"} 
     $10 > 0.05 && $10 < 0.95 || $11 > 0.05 && $11 < 0.95 || $12 > 0.05 && $12 < 0.95|| $13 > 0.05 && $13 < 0.95 {n2="C"} 
     $14 > 0.05 && $14 < 0.95 || $15 > 0.05 && $15 < 0.95 || $16 > 0.05 && $16 < 0.95|| $17 > 0.05 && $17 < 0.95 {n3="G"} 
     $18 > 0.05 && $18 < 0.95 || $19 > 0.05 && $19 < 0.95 || $20 > 0.05 && $20 < 0.95|| $21 > 0.05 && $21 < 0.95 {n4="T"} 
     n1 == "A" && n2 == "C" && n3 == "G" {nc="V"; n1=""; n2=""; n3=""} 
     n1 == "A" && n2 == "C" && n4 == "T" {nc="H"; n1=""; n2=""; n4=""} 
     n1 == "A" && n3 == "G" && n4 == "T" {nc="D"; n1=""; n3=""; n4=""} 
     n2 == "C" && n3 == "G" && n4 == "T" {nc="B"; n2=""; n3=""; n4=""} 
     n1 == "A" && n2 == "C" {nc="M"; n1=""; n2=""} 
     n1 == "A" && n3 == "G" {nc="R"; n1=""; n3=""} 
     n1 == "A" && n4 == "T" {nc="W"; n1=""; n4=""} 
     n2 == "C" && n3 == "G" {nc="S"; n2=""; n3=""} 
     n2 == "C" && n4 == "T" {nc="Y"; n2=""; n4=""} 
     n3 == "G" && n4 == "T" {nc="K"; n3=""; n4=""} 
     {print $1, $2, nc}' \
  ${INGROUP}_outG.ACGTpersite.p_anc.txt \
  > ${INGROUP}_outG.ACGTpersite.p_anc.txtF







