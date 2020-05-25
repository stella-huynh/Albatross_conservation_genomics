################################################################################
# Run BLASTP (using blastall or blast+)

################################################################################
# Extract unique best hit of each query (SP1) onto database (SP2)

nohup bash run_extract_rbhits.sh list_proteins_${SP1}.txt Blastallp_${SP1}_${SP2}.txt Blastallp_${SP1}_${SP2}.besthits.txt &



################################################################################
# Keep query/database pairs with identical reciprocal best hit


SP1="BFAL"; SP2="STAL"
awk -F'[ \t]' '{OFS="\t"} NR==FNR {a[$1 FS $2]=$1; b[$1 FS $2]=$2 ; next} {idx=$2 FS $1; if (idx in a) print a[idx], b[idx], $2, $1}' \
	Blastallp_${SP1}_${SP2}.besthits.txt \
	Blastallp_${SP2}_${SP1}.besthits.txt \
	> Blastallp_${SP1}_${SP2}.ubesthits.txt



################################################################################
# Check if hits are indeed the same between reciprocal blasts


for PAIR in BFAL_LAAL BFAL_STAL BFAL_WAAL BFAL_WAL LAAL_STAL LAAL_WAAL LAAL_WAL STAL_WAAL STAL_WAL WAAL_WAL 
do 
	file="Blastallp_${PAIR}.ubesthits.txt"
	check=$(diff <(awk -F'[ \t]' '{print $1,$2,$13,$14}' $file | cut -d" " -f1) <(awk -F'[ \t]' '{print $1,$2,$13,$14}' $file | cut -d" " -f4))
	if [[ $check = "" ]]; then echo "ok"; else echo "not ok" ; fi
done

for PAIR in BFAL_LAAL BFAL_STAL BFAL_WAAL BFAL_WAL LAAL_STAL LAAL_WAAL LAAL_WAL STAL_WAAL STAL_WAL WAAL_WAL 
do 
	file="Blastallp_${PAIR}.ubesthits.txt"
	check=$(diff <(awk -F'[ \t]' '{print $1,$2,$13,$14}' $file | cut -d" " -f2) <(awk -F'[ \t]' '{print $1,$2,$13,$14}' $file | cut -d" " -f3))
	if [[ $check = "" ]]; then echo "ok"; else echo "not ok" ; fi 
done


################################################################################
# Keep reciprocal best hits that are found between all 5 species

awk '{FS=OFS="\t"} NR==FNR {a[$2 FS $3]=$0; next} { idx=$1 FS $2 ; if(idx in a) print a[idx]}' Blastallp_BFAL_all_ubesthits.txt <(sed 's/ /\t/g' Blastallp_LAAL_STAL.ubesthits.txt) | awk '{FS=OFS="\t"} NR==FNR {a[$2 FS $4]=$0; next} { idx=$1 FS $2 ; if(idx in a) print a[idx]}' - <(sed 's/ /\t/g' Blastallp_LAAL_WAAL.ubesthits.txt) | awk '{FS=OFS="\t"} NR==FNR {a[$2 FS $5]=$0; next} { idx=$1 FS $2 ; if(idx in a) print a[idx]}' - <(sed 's/ /\t/g' Blastallp_LAAL_WAL.ubesthits.txt) | awk '{FS=OFS="\t"} NR==FNR {a[$3 FS $4]=$0; next} { idx=$1 FS $2 ; if(idx in a) print a[idx]}' - <(sed 's/ /\t/g' Blastallp_STAL_WAAL.ubesthits.txt) | awk '{FS=OFS="\t"} NR==FNR {a[$3 FS $5]=$0; next} { idx=$1 FS $2 ; if(idx in a) print a[idx]}' - <(sed 's/ /\t/g' Blastallp_STAL_WAL.ubesthits.txt) | awk '{FS=OFS="\t"} NR==FNR {a[$4 FS $5]=$0; next} { idx=$1 FS $2 ; if(idx in a) print a[idx]}' - <(sed 's/ /\t/g' Blastallp_WAAL_WAL.ubesthits.txt) > Blastallp_ALB.ubesthist.txt

#GFF markers type
  #maker-
  #exonerate_protein2genome-
  #snap-
  #augustus-
  #augustus_masked-
#Extract and rename orthologs in GFF3 file
#seqkit grep -r -p $LOC --id-regexp "^(\S\s\S)" tmp.fasta | seqkit replace -p "^.+\s" -r "${SP}_" | head -n1


################################################################################
# Extract corresponding sequence data and separate by orthologous gene (ie. one fasta file per ortholog)

nohup ./run_parse_orthoGenes.sh BFAL 1 gene_set/ > nohup_parseOrtho_BFAL.out 2>&1 &
nohup ./run_parse_orthoGenes.sh LAAL 2 gene_set/ > nohup_parseOrtho_LAAL.out 2>&1 &
nohup ./run_parse_orthoGenes.sh STAL 3 gene_set/ > nohup_parseOrtho_STAL.out 2>&1 &
nohup ./run_parse_orthoGenes.sh WAAL 4 gene_set/ > nohup_parseOrtho_WAAL.out 2>&1 &
nohup ./run_parse_orthoGenes.sh WAL 5 gene_set/ > nohup_parseOrtho_WAL.out 2>&1 &

#double-check one sequence of each species has effectively been retrieved in the output 
 #(12619 orthologs retained, should find 5 FASTA headers in each gene file)
 #for i in {1..12619}; do grep -c "^>" 05_RBlast/orthologs/multiple_hits/gene_set/Blastallp_ALB_ORTH${i}.fasta >> tmp.txt ; done
 col=0
 for SP in BFAL LAAL STAL WAAL WAL; do
 
   # alt 1 # check all species are present once in each ortholog file
   for i in {1..12619}; do
	grep -c "^>${SP}" 05_RBlast/orthologs/multiple_hits/gene_set/Blastallp_ALB_ORTH${i}.fasta >> check_orth_${SP}.txt
   done
   occ=$(cat check_orth_${SP}.txt | paste -sd+ | bc) ; echo "${SP} : ${occ} orthologous sequences extracted."
   if [[ $occ != 12619 ]]; then grep -v -n 1 check_orth_${SP}.txt > check_orth_${SP}.err ; fi #retrieve $LOC lines that has failed
   
   # alt 2 # check FASTA header of species sequence is as expected from list "Blastallp_ALB.ubesthits.txt"
   for i in {1..12619}; do
	grep "^>${SP}" 05_RBlast/orthologs/multiple_hits/gene_set/Blastallp_ALB_ORTH${i}.fasta >> checkseq_orth_${SP}.txt
   done
   col=(( col +1 ))
   checkdiff=$( diff <(cut -f$col Blastallp_ALB.ubesthits.txt) <(sed "s/>${SP}_//g" checkseq_orth_${SP}.txt) )
   if [[ -z $checkdiff ]]; then echo "No differences observed with list of unique best hits for species : \"${SP}\"."
   else echo "WARNINGS! : differences are observed with list of unique best hits for species : \"${SP}\"." ; fi
   
   
 done


################################################################################
# Locally align sequences per orthologous gene

for i in {1..12619}
do
	oDIR="exon_set/Blastallp_ALB_ORTH$i" && mkdir -p -m775 ${oDIR}
	
	#locally align orthologous sequences using MAFFT
	mafft --nuc --adjustdirection --auto --op 2 --thread 4 \
	      --progress ${FULL_PATH}/${oDIR}/Blastallp_ALB_ORTH$i.mafft.fas \
	      ${oDIR}/Blastallp_ALB_ORTH$i.fasta \
	      > ${oDIR}/Blastallp_ALB_ORTH$i.mafft.fas && wait
	
	#trim alignment from problematic gaps
	trimAl -in ${oDIR}/Blastallp_ALB_ORTH$i.mafft.fas \
	       -out ${oDIR}/Blastallp_ALB_ORTH$i.mafft.trim.aln \
	       -htmlout ${oDIR}/Blastallp_ALB_ORTH$i.mafft.trim.html \
	       -automated1 && wait
done


################################################################################
# Check alignments completeness and remove the unresolved/small ones

#get alignment info on %id
  for file in `ls 05_alignments/exon_set/alnF/*`; do 
  	PREFIX=$(echo $file | rev | cut -d"." -f2- | rev ) ; infoalign $file $PREFIX.infoalign -only -heading -name -alignlength -idcount -simcount -diffcount -change 2>/dev/null &
  for file in `ls 05_alignments/exon_set/alnF_trim/*`; do 
  	PREFIX=$(echo $file | rev | cut -d"." -f2- | rev ) ; infoalign $file $PREFIX.infoalign -only -heading -name -alignlength -idcount -simcount -diffcount -change 2>/dev/null &

#cp files in new folders for ParGenes
  for file in $(find 05_alignments/exon_set/ -type f -name "*.mafft.fas") ; do 
	seqkit replace -p "_ORTH.+$" -r "" $file | sed 's/_R_//g' > 05_alignments/exon_set/alnF/$(echo $file | rev | cut -d"/" -f1 | rev); done &
  for file in $(find 05_alignments/exon_set/ -type f -name "*.mafft.trim.aln") ; do 
	seqkit replace -p "_ORTH.+$" -r "" $file | sed 's/_R_//g' > 05_alignments/exon_set/alnF_trim/$(echo $file | rev | cut -d"/" -f1 | rev); done &

#check alignment completeness
  for file in `ls 05_alignments/exon_set/alnF/*.mafft.fas` ; do 
  	seqkit fx2tab -ln $file >> check_exons_alnF.txt ; 
	nline=$(seqkit fx2tab -ln $file) ; if [[ $nline < 5 ]]; then echo $file >> check_exons_alnF2.txt ; fi ; done &
  for file in `ls 05_alignments/exon_set/alnF_trim/*.mafft.trim.aln` ; do 
  	seqkit fx2tab -ln $file >> check_exons_alnFtrim.txt 
  	nline=$(seqkit fx2tab -ln $file) ; if [[ $nline < 5 ]]; then echo $file >> check_exons_alnFtrim2.txt ; fi ; done &
wait
#Remove bad/poor alignments


################################################################################
# Infer species tree phylogeny using RAxML+ASTRAL (in parGenes) 

pargenes.py -a exon_set/alnF -o exon_set/alnF/parGenes \
	    -m --modeltest-criteria "AICc" --modeltest-perjob-cores 4 --use-astral -b 100 -d nt -c 5 \
	    &> nohup_exons_parGenes_alnF.out &

pargenes.py -a exon_set/alnF_trim -o exon_set/alnF_trim/parGenes \
	    -m --modeltest-criteria "AICc" --modeltest-perjob-cores 4 --use-astral -b 100 -d nt -c 5 \
	    &> nohup_exons_parGenes_alnFtrim.out &
wait







