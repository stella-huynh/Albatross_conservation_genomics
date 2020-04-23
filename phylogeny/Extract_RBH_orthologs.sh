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
	#create gene folder for local alignment
	mkdir -p -m775 gene_set/Blastallp_ALB_ORTH$i
	
	#reformat starting tree file
	echo "(WAL_ORTH${i},((BFAL_ORTH${i},LAAL_ORTH${i}),(STAL_ORTH${i},WAAL_ORTH${i})))" > gene_set/Blastallp_ALB_ORTH$i/start_tree${i}.newick
	
	#reformat fasta file
	cp gene_set/Blastallp_ALB_ORTH$i.fasta gene_set/Blastallp_ALB_ORTH$i/.
	for SP in BFAL LAAL STAL WAAL WAL; do
		sed -i "s/\(^>${SP}_\).*$/\1ORTH${i}/g" gene_set/Blastallp_ALB_ORTH$i/Blastallp_ALB_ORTH$i.fasta
	done
	
	#run iterative local alignment with PASTA + SATé (using MAFFT as aligner)
	run_pasta.py \
		--input gene_set/Blastallp_ALB_ORTH${i}.fasta \
		--alignment-suffix=.aln \
		--output-directory=gene_set/Blastallp_ALB_ORTH$i \
		--datatype=dna \
		--job gene_set/Blastallp_ALB_ORTH${i} \
		--raxml-search-after \
		--temporaries=/data/huynhs \
		--treefile=gene_set/Blastallp_ALB_ORTH${i}/start_tree${i}.newick \
		--num-cpus=${NCORE} \
		--iter-limit=10 \
		--aligner=mafft \
		--merger=muscle \
		--tree-estimator= &

done













