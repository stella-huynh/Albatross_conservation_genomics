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

nohup ./run_parse_orthoGenes.sh BFAL 1 gene_set/ > nohup_parseOrtho_BFAL.out 2>&1 &
nohup ./run_parse_orthoGenes.sh LAAL 2 gene_set/ > nohup_parseOrtho_LAAL.out 2>&1 &
nohup ./run_parse_orthoGenes.sh STAL 3 gene_set/ > nohup_parseOrtho_STAL.out 2>&1 &
nohup ./run_parse_orthoGenes.sh WAAL 4 gene_set/ > nohup_parseOrtho_WAAL.out 2>&1 &
nohup ./run_parse_orthoGenes.sh WAL 5 gene_set/ > nohup_parseOrtho_WAL.out 2>&1 &





