#!/bin/bash


POSITIONAL=()

while [[ $# -gt 0 ]]
do
        key="$1"

case $key in
        -f|--synteny_file) #full path of Satsuma synteny output file
        SYNFILE="$2"
        shift # past argument
        shift # past value
        ;;
        -i|--input_file) #full path of file with QUERY sequences to reorder along REF chromosomes
        INFILE="$2"
        shift # past argument
        shift # past value
        ;;
        -o|--output_file) #full path & name for output file
        OUTFILE="$2"
        shift # past argument
        shift # past value
        ;;
        -s|--species) #species name
        SPECIES="$2"
        shift # past argument
        shift # past value
        ;;
        -c|--column_scaff) #column in input file listing the QUERY sequences
        COL="$2"
        shift # past argument
        shift # past value
        ;;
        --default)
        DEFAULT=YES
        shift # past argument
        ;;
        *) # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;
esac

done

set -- "${POSITIONAL[@]}" # restore positional parameters


sDIR=$( echo ${SYNFILE} | rev | cut -d"/" -f 2- | rev ) #path directory of Satsuma file
sPREFIX="${sDIR}/satsuma_summary_${SPECIES}"



#1##  Format REF name & Extract lengths of each REF/QRY hits  ###

sed 's/_dna.*REF//' ${SYNFILE} | sed 's/^/GRCg6a_/' | awk 'BEGIN {FS=OFS="\t"} {print $0, $3-$2}' > ${sPREFIX}_len.txt
#columns: REF / Rstart / Rend / QRY / Qstart / Qend / QC / FRM / LEN , (see Satsuma manual)



#2##  Sum hit lengths for each REF/QRY pairs  ###

awk '{FS=OFS="\t"; a[$1"\t"$4] += $9; qc[$1"\t"$4] += $7; count[$1"\t"$4]++; frame[$1"\t"$4] = $8; min[$1"\t"$4] = $2; max[$1"\t"$4] = $3} { min[$1"\t"$4]=$2 < min[$1"\t"$4] ? $2 : min[$1"\t"$4]; max[$1"\t"$4]=$3 > max[$1"\t"$4] ? $3:max[$1"\t"$4] } END {for (i in a) {print i, a[i], qc[i]/count[i], frame[i], min[i], max[i] } }' ${sPREFIX}_len.txt | sort -k2.11n > ${sPREFIX}_sumlen.txt
#columns: REF / QRY / SUMLEN / meanQC / FRM / min_Rstart / max_Rend



#3##  Keep best hits (ie. total longest hit onto REF) for each QRY scaffold/chromosome  ###

awk 'BEGIN {FS=OFS="\t"} $3 > a[$2] {a[$2]=$3; b[$2]=$1} END {for(k in a) print b[k], k, a[k]}' ${sPREFIX}_sumlen.txt | sort -k2.11n > ${sPREFIX}_umatch.txt
#columns: REF / QRY / SUMLEN



#4##  Extract & Sort QRY scaffolds based on their hits along REF chromosomes  ###

rm -rf ${sPREFIX}_sortumatch.txt

while read REF QRY LEN
do
        awk -v REF="$REF" -v QRY="$QRY" \
        'BEGIN {FS=OFS="\t"} $1==REF && $4==QRY {print $0}' ${sPREFIX}_len.txt \
        | awk 'BEGIN {FS=OFS="\t"; getline; ref=$1; qry=$4; min=$2; max=$3; frm=$8} {(min>$2) ? min=$2:""; (max>$3) ?"":max=$3} END {print ref, qry, min, max, frm}' \
        | sort -V -k1.8 -k3n
done < ${sPREFIX}_umatch.txt >> ${sPREFIX}_sortumatch.txt
#columns: REF / QRY / Rstart / Rend / FRM
# The main table is now ready for use to rearrange QRY scaffolds in other files



#5##  Reorder INPUT file according to Satsuma synteny order of QRY onto REF  ###

awk -v COL="$COL" 'NR==FNR { FS=OFS="\t"; REF[$2]=$1; Rstart[$2]=$3; next } { print REF[$COL], Rstart[$COL], $0 }' \
		${sPREFIX}_sortumatch.txt \
		${INFILE} \
		| sort -V -k1.8 -k2n \
		> ${OUTFILE}
#Two columns are appended at the beginning of the INPUT file, with REF chromosome in 1st column and corresponding REF position in 2nd column.



