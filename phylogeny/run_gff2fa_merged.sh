#!/bin/bash

display_help() {
        echo -e "\nUsage: $0 [OPTIONS...] \n"
        echo "    -g/--gene-name      Name of gene to grep from GFF file. [REQUIRED]"
        echo "    -gff/--gff-gile     Name of GFF file. [REQUIRED]"
        echo "    -a/--annotation     Region of gene to extract: mRNA / exon / CDS. [REQUIRED]"
        echo "    -f/-fasta-file      Name of FASTA file to extract gene sequence from. [REQUIRED]"
        echo "    -o/--output-file    Name of OUTPUT file. [OPTIONAL]"
        echo "    -r/--overwrite      1= overwrite output file. 0= append to output file. [ DEFAULT= 1 ]"
        echo "    -app/--append-head  Add name/ID to output FASTA header. [OPTIONAL]"
        echo "    -w/--wrap-seq       Wrap FASTA sequence. [ DEFAULT= 60 ]"
        echo "    -q/--quiet          Mask extra information."
        echo "    -h/--help           Display help."
        echo -e "\n" && exit 1
}


POSITIONAL=()

while [[ $# -gt 0 ]]
do
        key="$1"
case $key in
        -g|--gene-name) #GENE name
        GENE="$2"
        shift # past argument
        shift # past value
        ;;
        -gff|--gff-file) #GFF file to convert
        GFF_FILE="$2"
        shift # past argument
        shift # past value
        ;;
        -a|--annotation) #ANNOTATION name
        ANNOT="$2"
        shift # past argument
        shift # past value
        ;;
        -f|--fasta-file) #FASTA file to convert
        FASTA_FILE="$2"
        shift # past argument
        shift # past value
        ;;
        -o|--output_file) #OUTPUT file name
        OUTFILE="$2"
        shift # past argument
        shift # past value
        ;;
        -r|--overwrite) #OUTPUT file name
        OVERWRITE="$2"
        shift # past argument
        shift # past value
        ;;
        -app|--append-header) #APPEND name/pattern to output FASTASEQ header
        SEQ_HEAD="$2"
        shift # past argument
        shift # past value
        ;;
        -w|--wrap-seq) #wrap sequence in output FASTA file
        WRAP="$2"
        shift # past argument
        shift # past value
        ;;
        -q|--quiet)
        QUIET=1
        shift # past argument
        ;;
        -h|--help)
        display_help #call funtion to display script commands
        exit 0
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



##################################################################


if [[ -z $OVERWRITE ]]; then OVERWRITE=1 ; fi
if [[ -z $WRAP ]]; then WRAP=60 ; fi
if [[ ! -z $SEQ_HEAD ]]; then SEQ_HEAD="${SEQ_HEAD}_" ; fi
if [[ $OVERWRITE == 1 ]]; then rm -f ${OUTFILE} ; fi
if [[ $QUIET == 1 ]]; then Q="--quiet" ; fi


#----------------------------------------------------------------#

if [[ ! $QUIET == 1 ]]; then

  echo -e "\n#---------------------------------------------#"
  echo -e "GENE name = ${GENE}"
  echo -e "ANNOT name = ${ANNOT}"
  echo -e "GFF file = ${GFF_FILE}"
  echo -e "FASTA file = ${FASTA_FILE}"
  echo -e "OUTPUT file = ${OUTFILE}"
  echo -e "OVERWRITE = ${OVERWRITE}"
  echo -e "APPEND to FASTA header = ${SEQ_HEAD}"
  echo -e "#---------------------------------------------#\n"

fi


#----------------------------------------------------------------#


  grep "$GENE" ${GFF_FILE} | grep -w "$ANNOT" \
        | convert2bed --do-not-sort --input=gff < /dev/stdin \
        | seqkit subseq $Q --bed /dev/stdin ${FASTA_FILE} \
        | seqkit replace -p "^(.+):\+ " | seqkit replace -p "^(.+):\- " \
        | seqkit replace -p ":.+$" -r "" \
        | awk '/^>/ {if(prev!=$0) {prev=$0;printf("%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' \
        | seqkit replace -w${WRAP} -p "^" -r "${SEQ_HEAD}" \
  >> ${OUTFILE}



##################################################################


