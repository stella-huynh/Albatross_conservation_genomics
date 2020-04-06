#!/bin/bash


POSITIONAL=()

while [[ $# -gt 0 ]]
do
        key="$1"

case $key in
        -i|--infile) #.PED file from PLINK
        INFILE="$2"
        shift # past argument
        shift # past value
        ;;
        -k|--kclusters) #number of K clusters to analyse
        K="$2"
        shift # past argument
        shift # past value
        ;;
        -p|--nparallel) #number of parallel processes
        nparallel="$2"
        shift # past argument
        shift # past value
        ;;
        -b|--bootstrap) #number of bootstraps
        nboots="$2"
        shift # past argument
        shift # past value
        ;;
        -c|--convergence_criteria) #number of bootstraps
        conv="$2"
        shift # past argument
        shift # past value
        ;;
esac

done

set -- "${POSITIONAL[@]}" # restore positional parameters




echo -e "\n###############################################"
echo -e "PARAM files are as following:"
echo -e ".PED file = ${INFILE}"
echo -e "K clusters = ${K}"
echo -e "Number of parallel processes = ${nparallel}"
echo -e "Number of bootstraps to run = ${nboots}"
echo -e "Convergence threshold = ${conv}"
echo -e "###############################################\n"



for i in $(seq 1 1 $K )
do


        ((j=j%nparallel)); ((j++==0)) && wait
        echo -e "\nRun admixture for ${INFILE} : [ $i ].\n"

        admixture -C ${conv} -B${nboots} $INFILE $i &


done

wait

