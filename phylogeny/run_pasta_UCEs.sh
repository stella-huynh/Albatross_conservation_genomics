#!/bin/bash



POSITIONAL=()

while [[ $# -gt 0 ]]
do
        key="$1"

case $key in
        -p|--nparallel) #number of THREADS
        nparallel="$2"
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



if [[ $nparallel == "" || $nparallel != *[[:digit:]]* ]]; then nparallel=1 ; fi


wDIR="/group/sbs_ssin/stella/albatross/05_reseq/initial_pipeline/06_ANGSDv2/05_RBlast/UCEs"
iDIR="${wDIR}/gene_set"

for file in `ls ${iDIR}/*.fasta`
do

        PREFIX=$(echo $file | rev | cut -d"/" -f1 | rev | sed 's/UCE_//' | sed 's/.fasta//')
        oDIR="${iDIR}/UCE_$PREFIX"

        ((j=j%nparallel)); ((j++==0)) && wait


        #create gene folder for local alignment
        mkdir -p -m775 ${oDIR}
        rm -f ${oDIR}/*

        #copy fasta file into folder
        seqkit replace -p "\s.+$" -r "" $file > ${oDIR}/UCE_${PREFIX}.fasta


        #reformat starting tree file
         #copied & modified from preliminary WGA analyses "03_wga/pcactus/seqFile.fa"
        cp ${wDIR}/ALB_startTree.newick ${oDIR}/UCE_${PREFIX}_startTree.newick
        for SP in BFAL LAAL STAL WAAL WAL; do
          sed -i "s/${SP}:/${SP}_${PREFIX}:/g" ${oDIR}/UCE_${PREFIX}_startTree.newick
        done

        #run iterative local alignment with PASTA + SATé (using MAFFT as aligner)
        echo -e "\n====>   Run pasta.py for UCE : ${PREFIX}  <===="

        nohup run_pasta.py \
                --input ${oDIR}/UCE_${PREFIX}.fasta \
                --output-directory ${oDIR} \
                --job ${oDIR} \
                --datatype dna \
                --job UCE_${PREFIX} \
                --raxml-search-after \
                --num-cpus 1 \
                --iter-limit 10 \
                --aligner mafft \
                --merger muscle \
                --tree-estimator fasttree \
                --treefile ${oDIR}/UCE_${PREFIX}_startTree.newick \
                &> ${oDIR}/UCE_${PREFIX}.log &

        sleep 1

done
wait
