#!/bin/bash


POSITIONAL=()

while [[ $# -gt 0 ]]
do
        key="$1"

case $key in
        -i|--bamlist) #list of single-sample BAM files
        BAMLIST="$2"
        shift # past argument
        shift # past value
        ;;
        -r|--region) #number of samples (ie. number of BAM files)
        REGION="$2"
        shift # past argument
        shift # past value
        ;;
        -m|--min_maf) #FASTA file of REFERENCE GENOME
        minMAF="$2"
        shift # past argument
        shift # past value
        ;;
        -s|--snp_pval) #FASTA file of REFERENCE GENOME
        snpPVAL="$2"
        shift # past argument
        shift # past value
        ;;
        -g|--genoCutoff) #PREFIX of OUTPUT file
        gCUTOFF="$2"
        shift # past argument
        shift # past value
        ;;
        -o|--outfile) #chromosomal REGION to process
        OUT="$2"
        shift # past argument
        shift # past value
        ;;
        -q|--outfile2) #chromosomal REGION to process
        OUT2="$2"
        shift # past argument
        shift # past value
        ;;
        -t|--n_threads) #number of THREADS
        NTHREADS="$2"
        shift # past argument
        shift # past value
        ;;
        -h1|--h_snp) #homozyg-snp
        hSNP="$2"
        shift # past argument
        shift # past value
        ;;
        -h2|--h_den) #homozyg-density
        hDEN="$2"
        shift # past argument
        shift # past value
        ;;
        -h3|--h_gap) #homozyg-gap
        hGAP="$2"
        shift # past argument
        shift # past value
        ;;
        -h4|--h_kb) #homozyg-window-threshold
        hKB="$2"
        shift # past argument
        shift # past value
        ;;
        -hw1|--hw_snp) #homozyg-window-snp
        hwSNP="$2"
        shift # past argument
        shift # past value
        ;;
        -hw2|--hw_het) #homozyg-window-heterozygosity
        hwHET="$2"
        shift # past argument
        shift # past value
        ;;
        -hw3|--hw_miss) #homozyg-window-missing
        hwMISS="$2"
        shift # past argument
        shift # past value
        ;;
        -hw4|--hw_thr) #homozyg-window-threshold
        hwTHR="$2"
        shift # past argument
        shift # past value
        ;;
        -pw|--pair_win) #homozyg-window-heterozygosity
        pWIN="$2"
        shift # past argument
        shift # past value
        ;;
        -ps|--pair_step) #homozyg-window-missing
        pSTEP="$2"
        shift # past argument
        shift # past value
        ;;
        -pr|--pair_rval) #homozyg-window-threshold
        pR="$2"
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



echo -e "\n###############################################"
echo -e "PARAM files are as following:"
echo -e "BAMLIST = ${BAMLIST}"
echo -e "REGION = ${REGION}"
echo -e "OUTFILE = ${OUT}"
echo -e "OUTFILE = ${OUT2}"
echo -e "NTHREADS = ${NTHREADS}"
echo -e "SNP param : minMaf = ${minMAF}, SNP_pval = ${snpPVAL}, geno_postCutoff = ${gCUTOFF}."
echo -e "ROH param : h_snp = ${hSNP}, h_den = ${hDEN}, h_gap = ${hGAP}, h_kb = ${hKB} , hw_snp = ${hwSNP}, hw_het = ${hwHET}, hw_miss = ${hwMISS}, hw_thr = ${hwTHR}."
echo -e "LD pruning param : pair_win = ${pWIN}, pair_step = ${pSTEP}, pair_r2val = ${pR}."
echo -e "###############################################\n"



##################################################################

if [ ! -s ${OUT}.tped ] || [ ! -s ${OUT}.tfam ]; then

        echo -e "\nA1. Run ANGSD for ${group}.\n"

  angsd   -bam ${BAMLIST} -rf ${REGION} -out ${OUT} -nThreads ${NTHREADS} \
          -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 0 -minMapQ 20 -minQ 20 -trim 0 \
          -baq 1 -C 50 \
          -doMajorMinor 1 -skipTriallelic 1 \
          -doMaf 1 -SNP_pval ${snpPVAL} \
          -doCounts 1 -minMaf ${minMAF} \
          -doGeno -32 -doPost 1 -postCutoff ${gCUTOFF} \
          -GL 1 \
          -doPlink 2
fi
wait


##################################################################

#if [ ! -s ${OUT2}.hom.summary ]; then

        echo -e "\nA2. Detect ROHs with PLINK for ${group}.\n"

  plink  --tped ${OUT}.tped --tfam ${OUT}.tfam --out ${OUT2} --allow-extra-chr \
         --nonfounders --maf ${minMAF} \
         --homozyg --homozyg-snp ${hSNP} --homozyg-density ${hDEN} --homozyg-gap ${hGAP} \
         --homozyg-window-snp ${hwSNP} --homozyg-window-het ${hwHET} --homozyg-window-missing ${hwMISS} \
         --homozyg-window-threshold ${hwTHR} --homozyg-kb ${hKB}
#fi


##################################################################

OUT3="${OUT}_w${pWIN}s${pSTEP}r${pR}"

if [ ! -s ${OUT3}.prune.in ] || [ ! -s ${OUT3}.prune.out ] ; then

        echo -e "\nA3. Prune-LD for dataset ${group}.\n"

  plink --tped ${OUT}.tped --tfam ${OUT}.tfam --out ${OUT3} --allow-extra-chr \
        --indep-pairwise ${pWIN} ${pSTEP} ${pR}

fi

