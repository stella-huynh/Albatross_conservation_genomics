for SP1 in EURHE FULGL
do
    for SP2 in BFAL LAAL STAL WAAL WAL
    do

        if [[ ${SP1} != ${SP2} ]]; then

          if [[ $SP1 == "BFAL" ]]; then
            SPECIES1="${iDIR}/${SP1}_genome_v2.fasta"
          elif [[ $SP1 == "LAAL" || $SP1 == "STAL" || $SP1 == "WAAL" || $SP1 == "WAL" ]]; then
            SPECIES1="${iDIR}/${SP1}_genome.fasta"
          elif [[ $SP1 == "GALGA" || $SP1 == "GAVST" || $SP1 == "FULGL" || $SP1 == "EURHE" ]]; then
            SPECIES1="${oDIR}/UCEs_unaligned_${SP1}.fasta"
          fi


          if [[ $SP2 == "BFAL" ]]; then
            SPECIES2="${iDIR}/${SP2}_genome_v2.fasta"
          elif [[ $SP2 == "LAAL" || $SP2 == "STAL" || $SP2 == "WAAL" || $SP2 == "WAL" ]]; then
            SPECIES2="${iDIR}/${SP2}_genome.fasta"
          elif [[ $SP2 == "GALGA" || $SP2 == "GAVST" || $SP2 == "FULGL" || $SP2 == "EURHE" ]]; then
            SPECIES2="${oDIR}/UCEs_unaligned_${SP2}.fasta"
          fi


          OUT="${oDIR}/Megablast_UCEs_${SP1}_${SP2}.txt"
          echo $OUT

          if [[ ! -s ${OUT} ]]; then

             ((j=j%nparallel)); ((j++==0)) && wait
             echo -e "\nA0. Run megablast for species: ${SP1} vs. ${SP2}.\n"

             blastn -task megablast \
                    -query ${SPECIES1} -db ${SPECIES2} \
                    -out ${OUT} -outfmt 6 \
                    -num_threads ${nt} \
                    -evalue 1e-10 -perc_identity 10 -penalty -3 -reward 2 \
                    -gapopen 5 -gapextend 2 -word_size 11 &

                    #-soft_masking true
                    #-F "m S" -s T -use_sw_tback
                    #-word_size (default=28 for megablast, 11 for blastn)
                    #-gapopen (default-0 for megablast, 5 for blastn), -gapextend (default=2)
                    #-reward (default=1 for megablast, 2 for blastn)
                    #-penalty (default=-2 for megablast, 3 for blastn)
                    #-perc_identity (default=0)
          fi


        fi

    done

done

