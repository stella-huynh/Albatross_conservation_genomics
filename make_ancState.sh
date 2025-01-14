#!/bin/bash

nparallel=2

## For ingroup and outgroup species ##

# step 1: get coverage per nucleotide per site across each species #

      angsd -ref $REF \
            -bam $BAMLIST \
            -rf list_allscafs_forANGSD.txt \
            -out 01_depth/recal/$PREFIX.ACGTpersite \
            -nThreads 20 \
            -nInd $NIND \
            -doCounts 1 -doDepth 1 -dumpCounts 3

# step 2: merge the .pos and .counts file together #

      paste <(zcat 01_depth/recal/$PREFIX.ACGTpersite.pos.gz) \
            <(zcat 01_depth/recal/$PREFIX.ACGTpersite.counts.gz) \
            > 01_depth/recal/$PREFIX.ACGTpersite.txt

# step 3: convert coverage into counts of individual per nucleotide per site #

      awk -v NIND=$NIND \
          'NR==1 {print $0} 
           NR>1 {printf "%s\t%d\t%d\t%.0f,%.0f,%.0f,%.0f\n", $1, $2, $3, $4/$3*NIND, $5/$3*NIND, $6/$3*NIND, $7/$3*NIND}' \
           01_depth/recal/${SP}.ACGTpersite.txt \
           > 01_depth/recal/${SP}.ACGTpersite.txtF

# step 4: split file by scaffolds #

      rm -f 01_depth/recal/tmp/*.tmp #to avoid duplicating tmp file contents
      awk -v SP=$SP '{FS=OFS="\t"} { split($1,a,":"); prevfile=ofile; 
                     ofile="01_depth/recal/tmp/"SP".ACGTpersite_"a[1]".tmp"; 
                     if(NR > 1 && ofile != prevfile) close(prevfile); print $0 >> ofile }' \ 
          01_depth/recal/${SP}.ACGTpersite.txtF &

          # merge ingroup species into one if inferring common ancestor of several ingroups #
          SP1="BFAL" ; SP2="LAAL" ; INGROUP="BFAL_LAAL"
          for scaff in `cat list_allscafs.bed | cut -f1`
          do
                ((j=j%nparallel)); ((j++==0)) && wait
                awk '{FS=OFS="\t"} FNR==NR {scaff[$1 SUBSEP $2]=$0; next}
                     {idx=$1 SUBSEP $2; if(idx in scaff) print scaff[$1 SUBSEP $2], $0}' \
                     01_depth/recal/tmp/${SP1}.ACGTpersite_${scaff}.tmp \
                     01_depth/recal/tmp/${SP2}.ACGTpersite_${scaff}.tmp \
                     | awk '{FS=OFS="\t"} {split($4,a,","); split($8,b,",")}
                            {print $1,$2, $3+$7, a[1]+b[1]","a[2]+b[2]","a[3]+b[3]","a[4]+b[4]}' \
                     > 01_depth/recal/tmp/${INGROUP}.ACGTpersite_${scaff}.tmp
          done

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

      list=$(for s in `cat list_allscafs.bed | cut -f1`; do (echo "01_depth/recal/tmp/${INGROUP}_outG.ACGTpersite_$s.tmp "); done)
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

      rm -f ${INGROUP}_outG.ACGTpersite.p_anc
      for file in `ls ${INGROUP}_outG.ACGTpersite.txtF-*.p_anc`
      do
            sed '1,8d' $file >> ${INGROUP}_outG.ACGTpersite.p_anc
      done

      paste <(cut -f1,2 ${INGROUP}_outG.ACGTpersite.txt) \
            ${INGROUP}_outG.ACGTpersite.p_anc \
            > ${INGROUP}_outG.ACGTpersite.p_anc.txt


# step 9: recode probabilities into allele calls #

      awk '{FS=OFS="\t"} 
           { 
            if($6 >= 0.05 || $7 >= 0.05 || $8 >= 0.05 || $9 >= 0.05) nA=1 ; 
            if($10 >= 0.05 || $11 >= 0.05 || $12 >= 0.05 || $13 >= 0.05) nC=2 ; 
            if($14 >= 0.05 || $15 >= 0.05 || $16 >= 0.05 || $17 >= 0.05) nG=4 ; 
            if($18 >= 0.05 || $19 >= 0.05 || $20 >= 0.05 || $21 >= 0.05) nT=16 ; 
            n = nA + nC + nG + nT ; 
            if(n == 1) nc="A" ; 
            else if(n == 2) nc="C" ; 
            else if(n == 4) nc="G" ; 
            else if(n == 16) nc="T" ; 
            else if(n == 3) nc="M" ; 
            else if(n == 5) nc="R" ; 
            else if(n == 17) nc="W" ; 
            else if(n == 6) nc="S" ; 
            else if(n == 18) nc="Y" ; 
            else if(n == 20) nc="K" ; 
            else if(n == 7) nc="V" ; 
            else if(n == 19) nc="H" ; 
            else if(n == 21) nc="D" ; 
            else if(n == 22) nc="B" ; 
            else nc="N" ; 
            print $1, $2, nc ; nA=0 ; nC=0 ; nG=0 ; nT=0
           }'
        ${INGROUP}_outG.ACGTpersite.p_anc.txt \
        > ${INGROUP}_outG.ACGTpersite.p_anc.txtF




###############################################################################
##   !!!!    TO BE CHECKED WITH ME BEFORE USING THE FOLLOWING PART    !!!!   ##
###############################################################################



# step 10: Fill in the missing sites in inferred ancestral genome (if present)

   ## Extract the missing regions and fill them ...
   awk 'OFS="\t"{print $1, $2-1, $2}' ${INGROUP}_outG.ACGTpersite.p_anc.txtF \
      | bedtools merge -i - \
      | bedtools complement -i - -g <( cut -f1,3 list_allscafs.bed )
      > ${INGROUP}_outG.ACGTpersite.p_anc.complement.bed
      
      # a. method 1 : with N's
            # write positions in list format and add N's in last column
            awk '{OFS="\t"; for (i = $2+1; i <= $3; ++i) print $1, i}' \
            ${INGROUP}_outG.ACGTpersite.p_anc.complement.bed \
            | sed "s/$/\tN/g" \
            > ${INGROUP}_outG.ACGTpersite.p_anc.complement_Ns.txt0

      # b. method 2 : with alleles from the closest selected outgroup species
            # extract DNA sequences of corresponding regions
            seqkit subseq -w0 --bed ${INGROUP}_outG.ACGTpersite.p_anc.complement.bed \
                 ${OUTG1}_genome.fasta \
                 > ${INGROUP}_outG.ACGTpersite.p_anc.complement.fasta
            # write positions in list format
            grep "^>" ${INGROUP}_outG.ACGTpersite.p_anc.complement.fasta \
              | sed 's/>//g' | sed 's/:.//g' \
              | awk -F'[_-]' '{OFS="\t"; for (i = $3; i <= $4; ++i) print $1"_"$2, i}' \
              > ${INGROUP}_outG.ACGTpersite.p_anc.complement.txt1 &
            # write nucleotides in list format
            grep -v "^>" ${INGROUP}_outG.ACGTpersite.p_anc.complement.fasta \
              | fold -w1 > ${INGROUP}_outG.ACGTpersite.p_anc.complement.txt2 &
            # combine list of positions and nucleotides next to each other (similar to cbind in R)
            paste ${INGROUP}_outG.ACGTpersite.p_anc.complement.txt1 \
               ${INGROUP}_outG.ACGTpersite.p_anc.complement.txt2 \
               > ${INGROUP}_outG.ACGTpersite.p_anc.complement.txt0

         
   ## Combine all positions (ancestral + missing) and sort by scaffold/chr and positions
      # ( change to "XXX.complement.txt0" or "XXX.complement_Ns.txt0" according to the method chosen )
      
   cat ${INGROUP}_outG.ACGTpersite.p_anc.txtF \
       ${INGROUP}_outG.ACGTpersite.p_anc.complement.txt0 \
       > ${INGROUP}_outG.ACGTpersite.p_anc.full.txtF

   sort -V ${INGROUP}_outG.ACGTpersite.p_anc.full.txtF \
           > ${INGROUP}_outG.ACGTpersite.p_anc.full.txtF0



# step 11: Convert ancestral genome to FASTA format

      # split by scaffold
      awk '{FS=OFS="\t"} { prevfile=ofile;
           ofile="05_uSFS/est-sfs/tmp/BFAL_LAAL_outG.noBAQ.ACGTpersite.p_anc.full_"$1".tmp";
           if(NR > 1 && ofile != prevfile) close(prevfile); print $0 >> ofile }' \
           ${INGROUP}_outG.ACGTpersite.p_anc.full.txtF0

      # convert list to fasta format for each scaffold/chr & append them sequentially to final ancestral FASTA file 
      rm -f ${INGROUP}_outG.ACGTpersite.ancGenome.fasta0
      for scaff in `cat list_allscafs.bed | cut -f1`
      do
              file="tmp/${INGROUP}_outG.ACGTpersite.p_anc.full_${scaff}.tmp"
              seq=$( cut -f3 $file | awk '{print}' ORS='' )

              echo -e ">${scaff}\n${seq}" >> ${INGROUP}_outG.ACGTpersite.ancGenome.fasta0
      done

      # rewrite the final ancestral FASTA file with interleaved sequences
      seqkit seq -w60 ${INGROUP}_outG.ACGTpersite.ancGenome.fasta0 \
                      > ${INGROUP}_outG.ACGTpersite.ancGenome.fasta

      # index the final ancestral FASTA file
      module load samtools/1.4
      samtools faidx ${INGROUP}_outG.ACGTpersite.ancGenome.fasta















