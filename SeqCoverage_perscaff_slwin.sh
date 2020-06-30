#!/bin/bash


for group in BFAL_m_allscafs_baQ BFAL_f_allscafs_baQ LAAL_m_allscafs_baQ LAAL_f_allscafs_baQ
do

  if [[ ${group} == "BFAL_f"* ]]; then        BAMLIST="BFAL_female.bamlist"
  elif [[ ${group} == "BFAL_m"* ]]; then      BAMLIST="BFAL_male.bamlist"
  elif [[ ${group} == "LAAL_f"* ]]; then      BAMLIST="LAAL_female.bamlist"
  elif [[ ${group} == "LAAL_m"* ]]; then      BAMLIST="LAAL_male.bamlist"
  fi

  REF="BFAL_genome_v2.fasta"
  REGION="list_autosomes_scafs_forANGSD.txt"
  
  
  ###  Get per-site read coverage using ANGSD  ###
  
  angsd -ref $REF \
        -bam $BAMLIST \
        -rf $REGION \
        -out $OUTFILE \
        -nThreads 10 \
        -remove_bads 1 -trim 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -only_proper_pair 0 \
        -baq 1 -C 50 \
        -doCounts 1 -doDepth 1 -dumpCounts 1


  ###  Summarize coverage obtained from ANGSD  ###

    tDIR="$HOME/tmp" && mkdir -p -m775 $tDIR
    OUT=$( echo $OUTFILE | sed 's/.pos$//' )
    PREFIX=$( basename $OUT )
  
    # Split file by scaffold
    sed '1d' $OUTFILE | awk -v tDIR=$tDIR -v PREFIX=$PREFIX '{ print > tDIR"/"PREFIX"_"$1".tmp" }'

    for file in 'ls ${tDIR}/${PREFIX}_*.tmp'
    do
  
        # Compute global depth per scaffold
        awk '{FS=OFS="\t"} {count[$1]++; sum[$1]+=$3} END {for(i in count) print i, count[i], sum[i], sum[i]/count[i]}' \
		    $file >> ${OUT}_globalDepth.txt

        # Compute sliding-window depth per scaffold
        SCAFF_SIZE=$( cut -f2 $file | tail -n1) ; WINSIZE=100000 ; STEPSIZE=50000
        for WIN in `seq 0 $STEPSIZE`
        do
            MIN="$WIN" ; MAX=$(( MIN + WINSIZE ))
            if [[ $MAX -gt $SCAFF_SIZE ]]; then MAX="$SCAFF_SIZE" ; fi
            MID=$(( MIN + (MAX-MIN)/2 ))
        
            if [[ $MIN -lt $SCAFF_SIZE ]]; then
                ((j=j%nparallel)); ((j++==0)) && wait
                
                awk -v MIN=$MIN -v MAX=$MAX -v MID=$MID \
				        '{FS=OFS="\t"} MIN <= $2 && $2 < MAX {count[$1]++; sum[$1]+=$3} END {for (i in count) print i, MIN, MAX, MID, count[i], sum[i], sum[i]/count[i]}' \
				        $file >> ${OUT}_${WINSIZE}kb_s${STEPSIZE}_winDepth.txt &
                
                sleep 1
            fi
            
        done

    done


done


