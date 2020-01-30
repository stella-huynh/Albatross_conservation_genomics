#!/bin/bash

SP=$1
BAMLIST=$2
REGION=$3
OUTFILE=$4
NIND=$5
NTHREADS=$6

angsd -bam ${BAMLIST} \
      -rf ${REGION} \
      -dumpCounts 1 \
      -out 01_maf/${OUTFILE} \
      -nThreads ${NTHREADS} -nInd ${NIND} \
      -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 0 -minMapQ 20 \
      -doCounts 1 -minQ 20 -trim 0


zcat 01_maf/${OUTFILE}.pos.gz > 01_maf/${OUTFILE}.pos

bash posDepth_global.sh -i 01_maf/${OUTFILE}.pos \
                        -o 06_misc/${OUTFILE}.posU \
                        -d 06_misc/tmp -t ${NTHREADS}


Rscript Rplot_readCov_slwin2.R ${SP}

