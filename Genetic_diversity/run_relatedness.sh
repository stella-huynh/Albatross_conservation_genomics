#!/bin/bash

SP=$1
BAMLIST=$2
REGION=$3
NIND=$4
OUTFILE=$5
NTHREADS=$6

angsd -ref BFAL_genome.fasta \
      -bam ${BAMLIST} \
      -rf ${REGION} \
      -out ${OUTFILE} \
      -nThreads ${NTHREADS} \
      -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 0 -minMapQ 20 -minQ 20 -trim 0 \
      -doMajorMinor 1 -skipTriallelic 1 \
      -doCounts 1 -doMaf 1 -SNP_pval 1e-6 \
      -doGeno -32 -doPost 1 \
      -GL 1 -doGlf 3

zcat ${OUTFILE}.mafs.gz | cut -f6 | sed 1d > ${OUTFILE}.mafs

ngsRelate -g ${OUTFILE}.glf.gz -f ${oDIR}/${OUTFILE}.mafs \
          -n ${NIND} -O ${oDIR}/${OUTFILE}.related


# with a little bit more filters (-l 0.05 eq.to -minMaf 0.05)?
ngsRelate -g ${OUTFILE}.glf.gz -f ${oDIR}/${OUTFILE}.mafs \
          -n ${NIND} -p {NTHREADS} -l 0.05 -O ${oDIR}/${OUTFILE}_minmaf.related
