#!/bin/bash

PREFIX=$1
BAMLIST=$2
REGION=$3
NIND=$4
NTHREADS=$5

# similar settings to "run_relatedness.sh"
angsd -ref BFAL_genome.fasta \
      -bam ${BAMLIST} \
      -rf ${REGION} \
      -out ${PREFIX} \
      -nThreads ${NTHREADS} \
      -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 0 -minMapQ 20 -minQ 20 -trim 0 \
      -doCounts -1 -doMaf -1 -SNP_pval 1e-6 \
      -doMajorMinor 1 -skipTriallelic 1 \
      -doGeno 32 -doPost 1 \
      -GL 1 -doGlf -3


zcat ${PREFIX}.glf.pos.gz > ${PREFIX}.glf.pos
zcat ${PREFIX}.glf.gz > ${PREFIX}.glf

N_SITES=`cat ${PREFIX}.glf.pos | wc -l`

ngsLD --pos ${PREFIX}.glf.pos \
      --geno ${PREFIX}.glf \
      --n_ind ${NIND} \
      --n_sites ${N_SITES} \
      --probs --log_scale \
      --out ${PREFIX}.ld	# Output files of ~500Gb !


ngsLD --pos ${PREFIX}.glf.pos \
      --geno ${PREFIX}.glf \
      --n_ind ${NIND} \
      --n_sites ${N_SITES} \
      --probs --log_scale --extend_out \
      --min_maf 0.05 --call_geno --call_thresh 0.10 --ignore_miss_data \
      --n_threads ${NTHREADS} \
      --out ${PREFIX}_filt.ld	# Add more filters...


