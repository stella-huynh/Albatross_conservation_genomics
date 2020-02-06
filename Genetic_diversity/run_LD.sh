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
      -doGeno -32 -doPost 1 -geno_minDepth 1 \
      -GL 1 -doGlf 3


zcat ${PREFIX}.glf.pos.gz > ${PREFIX}.glf.pos
zcat ${PREFIX}.glf.gz > ${PREFIX}.glf

N_SITES=`cat ${PREFIX}.glf.pos | wc -l`

ngsLD --pos ${PREFIX}.glf.pos \
      --geno ${PREFIX}.glf \
      --probs --log_scale \
      --n_ind ${NIND} --n_sites ${N_SITES} \
      --n_threads ${NTHREADS} \
      --out ${PREFIX}.ld	# Output files of ~500Gb with default max dist 100kb!


ngsLD --pos ${PREFIX}.glf.pos \
      --geno ${PREFIX}.glf \
      --probs --log_scale \
      --n_ind ${NIND} --n_sites ${N_SITES} \
      --rnd_sample 0.05 --max-kb_dist 10000000 \
      --min_maf 0.05 --call_geno --call_thresh 0.10 --ignore_miss_data \
      --n_threads ${NTHREADS} \
      --out ${PREFIX}_filt_samp05.ld	# Add more filters...


