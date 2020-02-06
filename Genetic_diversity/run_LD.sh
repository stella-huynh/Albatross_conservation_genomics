#!/bin/bash

## This script was run for files with PREFIX:
## - BFAL_autosome / LAAL_autosome / BFAL_allscafs / LAAL_allscafs

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


## Use scripts provided with ngsLD to plot LD decay ##

fit_LDdecay.R  --ld_files ${oDIR}/ld_file_${SP}.list \
               --col 7 --ld r2 --fit_boot 10 \
               --out ${oDIR}/${SP}_autosome_filt_LDdecay.pdf &


# Use scripts provided with ngsLD to plot extract LD-unlinked SNPs ##

prune_graph.pl --in_file ${oDIR}/${SP}_autosome_filt_samp05.ld \
               --field_dist 7 --field_weight 11 \
               --max_kb_dist 1 --min_weight 0.5 \
               --out ${oDIR}/${SP}_autosome_filt_unlinked.id &

