#!/bin/bash

rDIR=$1		# full path to reference genome BFAL
iDIR=$2		# full path to input files
oDIR=$3		# full path to store output files
BAMLIST=$4      # BFAL_LAAL_31.bamlist, BFAL_16.bamlist, LAAL_15.bamlist
NIND=$5         # 31, 16, 15  (for ngsCovar)
REGION=$6       # list_autosomes_scafs_forANGSD.txt / list_zchrom_scafs_forANGSD.txt
OUTFILE=$7
NCORES=$8


# Run ANGSD

angsd -ref ${rDIR}/BFAL_genome.fasta \
        -bam ${iDIR}/${BAMLIST} \
        -rf ${oDIR}/${REGION} \
        -out ${oDIR}/${SP} \
        -nThreads ${NCORES} \
        -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 0 -minMapQ 20 -minQ 20 -trim 0 \
        -GL 1 -doGeno 32 -doPost 1 \
        -doMaf 1 -SNP_pval 1e-6 \
        -doCounts 1 -doMajorMinor 1 -skipTriallelic 1

# Run ngsCovar

zcat ${oDIR}/${SP}.geno.gz > ${oDIR}/${SP}.geno
zcat ${oDIR}/${SP}.mafs.gz > ${oDIR}/${SP}.mafs


N_SITES=`cat ${oDIR}/${SP}.mafs | tail -n+2 | wc -l`

ngsCovar -probfile ${oDIR}/${SP}.geno \
         -outfile ${oDIR}/${SP}.covar \
         -nind ${NIND} \
         -minmaf 0.05 \
         -nsites ${N_SITES}


Rscript -e 'write.table(cbind(seq(1,34),
                              c("BFAL_10_20","BFAL_11_1","BFAL_12_BFAL674","BFAL_13_HEW052","BFAL_14_BFAL680","BFAL_15_18","BFAL_1_BFAL676","BFAL_220bp_700","BFAL_2_3","BFAL_3_BFAL698","BFAL_4_HEW236","BFAL_5_HEW077","BFAL_6_HEW182","BFAL_7_BFAL691","BFAL_8_HEW072","BFAL_9_4",
                                "LAAL_16_HEW116","LAAL_17_HEW162","LAAL_18_HEW169","LAAL_19_HEW173","LAAL_20_HEW176","LAAL_21_HEW179","LAAL_22_HEW188","LAAL_23_HEW129","LAAL_24_HEW133","LAAL_25_HEW143","LAAL_26_HEW158","LAAL_27_HEW192","LAAL_28_HEW204","LAAL_29_HEW218","LAAL_30_HEW231",
                                "STAL_220bp","WAAL_220bp","WAL_220bp"),
                              c("BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL","BFAL",
                                "LAAL","LAAL","LAAL","LAAL","LAAL","LAAL","LAAL","LAAL","LAAL","LAAL","LAAL","LAAL","LAAL","LAAL","LAAL",
                                "STAL","WAAL","WAL")),
                        row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="01_maf/all_34albatross_annotations.clst", quote=F)'

grep -e "BFAL" -e "LAAL" 01_maf/all_34albatross_annotations.clst > 01_maf/BFAL_LAAL_annotations.clst
grep "BFAL" 01_maf/all_34albatross_annotations.clst > 01_maf/BFAL_annotations.clst
grep "LAAL" 01_maf/all_34albatross_annotations.clst > 01_maf/LAAL_annotations.clst


Rscript plot_PCA_ed.R -i 01_maf/${SP}_autosomes.covar \
                      -c 1-2 -a 01_maf/${SP}_annotations.clst \
                      -o 01_maf/${SP}_autosomes_pca_c1-2.pdf



