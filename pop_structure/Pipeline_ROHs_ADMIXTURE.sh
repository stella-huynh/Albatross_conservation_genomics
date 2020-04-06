## PATHs ##

REF="/group/sbs_ssin/stella/albatross/00_genomes/01_WGS/01_dna/BFAL_genome_v2.fasta"
wDIR="/group/sbs_ssin/stella/albatross/05_reseq/initial_pipeline/06_ANGSDv2"
bDIR="/group/sbs_ssin/stella/albatross/05_reseq/initial_pipeline//05_GATKv2"
oDIR="${wDIR}/01_maf/ROH"


## SPECIES BAMLISTS ##

if [[ ${group} == "BFAL"* ]]; then         BAMLIST="${bDIR}/BFAL_16.bamlist"      ; NIND=16 ; mxDEPTH=160
elif [[ ${group} == "LAAL"* ]]; then         BAMLIST="${bDIR}/LAAL_15.bamlist"      ; NIND=15 ; mxDEPTH=150
fi

## REGION FILE ##

REGION="${wDIR}/list_autosomes_scafs_forANGSD.txt"

## FILTERINGS ##

SNPfilt=1
pval=1e-6 ;  mmaf=0.05 ;  pcut="0.95" ;  FILT="95pp_1e6snp_05maf"
OUT="${oDIR}/${group}_${FILT}"

ROHfilt=1
OUT2="${OUT}_w50h5m5"
hsnp=50 ; hden=50 ; hgap=1000 ; hwsnp=50 ; hwhet=5 ; hwmiss=5 ; hwthr=0.05 ; pwin=50 ; pstep=10 ; pr2=0.2 ; hkb=100



################################################################################################
## BASH RUN ## 

run_ROHs.sh  -ref ${REF} -i ${BAMLIST} -r ${REGION} -o ${OUT} -q ${OUT2} -t ${nt} \
             -m ${mmaf} -s ${pval} -g ${pcut} --min_depth ${NIND} --max_depth ${mxDEPTH} \
             -h1 ${hsnp} -h2 ${hden} -h3 ${hgap} -h4 ${hkb} \
             -hw1 ${hwsnp} -hw2 ${hwhet} -hw3 ${hwmiss} -hw4 ${hwthr} \
             -pw ${pwin} -ps ${pstep} -pr ${pr2} &
 
################################################################################################

# !! Modify the <.tfam> files for appropriate IIDs and SEX code !! #

################################################################################################
# ADMIXTURE ANALYSES pipeline #

## Convert output file (.tped format) to .ped format for run with ADMIXTURE software ##

for SP in BFAL LAAL
do

    # Prune SNPs for high-LD (r2=0.2) & very strong Hardy-Weinberg deviation (hwe_pval<1e-50)
    plink --tped 01_maf/ROH/${SP}_autosome_1e6snp_05maf.tped \
          --tfam 01_maf/ROH/${SP}_autosome_1e6snp_05maf.tfam \
          --allow-extra-chr --nonfounders \
          --out 01_maf/ROH/${SP}_autosome_1e6snp_05maf_ghwe_LD02 \
          --maf 0.05 --geno 0.1 --hwe 1e-50 --indep-pairwise 50 10 0.2 

    # Extract data for retained SNPs
    grep -wF -f ${SP}_autosome_1e6snp_05maf_ghwe_LD02.prune.in \
                ${SP}_autosome_1e6snp_05maf.tped \
                > ${SP}_autosome_1e6snp_05maf_ghwe_LD02.tped
    
    cp ${SP}_autosome_1e6snp_05maf.tfam ${SP}_autosome_1e6snp_05maf_ghwe_LD02.tfam
    
    # Convert into .PED format
    plink --tfile ${SP}_autosome_1e6snp_05maf_ghwe_LD02 \
          --recode 12 --allow-extra-chr \
          --out ${SP}_autosome_1e6snp_05maf_ghwe_LD02
    
    sed -i 's/scaffold_//g' ${SP}_autosome_1e6snp_05maf_ghwe_LD02.map
    
    # Run ADMIXTURE
    bash run_admixture.sh -i ${SP}_autosome_1e6snp_05maf_ghwe_LD01 -k 10 -p 2 -b 1000 -c 0.00001 &

done










